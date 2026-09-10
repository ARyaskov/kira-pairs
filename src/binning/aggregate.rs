//! Hash aggregation with spill-to-disk and merge-reduce.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
use std::sync::Arc;

use crate::error::{KiraError, Result};
use crate::io::temp::TempManager;

const ENTRY_BYTES: u64 = 24;
// HashMap overhead is roughly 2x the payload.
const ENTRY_COST: u64 = ENTRY_BYTES * 3;

/// Streaming aggregator for one resolution.
pub struct Aggregator {
    resolution: u64,
    map: HashMap<(u64, u64), u64>,
    cap: usize,
    runs: Vec<PathBuf>,
    temp: Arc<TempManager>,
}

impl Aggregator {
    /// Aggregator with an in-memory budget in bytes.
    pub fn new(resolution: u64, memory_bytes: u64, temp: Arc<TempManager>) -> Self {
        let cap = ((memory_bytes / ENTRY_COST) as usize).max(1024);
        Self {
            resolution,
            map: HashMap::with_capacity(cap.min(1 << 20)),
            cap,
            runs: Vec::new(),
            temp,
        }
    }

    /// Resolution in bp.
    pub fn resolution(&self) -> u64 {
        self.resolution
    }

    /// Count one contact.
    #[inline]
    pub fn add(&mut self, b1: u64, b2: u64) -> Result<()> {
        *self.map.entry((b1, b2)).or_insert(0) += 1;
        if self.map.len() >= self.cap {
            self.spill()?;
        }
        Ok(())
    }

    fn sorted_entries(&mut self) -> Vec<(u64, u64, u64)> {
        let mut v: Vec<(u64, u64, u64)> = self.map.drain().map(|((a, b), n)| (a, b, n)).collect();
        v.sort_unstable();
        v
    }

    fn spill(&mut self) -> Result<()> {
        let entries = self.sorted_entries();
        let path = self.temp.next_path("bins");
        let f = File::create(&path).map_err(|e| KiraError::io(&path, e))?;
        let mut w = BufWriter::with_capacity(1 << 20, f);
        let mut buf = Vec::with_capacity(entries.len() * 24);
        for (a, b, n) in &entries {
            buf.extend_from_slice(&a.to_le_bytes());
            buf.extend_from_slice(&b.to_le_bytes());
            buf.extend_from_slice(&n.to_le_bytes());
        }
        w.write_all(&buf).map_err(|e| KiraError::io(&path, e))?;
        w.flush().map_err(|e| KiraError::io(&path, e))?;
        self.temp.add_bytes(buf.len() as u64);
        log::info!(
            "spilled {} bin entries at {} bp",
            entries.len(),
            self.resolution
        );
        self.runs.push(path);
        Ok(())
    }

    /// Finish and return a sorted, reduced table.
    pub fn finish(mut self) -> Result<BinTable> {
        let memory = self.sorted_entries();
        let mut sources: Vec<Source> = Vec::with_capacity(self.runs.len() + 1);
        for p in &self.runs {
            sources.push(Source::File(FileSource::open(p)?));
        }
        sources.push(Source::Memory(MemorySource {
            entries: memory,
            pos: 0,
        }));
        Ok(BinTable {
            sources,
            heads: Vec::new(),
            initialised: false,
            temp: Arc::clone(&self.temp),
            runs: std::mem::take(&mut self.runs),
        })
    }
}

struct FileSource {
    path: PathBuf,
    inp: BufReader<File>,
    cur: Option<(u64, u64, u64)>,
}

impl FileSource {
    fn open(path: &PathBuf) -> Result<Self> {
        let f = File::open(path).map_err(|e| KiraError::io(path, e))?;
        let mut s = Self {
            path: path.clone(),
            inp: BufReader::with_capacity(1 << 20, f),
            cur: None,
        };
        s.advance()?;
        Ok(s)
    }

    fn advance(&mut self) -> Result<()> {
        let mut buf = [0u8; 24];
        let mut n = 0;
        while n < 24 {
            match self.inp.read(&mut buf[n..]) {
                Ok(0) => break,
                Ok(k) => n += k,
                Err(e) if e.kind() == std::io::ErrorKind::Interrupted => {}
                Err(e) => return Err(KiraError::io(&self.path, e)),
            }
        }
        self.cur = match n {
            0 => None,
            24 => Some((
                u64::from_le_bytes(buf[..8].try_into().unwrap_or([0; 8])),
                u64::from_le_bytes(buf[8..16].try_into().unwrap_or([0; 8])),
                u64::from_le_bytes(buf[16..].try_into().unwrap_or([0; 8])),
            )),
            _ => {
                return Err(KiraError::RunFile {
                    path: self.path.clone(),
                    message: "truncated bin run".into(),
                });
            }
        };
        Ok(())
    }
}

struct MemorySource {
    entries: Vec<(u64, u64, u64)>,
    pos: usize,
}

enum Source {
    File(FileSource),
    Memory(MemorySource),
}

impl Source {
    fn head(&self) -> Option<(u64, u64, u64)> {
        match self {
            Source::File(f) => f.cur,
            Source::Memory(m) => m.entries.get(m.pos).copied(),
        }
    }

    fn advance(&mut self) -> Result<()> {
        match self {
            Source::File(f) => f.advance(),
            Source::Memory(m) => {
                m.pos += 1;
                Ok(())
            }
        }
    }
}

/// Sorted `(bin1, bin2, count)` rows, merge-reduced lazily.
pub struct BinTable {
    sources: Vec<Source>,
    heads: Vec<Option<(u64, u64, u64)>>,
    initialised: bool,
    temp: Arc<TempManager>,
    runs: Vec<PathBuf>,
}

impl BinTable {
    /// Next row, summing equal keys across sources.
    pub fn next_row(&mut self) -> Result<Option<(u64, u64, u64)>> {
        if !self.initialised {
            self.heads = self.sources.iter().map(Source::head).collect();
            self.initialised = true;
        }
        // Find the minimal key among heads.
        let mut min: Option<(u64, u64)> = None;
        for h in self.heads.iter().flatten() {
            let k = (h.0, h.1);
            if min.is_none_or(|m| k < m) {
                min = Some(k);
            }
        }
        let Some(key) = min else {
            return Ok(None);
        };
        let mut total = 0u64;
        for i in 0..self.sources.len() {
            if let Some(h) = self.heads[i]
                && (h.0, h.1) == key
            {
                total += h.2;
                self.sources[i].advance()?;
                self.heads[i] = self.sources[i].head();
            }
        }
        Ok(Some((key.0, key.1, total)))
    }

    /// Clone for tests: only valid before iteration starts and for
    /// memory-only tables.
    #[cfg(test)]
    pub(crate) fn clone_for_test(&self) -> BinTable {
        let sources = self
            .sources
            .iter()
            .map(|s| match s {
                Source::Memory(m) => Source::Memory(MemorySource {
                    entries: m.entries.clone(),
                    pos: m.pos,
                }),
                Source::File(f) => Source::File(FileSource::open(&f.path).unwrap()),
            })
            .collect();
        BinTable {
            sources,
            heads: Vec::new(),
            initialised: false,
            temp: Arc::clone(&self.temp),
            runs: Vec::new(),
        }
    }
}

impl Drop for BinTable {
    fn drop(&mut self) {
        for p in &self.runs {
            self.temp.remove(p);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spills_and_reduces() {
        let temp = Arc::new(TempManager::new(None, "kira-bin-test").unwrap());
        let mut agg = Aggregator::new(10, 1, Arc::clone(&temp));
        agg.cap = 100;
        for round in 0..5 {
            for i in 0..250u64 {
                agg.add(i % 37, i % 11 + round).unwrap();
            }
        }
        assert!(!agg.runs.is_empty());
        let mut t = agg.finish().unwrap();
        let mut total = 0;
        let mut prev = None;
        while let Some((a, b, n)) = t.next_row().unwrap() {
            if let Some(p) = prev {
                assert!(p < (a, b));
            }
            prev = Some((a, b));
            total += n;
        }
        assert_eq!(total, 1250);
    }
}
