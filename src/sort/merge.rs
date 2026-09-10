//! K-way merging of sorted sources with a loser tree.

use crate::error::Result;
use crate::pairs::record::PairKey;
use crate::sort::key::{ParsedChunk, SortEntry, SortKeyContext};
use crate::sort::run::RunReader;

/// A sorted sequence of records that can be consumed one at a time.
pub trait RunSource {
    /// True when exhausted.
    fn is_done(&self) -> bool;
    /// Current key (valid when not done).
    fn key(&self) -> &PairKey;
    /// Current line (valid when not done).
    fn line(&self) -> &[u8];
    /// Advance to the next record; false when exhausted.
    fn advance(&mut self) -> Result<bool>;
}

impl RunSource for RunReader {
    #[inline]
    fn is_done(&self) -> bool {
        RunReader::is_done(self)
    }
    #[inline]
    fn key(&self) -> &PairKey {
        RunReader::key(self)
    }
    #[inline]
    fn line(&self) -> &[u8] {
        RunReader::line(self)
    }
    #[inline]
    fn advance(&mut self) -> Result<bool> {
        RunReader::advance(self)
    }
}

/// An in-memory sorted run.
pub struct MemorySource {
    chunks: Vec<ParsedChunk>,
    entries: Vec<SortEntry>,
    pos: usize,
    started: bool,
}

impl MemorySource {
    /// Wrap sorted entries and the chunks they reference.
    pub fn new(chunks: Vec<ParsedChunk>, entries: Vec<SortEntry>) -> Self {
        Self {
            chunks,
            entries,
            pos: 0,
            started: false,
        }
    }

    /// Whether the first record has already been handed out.
    pub fn started(&self) -> bool {
        self.started
    }

    /// Record that the first record has been handed out.
    pub fn mark_started(&mut self) {
        self.started = true;
    }

    /// Number of records.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True when there are no records.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

impl RunSource for MemorySource {
    #[inline]
    fn is_done(&self) -> bool {
        self.pos >= self.entries.len()
    }
    #[inline]
    fn key(&self) -> &PairKey {
        &self.entries[self.pos].key
    }
    #[inline]
    fn line(&self) -> &[u8] {
        let e = &self.entries[self.pos];
        self.chunks[e.chunk as usize].line(e)
    }
    #[inline]
    fn advance(&mut self) -> Result<bool> {
        self.pos += 1;
        Ok(self.pos < self.entries.len())
    }
}

/// Either a run file or an in-memory run.
pub enum Source {
    /// Run file on disk.
    File(RunReader),
    /// In-memory run.
    Memory(MemorySource),
}

impl RunSource for Source {
    #[inline]
    fn is_done(&self) -> bool {
        match self {
            Self::File(r) => r.is_done(),
            Self::Memory(m) => m.is_done(),
        }
    }
    #[inline]
    fn key(&self) -> &PairKey {
        match self {
            Self::File(r) => r.key(),
            Self::Memory(m) => m.key(),
        }
    }
    #[inline]
    fn line(&self) -> &[u8] {
        match self {
            Self::File(r) => r.line(),
            Self::Memory(m) => m.line(),
        }
    }
    #[inline]
    fn advance(&mut self) -> Result<bool> {
        match self {
            Self::File(r) => r.advance(),
            Self::Memory(m) => m.advance(),
        }
    }
}

/// Loser-tree k-way merger.
pub struct Merger<S: RunSource> {
    sources: Vec<S>,
    tree: Vec<usize>,
    ctx: SortKeyContext,
    pending: Option<usize>,
    k: usize,
}

const EMPTY: usize = usize::MAX;

impl<S: RunSource> Merger<S> {
    /// Build a merger over sources (each positioned on its first record).
    pub fn new(sources: Vec<S>, ctx: SortKeyContext) -> Self {
        let k = sources.len();
        let mut m = Self {
            sources,
            tree: vec![EMPTY; k.max(1)],
            ctx,
            pending: None,
            k,
        };
        m.build();
        m
    }

    /// Number of sources.
    pub fn n_sources(&self) -> usize {
        self.k
    }

    #[inline]
    fn less(&self, a: usize, b: usize) -> bool {
        let sa = &self.sources[a];
        let sb = &self.sources[b];
        match (sa.is_done(), sb.is_done()) {
            (true, true) => a < b,
            (true, false) => false,
            (false, true) => true,
            (false, false) => self
                .ctx
                .cmp_full(sa.key(), sa.line(), sb.key(), sb.line())
                .is_lt(),
        }
    }

    fn build(&mut self) {
        if self.k == 0 {
            return;
        }
        if self.k == 1 {
            self.tree[0] = 0;
            return;
        }
        for leaf in 0..self.k {
            let mut w = leaf;
            let mut node = (leaf + self.k) / 2;
            while node >= 1 {
                let occupant = self.tree[node];
                if occupant == EMPTY {
                    self.tree[node] = w;
                    w = EMPTY;
                    break;
                }
                if self.less(occupant, w) {
                    self.tree[node] = w;
                    w = occupant;
                }
                node /= 2;
            }
            if w != EMPTY {
                self.tree[0] = w;
            }
        }
    }

    #[inline]
    fn replay(&mut self, leaf: usize) {
        if self.k == 1 {
            return;
        }
        let mut w = leaf;
        let mut node = (leaf + self.k) / 2;
        while node >= 1 {
            let occupant = self.tree[node];
            if self.less(occupant, w) {
                self.tree[node] = w;
                w = occupant;
            }
            node /= 2;
        }
        self.tree[0] = w;
    }

    /// Next record in merged order (a lending iterator).
    #[inline]
    pub fn next_record(&mut self) -> Result<Option<(&PairKey, &[u8])>> {
        if self.k == 0 {
            return Ok(None);
        }
        if let Some(p) = self.pending.take() {
            self.sources[p].advance()?;
            self.replay(p);
        }
        let w = self.tree[0];
        if self.sources[w].is_done() {
            return Ok(None);
        }
        self.pending = Some(w);
        let s = &self.sources[w];
        Ok(Some((s.key(), s.line())))
    }

    /// Consume the whole stream through a callback.
    pub fn for_each<F: FnMut(&PairKey, &[u8]) -> Result<()>>(&mut self, mut f: F) -> Result<()> {
        while let Some((k, l)) = self.next_record()? {
            f(k, l)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chroms::ChromDict;
    use std::sync::Arc;

    fn make_source(keys: &[(u32, u64, u64)]) -> MemorySource {
        let mut data = Vec::new();
        let mut entries = Vec::new();
        for (i, (c, p1, seq)) in keys.iter().enumerate() {
            let line = format!("r{seq}\tc{c}\t{p1}\tc{c}\t0\t+\t+\tUU");
            let off = data.len() as u32;
            data.extend_from_slice(line.as_bytes());
            data.push(b'\n');
            entries.push(SortEntry {
                key: PairKey {
                    seq: *seq,
                    pos1: *p1,
                    pos2: 0,
                    chrom1: *c,
                    chrom2: *c,
                    pair_type: *b"UU\0\0\0\0\0\0",
                    strand1: b'+',
                    strand2: b'+',
                    flags: 0,
                    pair_type_len: 2,
                },
                chunk: 0,
                off,
                len: line.len() as u32,
            });
            let _ = i;
        }
        MemorySource::new(
            vec![ParsedChunk {
                data,
                entries: Vec::new(),
            }],
            entries,
        )
    }

    #[test]
    fn merges_in_order_for_various_k() {
        let dict = ChromDict::with_names(["c0", "c1", "c2"]);
        for k in [1usize, 2, 3, 5, 8, 13] {
            let mut sources = Vec::new();
            let mut expected = Vec::new();
            for s in 0..k {
                let keys: Vec<(u32, u64, u64)> = (0..50u64)
                    .map(|i| ((i % 3) as u32, i * 7 % 23, (s as u64) * 1000 + i))
                    .collect();
                let mut sorted = keys.clone();
                sorted.sort_by_key(|(c, p, seq)| (*c, *p, *seq));
                expected.extend(sorted.iter().cloned());
                sources.push(make_source(&sorted));
            }
            expected.sort_by_key(|(c, p, seq)| (*c, *p, *seq));
            let ctx = SortKeyContext::new(Arc::new(dict.ranks()), Arc::new(Vec::new()));
            let mut m = Merger::new(sources, ctx);
            let mut got = Vec::new();
            while let Some((key, line)) = m.next_record().unwrap() {
                assert!(line.starts_with(format!("r{}\t", key.seq).as_bytes()));
                got.push((key.chrom1, key.pos1, key.seq));
            }
            assert_eq!(got, expected, "k={k}");
        }
    }

    #[test]
    fn empty_sources() {
        let dict = ChromDict::with_names(["c0"]);
        let ctx = SortKeyContext::new(Arc::new(dict.ranks()), Arc::new(Vec::new()));
        let mut m: Merger<MemorySource> = Merger::new(vec![], ctx.clone());
        assert!(m.next_record().unwrap().is_none());
        let mut m = Merger::new(
            vec![
                make_source(&[]),
                make_source(&[(0, 1, 1)]),
                make_source(&[]),
            ],
            ctx,
        );
        assert_eq!(m.next_record().unwrap().unwrap().0.seq, 1);
        assert!(m.next_record().unwrap().is_none());
    }
}
