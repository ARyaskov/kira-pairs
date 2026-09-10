//! Streaming, pairtools-compatible duplicate detection on block-sorted,
//! upper-triangular pairs.
//!
//! Two records are duplicate candidates when they share `chrom1`, `chrom2`,
//! `strand1`, `strand2` (and any `--extra-col-pair` columns) and their
//! positions are within `max_mismatch` under the chosen metric:
//!
//! * `max`: `max(|pos1_a - pos1_b|, |pos2_a - pos2_b|) <= max_mismatch`
//! * `sum`: `|pos1_a - pos1_b| + |pos2_a - pos2_b| <= max_mismatch`
//!
//! The default clustering is transitive (connected components, the
//! pairtools `scipy`/`sklearn` backends); the first record of a cluster in
//! input order is kept and the rest are duplicates. Greedy clustering
//! (pairtools `cython` backend) marks a record as a duplicate when it is
//! within reach of an earlier *non-duplicate* record.
//!
//! The algorithm is a sweep over the sorted stream: only records whose
//! `pos1` lies within `max_mismatch` of the current position are indexed,
//! so the work per record is proportional to the local density.

pub mod cluster;
pub mod window;

use std::collections::{HashMap, HashSet, VecDeque};

use crate::chroms::ChromDict;
use crate::error::{KiraError, Location, Result};
use crate::pairs::record::{PairKey, split_fields};
use cluster::UnionFind;
use window::{Neighbor, Window};

/// Distance metric for duplicate detection.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Method {
    /// Maximum of the per-side distances.
    #[default]
    Max,
    /// Sum of the per-side distances.
    Sum,
}

/// How duplicate groups are formed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Clustering {
    /// Connected components of the "within reach" graph (pairtools scipy).
    #[default]
    Transitive,
    /// Compare only against earlier non-duplicates (pairtools cython).
    Greedy,
}

/// Dedup parameters.
#[derive(Debug, Clone)]
pub struct DedupConfig {
    /// Maximum mismatch in bp.
    pub max_mismatch: u64,
    /// Metric.
    pub method: Method,
    /// Clustering mode.
    pub clustering: Clustering,
    /// Record the parent read ID for duplicates.
    pub keep_parent_id: bool,
    /// Extra column pairs that must match (0-based column indices).
    pub extra_col_pairs: Vec<(usize, usize)>,
    /// Column index of `readID` (for parent ids).
    pub readid_col: usize,
    /// Input name for diagnostics.
    pub input_name: String,
}

impl Default for DedupConfig {
    fn default() -> Self {
        Self {
            max_mismatch: 3,
            method: Method::Max,
            clustering: Clustering::Transitive,
            keep_parent_id: false,
            extra_col_pairs: Vec::new(),
            readid_col: 0,
            input_name: "-".into(),
        }
    }
}

/// Classification of an emitted record.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Outcome {
    /// At least one side unmapped: not considered for deduplication.
    Unmapped,
    /// Kept (the representative of its cluster).
    Unique,
    /// A duplicate; carries the parent read ID when requested.
    Duplicate {
        /// Read ID of the kept record (empty unless `keep_parent_id`).
        parent: Vec<u8>,
    },
}

/// A record handed back by the deduper, in input order.
#[derive(Debug)]
pub struct Emitted<'a> {
    /// Key of the record.
    pub key: &'a PairKey,
    /// Original line bytes.
    pub line: &'a [u8],
    /// Classification.
    pub outcome: Outcome,
}

/// Counters collected by the deduper.
#[derive(Debug, Default, Clone)]
pub struct DedupMetrics {
    /// Records seen.
    pub records: u64,
    /// Records with an unmapped side.
    pub unmapped: u64,
    /// Duplicates found.
    pub duplicates: u64,
    /// Largest number of records held back awaiting cluster closure.
    pub peak_pending: u64,
    /// Largest positional window size.
    pub peak_window: u64,
}

#[derive(Debug)]
struct Pending {
    key: PairKey,
    off: usize,
    len: usize,
    final_outcome: Option<Outcome>,
}

/// Streaming deduper.
pub struct Deduper {
    cfg: DedupConfig,
    unmapped_id: Option<u32>,
    pending: VecDeque<Pending>,
    base: u64,
    arena: Vec<u8>,
    uf: UnionFind,
    window: Window,
    block: (u32, u32),
    block_id: u64,
    block_max_pos1: u64,
    seen_blocks: HashSet<(u32, u32)>,
    parent_names: HashMap<u64, (Vec<u8>, u64)>,
    metrics: DedupMetrics,
    ends: Vec<u32>,
    neighbors: Vec<Neighbor>,
}

impl Deduper {
    /// Create a deduper. `dict` resolves the unmapped chromosome id.
    pub fn new(cfg: DedupConfig, dict: &ChromDict) -> Self {
        let unmapped_id = dict.get(crate::chroms::UNMAPPED_CHROM);
        Self {
            window: Window::new(cfg.max_mismatch),
            cfg,
            unmapped_id,
            pending: VecDeque::new(),
            base: 0,
            arena: Vec::with_capacity(1 << 16),
            uf: UnionFind::new(),
            block: (u32::MAX, u32::MAX),
            block_id: 0,
            block_max_pos1: 0,
            seen_blocks: HashSet::new(),
            parent_names: HashMap::new(),
            metrics: DedupMetrics::default(),
            ends: Vec::with_capacity(32),
            neighbors: Vec::with_capacity(16),
        }
    }

    /// Update the unmapped id if the dictionary learned it after creation.
    pub fn refresh_unmapped(&mut self, dict: &ChromDict) {
        if self.unmapped_id.is_none() {
            self.unmapped_id = dict.get(crate::chroms::UNMAPPED_CHROM);
        }
    }

    /// Metrics so far.
    pub fn metrics(&self) -> &DedupMetrics {
        &self.metrics
    }

    #[inline]
    fn is_unmapped(&self, key: &PairKey) -> bool {
        match self.unmapped_id {
            Some(u) => key.chrom1 == u || key.chrom2 == u,
            None => false,
        }
    }

    /// Signature of the non-positional match fields (strands + extra cols).
    fn signature(&mut self, key: &PairKey, line: &[u8]) -> u64 {
        // FNV-1a over strands and extra column bytes.
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut mix = |b: u8| {
            h ^= u64::from(b);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        };
        mix(key.strand1);
        mix(0xff);
        mix(key.strand2);
        if !self.cfg.extra_col_pairs.is_empty() {
            split_fields(line, &mut self.ends);
            for &(a, b) in &self.cfg.extra_col_pairs {
                for c in [a, b] {
                    mix(0xfe);
                    for &byte in field(line, &self.ends, c) {
                        mix(byte);
                    }
                }
            }
        }
        h
    }

    /// Concatenated extra-column values of a line (both columns of every
    /// pair, NUL separated), or `None` when no extra columns are configured.
    fn extra_values(&mut self, line: &[u8]) -> Option<Box<[u8]>> {
        if self.cfg.extra_col_pairs.is_empty() {
            return None;
        }
        split_fields(line, &mut self.ends);
        let mut v = Vec::new();
        for &(a, b) in &self.cfg.extra_col_pairs {
            for c in [a, b] {
                v.extend_from_slice(field(line, &self.ends, c));
                v.push(0);
            }
        }
        Some(v.into_boxed_slice())
    }

    #[inline]
    fn within(&self, pos1: u64, pos2: u64, other_pos1: u64, other_pos2: u64) -> bool {
        let d1 = pos1.abs_diff(other_pos1);
        let d2 = pos2.abs_diff(other_pos2);
        match self.cfg.method {
            Method::Max => d1.max(d2) <= self.cfg.max_mismatch,
            Method::Sum => d1.saturating_add(d2) <= self.cfg.max_mismatch,
        }
    }

    fn pending_line(&self, idx: u64) -> &[u8] {
        let p = &self.pending[(idx - self.base) as usize];
        &self.arena[p.off..p.off + p.len]
    }

    fn readid_of(&mut self, line: &[u8]) -> Vec<u8> {
        split_fields(line, &mut self.ends);
        field(line, &self.ends, self.cfg.readid_col).to_vec()
    }

    /// Push one record. Records that became final are handed to `sink` in
    /// input order.
    pub fn push<F>(&mut self, key: &PairKey, line: &[u8], lineno: u64, sink: &mut F) -> Result<()>
    where
        F: FnMut(Emitted<'_>) -> Result<()>,
    {
        self.metrics.records += 1;
        let idx = self.base + self.pending.len() as u64;
        if self.is_unmapped(key) {
            self.metrics.unmapped += 1;
            self.append_pending(key, line, Some(Outcome::Unmapped));
            return self.drain(false, sink);
        }
        let block = (key.chrom1, key.chrom2);
        if block != self.block {
            if !self.seen_blocks.insert(block) {
                return Err(self.not_sorted(lineno, "chromosome pair block appears twice"));
            }
            // Close everything from the previous block.
            self.block = block;
            self.block_id += 1;
            self.block_max_pos1 = 0;
            self.window.clear();
        } else if key.pos1 < self.block_max_pos1 {
            return Err(self.not_sorted(lineno, "pos1 decreases within a chromosome pair block"));
        }
        self.block_max_pos1 = key.pos1;
        self.window.expire(key.pos1);

        let sig = self.signature(key, line);
        let extra = self.extra_values(line);
        self.uf.push(idx, key.pos1, self.block_id);

        // Find neighbours in the window.
        let mut first_parent: Option<(u64, usize)> = None;
        let greedy = self.cfg.clustering == Clustering::Greedy;
        let mut neighbors = std::mem::take(&mut self.neighbors);
        self.window.candidates(key.pos2, sig, &mut neighbors);
        for n in &neighbors {
            if greedy && n.is_dup {
                continue;
            }
            if !self.within(key.pos1, key.pos2, n.pos1, n.pos2) {
                continue;
            }
            if extra.is_some() && self.window.entry(n.slot).extra != extra {
                continue;
            }
            if greedy {
                if first_parent.is_none() {
                    first_parent = Some((n.idx, n.slot));
                }
            } else {
                self.uf.union(idx, n.idx);
            }
        }
        neighbors.clear();
        self.neighbors = neighbors;

        let final_outcome = if greedy {
            match first_parent {
                Some((_p, slot)) => {
                    let parent = self
                        .window
                        .entry(slot)
                        .readid
                        .as_ref()
                        .map(|r| r.to_vec())
                        .unwrap_or_default();
                    Some(Outcome::Duplicate { parent })
                }
                None => Some(Outcome::Unique),
            }
        } else {
            None
        };
        let is_dup = matches!(final_outcome, Some(Outcome::Duplicate { .. }));
        let readid = if greedy && self.cfg.keep_parent_id && !is_dup {
            Some(self.readid_of(line).into_boxed_slice())
        } else {
            None
        };
        self.window.insert(window::Entry {
            idx,
            pos1: key.pos1,
            pos2: key.pos2,
            sig,
            is_dup,
            extra,
            readid,
        });
        self.metrics.peak_window = self.metrics.peak_window.max(self.window.len() as u64);
        self.append_pending(key, line, final_outcome);
        self.drain(false, sink)
    }

    fn not_sorted(&self, lineno: u64, what: &str) -> KiraError {
        KiraError::NotSorted {
            message: format!("{what}; run `kira-pairs sort` first"),
            location: Location::file(&self.cfg.input_name).at_line(lineno),
        }
    }

    fn append_pending(&mut self, key: &PairKey, line: &[u8], final_outcome: Option<Outcome>) {
        let off = self.arena.len();
        self.arena.extend_from_slice(line);
        self.pending.push_back(Pending {
            key: *key,
            off,
            len: line.len(),
            final_outcome,
        });
        self.metrics.peak_pending = self.metrics.peak_pending.max(self.pending.len() as u64);
    }

    /// Emit all pending records whose status is final.
    fn drain<F>(&mut self, finishing: bool, sink: &mut F) -> Result<()>
    where
        F: FnMut(Emitted<'_>) -> Result<()>,
    {
        while let Some(head) = self.pending.front() {
            let idx = self.base;
            let outcome = match &head.final_outcome {
                Some(o) => o.clone(),
                None => {
                    // Transitive: final once the cluster can no longer grow.
                    let root = self.uf.find(idx);
                    let closed = finishing
                        || self.uf.block_of(root) < self.block_id
                        || self
                            .uf
                            .last_pos1(root)
                            .saturating_add(self.cfg.max_mismatch)
                            < self.block_max_pos1;
                    if !closed {
                        break;
                    }
                    if root == idx {
                        if self.cfg.keep_parent_id && self.uf.size(root) > 1 {
                            let line = self.pending_line(idx).to_vec();
                            let name = self.readid_of(&line);
                            self.parent_names
                                .insert(root, (name, self.uf.size(root) - 1));
                        }
                        Outcome::Unique
                    } else {
                        let parent = if self.cfg.keep_parent_id {
                            match self.parent_names.get_mut(&root) {
                                Some((name, remaining)) => {
                                    let n = name.clone();
                                    *remaining -= 1;
                                    if *remaining == 0 {
                                        self.parent_names.remove(&root);
                                    }
                                    n
                                }
                                None => Vec::new(),
                            }
                        } else {
                            Vec::new()
                        };
                        Outcome::Duplicate { parent }
                    }
                }
            };
            if matches!(outcome, Outcome::Duplicate { .. }) {
                self.metrics.duplicates += 1;
            }
            let p = self.pending.pop_front().unwrap_or_else(|| unreachable!());
            let line = &self.arena[p.off..p.off + p.len];
            sink(Emitted {
                key: &p.key,
                line,
                outcome,
            })?;
            self.base += 1;
        }
        // Reclaim memory when nothing is pending.
        if self.pending.is_empty() {
            self.arena.clear();
            self.uf.truncate_before(self.base);
        } else if self.arena.len() > (64 << 20) {
            // Long-lived pending queue: compact the arena.
            let first = self.pending.front().map(|p| p.off).unwrap_or(0);
            if first > 0 {
                self.arena.drain(..first);
                for p in &mut self.pending {
                    p.off -= first;
                }
            }
            self.uf.truncate_before(self.base);
        }
        Ok(())
    }

    /// Flush all remaining records.
    pub fn finish<F>(&mut self, sink: &mut F) -> Result<()>
    where
        F: FnMut(Emitted<'_>) -> Result<()>,
    {
        self.drain(true, sink)
    }
}

#[inline]
fn field<'a>(line: &'a [u8], ends: &[u32], i: usize) -> &'a [u8] {
    match ends.get(i) {
        Some(&e) => {
            let s = if i == 0 { 0 } else { ends[i - 1] as usize + 1 };
            &line[s..e as usize]
        }
        None => b"",
    }
}

/// Rewrite the pair type column of `line` to `DD` (pairtools `--mark-dups`).
pub fn mark_dd(line: &[u8], pair_type_col: Option<usize>, out: &mut Vec<u8>) {
    out.clear();
    let Some(col) = pair_type_col else {
        out.extend_from_slice(line);
        return;
    };
    let mut ends = Vec::with_capacity(16);
    split_fields(line, &mut ends);
    if col >= ends.len() {
        out.extend_from_slice(line);
        return;
    }
    let start = if col == 0 {
        0
    } else {
        ends[col - 1] as usize + 1
    };
    let end = ends[col] as usize;
    out.extend_from_slice(&line[..start]);
    out.extend_from_slice(b"DD");
    out.extend_from_slice(&line[end..]);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pairs::columns::ColumnMap;
    use crate::pairs::record::parse_line;

    fn run(lines: &[&str], cfg: DedupConfig) -> Vec<(String, Outcome)> {
        let dict = ChromDict::new();
        dict.intern(b"!");
        let cols = ColumnMap::standard();
        let mut ends = Vec::new();
        let mut d = Deduper::new(cfg, &dict);
        let mut out = Vec::new();
        let mut sink = |e: Emitted<'_>| {
            out.push((String::from_utf8_lossy(e.line).into_owned(), e.outcome));
            Ok(())
        };
        for (i, l) in lines.iter().enumerate() {
            let key = parse_line(
                l.as_bytes(),
                &cols,
                &dict,
                i as u64,
                &mut ends,
                Location::default,
            )
            .unwrap();
            d.push(&key, l.as_bytes(), i as u64 + 1, &mut sink).unwrap();
        }
        d.finish(&mut sink).unwrap();
        out
    }

    fn line(id: &str, c1: &str, p1: u64, c2: &str, p2: u64, s: &str) -> String {
        format!("{id}\t{c1}\t{p1}\t{c2}\t{p2}\t{}\t{}\tUU", &s[..1], &s[1..])
    }

    fn outcomes(v: &[(String, Outcome)]) -> Vec<&'static str> {
        v.iter()
            .map(|(_, o)| match o {
                Outcome::Unmapped => "X",
                Outcome::Unique => "U",
                Outcome::Duplicate { .. } => "D",
            })
            .collect()
    }

    #[test]
    fn exact_and_near_duplicates() {
        let lines = [
            line("a", "chr1", 100, "chr1", 200, "++"),
            line("b", "chr1", 100, "chr1", 200, "++"),
            line("e", "chr1", 100, "chr1", 200, "+-"),
            line("c", "chr1", 103, "chr1", 203, "++"),
            line("d", "chr1", 104, "chr1", 200, "++"),
            line("f", "chr1", 500, "chr2", 200, "++"),
        ];
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        let out = run(&refs, DedupConfig::default());
        // a kept, b exact dup, e differs by strand, c within 3 (max), d within 3 of c (transitive), f other block.
        assert_eq!(outcomes(&out), vec!["U", "D", "U", "D", "D", "U"]);
        let out = run(
            &refs,
            DedupConfig {
                method: Method::Sum,
                ..Default::default()
            },
        );
        // sum: c is 3+3=6 away from a and b -> not dup; d is 1+3=4 from c, 4 from a -> not dup.
        assert_eq!(outcomes(&out), vec!["U", "D", "U", "U", "U", "U"]);
    }

    #[test]
    fn transitive_chain_versus_greedy() {
        let lines = [
            line("a", "chr1", 100, "chr1", 100, "++"),
            line("b", "chr1", 103, "chr1", 100, "++"),
            line("c", "chr1", 106, "chr1", 100, "++"),
        ];
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        let out = run(
            &refs,
            DedupConfig {
                keep_parent_id: true,
                ..Default::default()
            },
        );
        assert_eq!(outcomes(&out), vec!["U", "D", "D"]);
        assert!(matches!(&out[2].1, Outcome::Duplicate { parent } if parent == b"a"));
        let out = run(
            &refs,
            DedupConfig {
                clustering: Clustering::Greedy,
                keep_parent_id: true,
                ..Default::default()
            },
        );
        assert_eq!(outcomes(&out), vec!["U", "D", "U"]);
        assert!(matches!(&out[1].1, Outcome::Duplicate { parent } if parent == b"a"));
    }

    #[test]
    fn retroactive_merge_keeps_input_order() {
        // A and B are not neighbours, C joins both: B becomes a duplicate.
        let lines = [
            line("a", "chr1", 100, "chr1", 100, "++"),
            line("b", "chr1", 100, "chr1", 106, "++"),
            line("c", "chr1", 101, "chr1", 103, "++"),
            line("d", "chr1", 900, "chr1", 100, "++"),
        ];
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        let out = run(
            &refs,
            DedupConfig {
                keep_parent_id: true,
                ..Default::default()
            },
        );
        assert_eq!(outcomes(&out), vec!["U", "D", "D", "U"]);
        let ids: Vec<&str> = out
            .iter()
            .map(|(l, _)| l.split('\t').next().unwrap())
            .collect();
        assert_eq!(ids, vec!["a", "b", "c", "d"]);
        assert!(matches!(&out[1].1, Outcome::Duplicate { parent } if parent == b"a"));
    }

    #[test]
    fn unmapped_and_unsorted() {
        let lines = [
            "u\t!\t0\tchr1\t5\t-\t+\tNU".to_string(),
            line("a", "chr1", 100, "chr1", 100, "++"),
            line("b", "chr1", 90, "chr1", 100, "++"),
        ];
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        let dict = ChromDict::new();
        dict.intern(b"!");
        let cols = ColumnMap::standard();
        let mut ends = Vec::new();
        let mut d = Deduper::new(DedupConfig::default(), &dict);
        let mut sink = |_e: Emitted<'_>| Ok(());
        let mut err = None;
        for (i, l) in refs.iter().enumerate() {
            let key = parse_line(
                l.as_bytes(),
                &cols,
                &dict,
                i as u64,
                &mut ends,
                Location::default,
            )
            .unwrap();
            if let Err(e) = d.push(&key, l.as_bytes(), i as u64 + 1, &mut sink) {
                err = Some(e);
                break;
            }
        }
        let e = err.expect("unsorted input must be rejected");
        assert!(e.to_string().contains("not sorted"), "{e}");
        assert!(e.to_string().contains("line 3"), "{e}");
    }

    #[test]
    fn extra_columns_must_match() {
        let lines = [
            "a\tchr1\t100\tchr1\t200\t+\t+\tUU\t1\t1",
            "b\tchr1\t100\tchr1\t200\t+\t+\tUU\t1\t1",
            "c\tchr1\t100\tchr1\t200\t+\t+\tUU\t0\t1",
        ];
        let out = run(
            &lines,
            DedupConfig {
                extra_col_pairs: vec![(8, 9)],
                ..Default::default()
            },
        );
        assert_eq!(outcomes(&out), vec!["U", "D", "U"]);
    }

    #[test]
    fn marks_dd() {
        let mut out = Vec::new();
        mark_dd(b"r\tchr1\t1\tchr1\t2\t+\t-\tUU\tx", Some(7), &mut out);
        assert_eq!(out, b"r\tchr1\t1\tchr1\t2\t+\t-\tDD\tx");
        mark_dd(b"r\tchr1\t1\tchr1\t2\t+\t-", Some(7), &mut out);
        assert_eq!(out, b"r\tchr1\t1\tchr1\t2\t+\t-");
    }

    #[test]
    fn large_window_uses_index() {
        // 5000 records at the same pos1 with distinct pos2: none are dups,
        // except exact repeats.
        let mut lines = Vec::new();
        for i in 0..5000u64 {
            lines.push(line(
                &format!("r{i}"),
                "chr1",
                1000,
                "chr1",
                10 + i * 10,
                "++",
            ));
            if i % 100 == 0 {
                lines.push(line(
                    &format!("d{i}"),
                    "chr1",
                    1000,
                    "chr1",
                    10 + i * 10 + 2,
                    "++",
                ));
            }
        }
        let refs: Vec<&str> = lines.iter().map(String::as_str).collect();
        let out = run(&refs, DedupConfig::default());
        let n_dups = out
            .iter()
            .filter(|(_, o)| matches!(o, Outcome::Duplicate { .. }))
            .count();
        assert_eq!(n_dups, 50);
    }
}
