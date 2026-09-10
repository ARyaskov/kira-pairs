//! The statistics accumulator and its derived summaries.

use std::collections::HashMap;

use crate::chroms::{ChromDict, ChromSizes};
use crate::pairs::record::PairKey;
use crate::stats::lambertw::lambert_w0;

/// Strand combinations in pairtools' `dist_freq` order.
pub const DIRS: [&str; 4] = ["+-", "-+", "--", "++"];
/// Strand order used by pairtools' convergence summary.
pub const CONVERGENCE_STRANDS: [&str; 4] = ["++", "--", "-+", "+-"];
/// Thresholds (bp) for the `cis_Nkb+` counters.
pub const CIS_KB: [u64; 6] = [1, 2, 4, 10, 20, 40];
const CONVERGENCE_REL_DIFF: f64 = 0.05;

/// Logarithmic genomic-distance bins (pairtools `PairCounter._dist_bins`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DistBins {
    edges: Vec<u64>,
}

impl DistBins {
    /// Default: 8 bins per decade from 10^0 to 10^9.
    pub const DEFAULT_PER_DECADE: usize = 8;

    /// `unique([0] + round(10 ** arange(0, 9.001, 1/n)))`.
    pub fn new(n_per_decade: usize) -> Self {
        let n = n_per_decade.max(1);
        let step = 1.0 / n as f64;
        let stop = 9.0 + 0.001;
        let count = ((stop - 0.0) / step).ceil() as usize;
        let mut edges = vec![0u64];
        for i in 0..count {
            let x = i as f64 * step;
            if x >= stop {
                break;
            }
            let v = 10f64.powf(x).round_ties_even();
            edges.push(v as u64);
        }
        edges.sort_unstable();
        edges.dedup();
        Self { edges }
    }

    /// Bins from explicit left edges (e.g. parsed from a stats file).
    pub fn from_edges(mut edges: Vec<u64>) -> Self {
        edges.sort_unstable();
        edges.dedup();
        Self { edges }
    }

    /// Left edges.
    pub fn edges(&self) -> &[u64] {
        &self.edges
    }

    /// Number of bins.
    pub fn len(&self) -> usize {
        self.edges.len()
    }

    /// True when there are no bins.
    pub fn is_empty(&self) -> bool {
        self.edges.is_empty()
    }

    /// Index of the bin containing `dist` (`searchsorted(right) - 1`).
    #[inline]
    pub fn index(&self, dist: u64) -> usize {
        // partition_point gives the first edge > dist.
        self.edges.partition_point(|&e| e <= dist).saturating_sub(1)
    }
}

impl Default for DistBins {
    fn default() -> Self {
        Self::new(Self::DEFAULT_PER_DECADE)
    }
}

/// Small insertion-ordered counter for pair types.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PairTypeCounts {
    entries: Vec<(Vec<u8>, u64)>,
}

impl PairTypeCounts {
    #[inline]
    fn add(&mut self, pt: &[u8], n: u64) {
        for (k, v) in &mut self.entries {
            if k.as_slice() == pt {
                *v += n;
                return;
            }
        }
        self.entries.push((pt.to_vec(), n));
    }

    /// Sorted `(pair_type, count)` pairs.
    pub fn sorted(&self) -> Vec<(String, u64)> {
        let mut v: Vec<(String, u64)> = self
            .entries
            .iter()
            .map(|(k, n)| (String::from_utf8_lossy(k).into_owned(), *n))
            .collect();
        v.sort();
        v
    }
}

/// Incremental pair statistics.
#[derive(Debug, Clone)]
pub struct StatsAccumulator {
    bins: DistBins,
    unmapped_id: Option<u32>,
    /// Total records.
    pub total: u64,
    /// Both sides unmapped.
    pub total_unmapped: u64,
    /// Exactly one side mapped.
    pub total_single_sided_mapped: u64,
    /// Both sides mapped.
    pub total_mapped: u64,
    /// Mapped duplicates.
    pub total_dups: u64,
    /// Mapped non-duplicates.
    pub total_nodups: u64,
    /// Cis non-duplicates.
    pub cis: u64,
    /// Trans non-duplicates.
    pub trans: u64,
    /// `cis_1kb+ .. cis_40kb+`.
    pub cis_kb: [u64; 6],
    pair_types: PairTypeCounts,
    chrom_freq: HashMap<(u32, u32), u64>,
    dist_freq: [Vec<u64>; 4],
    chromsizes: Option<ChromSizes>,
}

impl StatsAccumulator {
    /// Accumulator with the given distance bins. `unmapped_id` is the
    /// dictionary id of `!` (may be `None` until seen).
    pub fn new(bins: DistBins, unmapped_id: Option<u32>) -> Self {
        let n = bins.len();
        Self {
            bins,
            unmapped_id,
            total: 0,
            total_unmapped: 0,
            total_single_sided_mapped: 0,
            total_mapped: 0,
            total_dups: 0,
            total_nodups: 0,
            cis: 0,
            trans: 0,
            cis_kb: [0; 6],
            pair_types: PairTypeCounts::default(),
            chrom_freq: HashMap::new(),
            dist_freq: [vec![0; n], vec![0; n], vec![0; n], vec![0; n]],
            chromsizes: None,
        }
    }

    /// Default bins, unmapped id resolved from the dictionary.
    pub fn with_dict(dict: &ChromDict) -> Self {
        Self::new(DistBins::default(), dict.get(crate::chroms::UNMAPPED_CHROM))
    }

    /// Distance bins in use.
    pub fn bins(&self) -> &DistBins {
        &self.bins
    }

    /// Set the unmapped chromosome id (if it became known later).
    pub fn set_unmapped_id(&mut self, id: Option<u32>) {
        self.unmapped_id = id;
    }

    /// Attach chromosome sizes (written as `chromsizes/*`).
    pub fn set_chromsizes(&mut self, cs: ChromSizes) {
        self.chromsizes = Some(cs);
    }

    /// Chromosome sizes, if attached.
    pub fn chromsizes(&self) -> Option<&ChromSizes> {
        self.chromsizes.as_ref()
    }

    #[inline]
    fn dir_index(s1: u8, s2: u8) -> usize {
        match (s1, s2) {
            (b'+', b'-') => 0,
            (b'-', b'+') => 1,
            (b'-', b'-') => 2,
            _ => 3,
        }
    }

    /// Observe a record. `pair_type` overrides the key's pair type (used
    /// when duplicates are marked `DD`); `is_dup` marks a mapped duplicate.
    #[inline]
    pub fn observe(&mut self, key: &PairKey, pair_type: Option<&[u8]>, is_dup: bool) {
        self.total += 1;
        let pt = pair_type.unwrap_or_else(|| key.pair_type_bytes());
        self.pair_types.add(pt, 1);
        let (u1, u2) = match self.unmapped_id {
            Some(u) => (key.chrom1 == u, key.chrom2 == u),
            None => (false, false),
        };
        if u1 && u2 {
            self.total_unmapped += 1;
        } else if !u1 && !u2 {
            self.total_mapped += 1;
            if is_dup {
                self.total_dups += 1;
            } else {
                self.total_nodups += 1;
                *self.chrom_freq.entry((key.chrom1, key.chrom2)).or_insert(0) += 1;
                if key.chrom1 == key.chrom2 {
                    self.cis += 1;
                    let dist = key.pos2.abs_diff(key.pos1);
                    let bin = self.bins.index(dist);
                    self.dist_freq[Self::dir_index(key.strand1, key.strand2)][bin] += 1;
                    for (i, kb) in CIS_KB.iter().enumerate() {
                        if dist >= kb * 1000 {
                            self.cis_kb[i] += 1;
                        }
                    }
                } else {
                    self.trans += 1;
                }
            }
        } else {
            self.total_single_sided_mapped += 1;
        }
    }

    /// Observe a record whose duplicate status is given by `pair_type == DD`
    /// (the `stats` command on a plain file).
    #[inline]
    pub fn observe_plain(&mut self, key: &PairKey) {
        self.observe(key, None, key.is_dd());
    }

    /// Add another accumulator's counts (bins must match).
    pub fn merge(&mut self, other: &StatsAccumulator) {
        debug_assert_eq!(self.bins, other.bins);
        self.total += other.total;
        self.total_unmapped += other.total_unmapped;
        self.total_single_sided_mapped += other.total_single_sided_mapped;
        self.total_mapped += other.total_mapped;
        self.total_dups += other.total_dups;
        self.total_nodups += other.total_nodups;
        self.cis += other.cis;
        self.trans += other.trans;
        for i in 0..6 {
            self.cis_kb[i] += other.cis_kb[i];
        }
        for (k, n) in &other.pair_types.entries {
            self.pair_types.add(k, *n);
        }
        for (k, n) in &other.chrom_freq {
            *self.chrom_freq.entry(*k).or_insert(0) += n;
        }
        for d in 0..4 {
            for (a, b) in self.dist_freq[d].iter_mut().zip(other.dist_freq[d].iter()) {
                *a += b;
            }
        }
        if self.chromsizes.is_none() {
            self.chromsizes = other.chromsizes.clone();
        }
    }

    /// Resolve chromosome ids to names and produce a snapshot suitable for
    /// formatting or merging with snapshots from other files.
    pub fn snapshot(&self, dict: &ChromDict) -> StatsSnapshot {
        let mut chrom_freq: Vec<((String, String), u64)> = self
            .chrom_freq
            .iter()
            .map(|((a, b), n)| {
                (
                    (
                        String::from_utf8_lossy(&dict.name(*a)).into_owned(),
                        String::from_utf8_lossy(&dict.name(*b)).into_owned(),
                    ),
                    *n,
                )
            })
            .collect();
        chrom_freq.sort();
        StatsSnapshot {
            bins: self.bins.clone(),
            total: self.total,
            total_unmapped: self.total_unmapped,
            total_single_sided_mapped: self.total_single_sided_mapped,
            total_mapped: self.total_mapped,
            total_dups: self.total_dups,
            total_nodups: self.total_nodups,
            cis: self.cis,
            trans: self.trans,
            cis_kb: self.cis_kb,
            pair_types: self.pair_types.sorted(),
            chrom_freq,
            dist_freq: self.dist_freq.clone(),
            chromsizes: self
                .chromsizes
                .as_ref()
                .map(|cs| cs.iter().map(|(n, s)| (n.to_string(), s)).collect())
                .unwrap_or_default(),
        }
    }
}

/// Name-resolved statistics, the unit of formatting and merging.
#[derive(Debug, Clone, PartialEq)]
pub struct StatsSnapshot {
    /// Distance bins.
    pub bins: DistBins,
    /// See [`StatsAccumulator`].
    pub total: u64,
    /// See [`StatsAccumulator`].
    pub total_unmapped: u64,
    /// See [`StatsAccumulator`].
    pub total_single_sided_mapped: u64,
    /// See [`StatsAccumulator`].
    pub total_mapped: u64,
    /// See [`StatsAccumulator`].
    pub total_dups: u64,
    /// See [`StatsAccumulator`].
    pub total_nodups: u64,
    /// See [`StatsAccumulator`].
    pub cis: u64,
    /// See [`StatsAccumulator`].
    pub trans: u64,
    /// See [`StatsAccumulator`].
    pub cis_kb: [u64; 6],
    /// Sorted pair type counts.
    pub pair_types: Vec<(String, u64)>,
    /// Sorted `((chrom1, chrom2), count)`.
    pub chrom_freq: Vec<((String, String), u64)>,
    /// Per-direction distance histograms (index order of [`DIRS`]).
    pub dist_freq: [Vec<u64>; 4],
    /// Chromosome sizes (name, length) in header order.
    pub chromsizes: Vec<(String, u64)>,
}

/// Derived summary values (pairtools `calculate_summaries`).
#[derive(Debug, Clone, PartialEq)]
pub struct Summary {
    /// `frac_cis`, `frac_cis_1kb+`, ... (7 values; `Value::Int(0)` when undefined).
    pub frac_cis: Vec<(String, crate::stats::Value)>,
    /// `frac_dups`.
    pub frac_dups: crate::stats::Value,
    /// `complexity_naive`.
    pub complexity_naive: crate::stats::Value,
    /// `dist_freq_convergence` sub-tree.
    pub convergence: Vec<(String, crate::stats::Value)>,
}

impl StatsSnapshot {
    /// Empty snapshot with the given bins.
    pub fn empty(bins: DistBins) -> Self {
        let n = bins.len();
        Self {
            bins,
            total: 0,
            total_unmapped: 0,
            total_single_sided_mapped: 0,
            total_mapped: 0,
            total_dups: 0,
            total_nodups: 0,
            cis: 0,
            trans: 0,
            cis_kb: [0; 6],
            pair_types: Vec::new(),
            chrom_freq: Vec::new(),
            dist_freq: [vec![0; n], vec![0; n], vec![0; n], vec![0; n]],
            chromsizes: Vec::new(),
        }
    }

    /// Sum with another snapshot (bins must match; chromsizes must match or
    /// be absent on one side).
    pub fn merge(&mut self, other: &StatsSnapshot) -> Result<(), String> {
        if self.bins != other.bins {
            return Err("distance bins differ between stats files".into());
        }
        self.total += other.total;
        self.total_unmapped += other.total_unmapped;
        self.total_single_sided_mapped += other.total_single_sided_mapped;
        self.total_mapped += other.total_mapped;
        self.total_dups += other.total_dups;
        self.total_nodups += other.total_nodups;
        self.cis += other.cis;
        self.trans += other.trans;
        for i in 0..6 {
            self.cis_kb[i] += other.cis_kb[i];
        }
        for (k, n) in &other.pair_types {
            match self.pair_types.iter_mut().find(|(a, _)| a == k) {
                Some((_, v)) => *v += n,
                None => self.pair_types.push((k.clone(), *n)),
            }
        }
        self.pair_types.sort();
        for (k, n) in &other.chrom_freq {
            match self.chrom_freq.iter_mut().find(|(a, _)| a == k) {
                Some((_, v)) => *v += n,
                None => self.chrom_freq.push((k.clone(), *n)),
            }
        }
        self.chrom_freq.sort();
        for d in 0..4 {
            for (a, b) in self.dist_freq[d].iter_mut().zip(other.dist_freq[d].iter()) {
                *a += b;
            }
        }
        if self.chromsizes.is_empty() {
            self.chromsizes = other.chromsizes.clone();
        } else if !other.chromsizes.is_empty() && self.chromsizes != other.chromsizes {
            return Err("cannot merge stats with different chromsizes".into());
        }
        Ok(())
    }

    /// Total per-direction cis counts.
    fn dist_by_dir(&self, dir: &str) -> &[u64] {
        let i = DIRS.iter().position(|d| *d == dir).unwrap_or(0);
        &self.dist_freq[i]
    }

    /// pairtools `calculate_summaries` + `find_dist_freq_convergence_distance`.
    pub fn summary(&self) -> Summary {
        use crate::stats::Value;
        let frac = |num: u64, den: u64| -> Value {
            if den > 0 {
                Value::Float(num as f64 / den as f64)
            } else {
                Value::Int(0)
            }
        };
        let mut frac_cis = vec![("frac_cis".to_string(), frac(self.cis, self.total_nodups))];
        for (i, kb) in CIS_KB.iter().enumerate() {
            frac_cis.push((
                format!("frac_cis_{kb}kb+"),
                frac(self.cis_kb[i], self.total_nodups),
            ));
        }
        let frac_dups = frac(self.total_dups, self.total_mapped);
        let complexity_naive = estimate_library_complexity(self.total_mapped, self.total_dups);

        // Convergence analysis.
        let nb = self.bins.len();
        let freqs: Vec<&[u64]> = CONVERGENCE_STRANDS
            .iter()
            .map(|d| self.dist_by_dir(d))
            .collect();
        let mut idx_max = [0usize; 4];
        for (si, f) in freqs.iter().enumerate() {
            for b in 0..nb {
                let avg = (0..4).map(|k| freqs[k][b] as f64).sum::<f64>() / 4.0;
                let dev = if avg == 0.0 {
                    0.0
                } else {
                    ((f[b] as f64 - avg) / avg).abs()
                };
                if dev > CONVERGENCE_REL_DIFF {
                    idx_max[si] = b;
                }
            }
        }
        let mut conv_idx = 0usize;
        let mut conv_strands = "??".to_string();
        let mut conv_dist = Value::Str("0".into());
        for (si, s) in CONVERGENCE_STRANDS.iter().enumerate() {
            if idx_max[si] > conv_idx {
                conv_idx = idx_max[si];
                conv_strands = s.to_string();
                conv_dist = if conv_idx + 1 < nb {
                    Value::Int(self.bins.edges()[conv_idx + 1])
                } else {
                    Value::Int(i64::MAX as u64)
                };
            }
        }
        let below: Vec<u64> = freqs.iter().map(|f| f[..=conv_idx].iter().sum()).collect();
        let above: Vec<u64> = freqs
            .iter()
            .map(|f| f[conv_idx + 1..].iter().sum())
            .collect();
        let below_all: u64 = below.iter().sum();
        let above_all: u64 = above.iter().sum();
        let mut conv: Vec<(String, Value)> = vec![
            ("convergence_dist".into(), conv_dist),
            (
                "strands_w_max_convergence_dist".into(),
                Value::Str(conv_strands),
            ),
            (
                "convergence_rel_diff_threshold".into(),
                Value::Float(CONVERGENCE_REL_DIFF),
            ),
        ];
        for (si, s) in CONVERGENCE_STRANDS.iter().enumerate() {
            conv.push((
                format!("n_cis_pairs_below_convergence_dist/{s}"),
                Value::Int(below[si]),
            ));
        }
        conv.push((
            "n_cis_pairs_below_convergence_dist_all_strands".into(),
            Value::Int(below_all),
        ));
        conv.push((
            "n_cis_pairs_above_convergence_dist_all_strands".into(),
            Value::Int(above_all),
        ));
        let norms = [
            ("cis", self.cis),
            ("total_mapped", self.total_mapped),
            ("total_nodups", self.total_nodups),
        ];
        let div = |a: u64, b: u64| -> f64 {
            if b == 0 {
                if a == 0 { f64::NAN } else { f64::INFINITY }
            } else {
                a as f64 / b as f64
            }
        };
        for (name, norm) in norms {
            let mut sum = 0.0;
            for (si, s) in CONVERGENCE_STRANDS.iter().enumerate() {
                let v = div(below[si], norm);
                sum += v;
                conv.push((
                    format!("frac_{name}_in_cis_below_convergence_dist/{s}"),
                    Value::Float(v),
                ));
            }
            conv.push((
                format!("frac_{name}_in_cis_below_convergence_dist_all_strands"),
                Value::Float(sum),
            ));
            conv.push((
                format!("frac_{name}_in_cis_above_convergence_dist_all_strands"),
                Value::Float(div(above_all, norm)),
            ));
        }
        Summary {
            frac_cis,
            frac_dups,
            complexity_naive,
            convergence: conv,
        }
    }
}

/// pairtools `estimate_library_complexity(nseq, ndup, 0)`.
pub fn estimate_library_complexity(nseq: u64, ndup: u64) -> crate::stats::Value {
    use crate::stats::Value;
    if nseq == 0 {
        return Value::Int(0);
    }
    let u = (nseq as f64 - ndup as f64) / nseq as f64;
    if u == 0.0 {
        return Value::Int(0);
    }
    if ndup == 0 {
        // scipy's lambertw evaluated exactly at the branch point -1/e yields
        // NaN in pairtools' formula; reproduce that output.
        return Value::Float(f64::NAN);
    }
    let x = -(-1.0 / u).exp() / u;
    let seq_to_complexity = lambert_w0(x) + 1.0 / u;
    Value::Float(nseq as f64 / seq_to_complexity)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::Value;

    #[test]
    fn default_bins_match_pairtools() {
        let b = DistBins::default();
        let e = b.edges();
        assert_eq!(&e[..8], &[0, 1, 2, 3, 4, 6, 7, 10]);
        assert_eq!(*e.last().unwrap(), 1_000_000_000);
        assert_eq!(b.index(0), 0);
        assert_eq!(b.index(1), 1);
        assert_eq!(b.index(5), 4);
        assert_eq!(b.index(2_000_000_000), e.len() - 1);
        let b4 = DistBins::new(4);
        assert_eq!(&b4.edges()[..6], &[0, 1, 2, 3, 6, 10]);
    }

    fn key(c1: u32, c2: u32, p1: u64, p2: u64, s: &[u8; 2], pt: &[u8; 2]) -> PairKey {
        let mut k = PairKey {
            seq: 0,
            pos1: p1,
            pos2: p2,
            chrom1: c1,
            chrom2: c2,
            pair_type: [0; 8],
            strand1: s[0],
            strand2: s[1],
            flags: 0,
            pair_type_len: 2,
        };
        k.pair_type[..2].copy_from_slice(pt);
        k
    }

    #[test]
    fn counts_like_pairtools() {
        let dict = ChromDict::with_names(["!", "chr1", "chr2"]);
        let mut s = StatsAccumulator::with_dict(&dict);
        s.observe_plain(&key(0, 0, 0, 0, b"--", b"NN"));
        s.observe_plain(&key(0, 1, 0, 5, b"-+", b"NU"));
        s.observe_plain(&key(1, 1, 100, 1600, b"+-", b"UU"));
        s.observe_plain(&key(1, 1, 100, 1600, b"+-", b"DD"));
        s.observe_plain(&key(1, 2, 100, 1600, b"++", b"UU"));
        s.observe_plain(&key(1, 1, 100, 100_100, b"++", b"UR"));
        assert_eq!(s.total, 6);
        assert_eq!(s.total_unmapped, 1);
        assert_eq!(s.total_single_sided_mapped, 1);
        assert_eq!(s.total_mapped, 4);
        assert_eq!(s.total_dups, 1);
        assert_eq!(s.total_nodups, 3);
        assert_eq!(s.cis, 2);
        assert_eq!(s.trans, 1);
        assert_eq!(s.cis_kb, [2, 1, 1, 1, 1, 1]);
        let snap = s.snapshot(&dict);
        assert_eq!(
            snap.pair_types,
            vec![
                ("DD".into(), 1),
                ("NN".into(), 1),
                ("NU".into(), 1),
                ("UR".into(), 1),
                ("UU".into(), 2)
            ]
        );
        assert_eq!(snap.chrom_freq.len(), 2);
        let sum = snap.summary();
        assert_eq!(
            sum.frac_cis[0],
            ("frac_cis".into(), Value::Float(2.0 / 3.0))
        );
        assert_eq!(sum.frac_dups, Value::Float(0.25));
        match sum.complexity_naive {
            Value::Float(c) => assert!(c > 4.0 && c < 20.0, "{c}"),
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn complexity_matches_scipy() {
        // pairtools 1.1.3: estimate_library_complexity(1000, 100) == 4660.793479934838
        match estimate_library_complexity(1000, 100) {
            Value::Float(c) => assert!((c - 4660.793479934838).abs() < 1e-6, "{c}"),
            other => panic!("{other:?}"),
        }
        // pairtools 1.1.3: estimate_library_complexity(4, 1) == 6.602185563963388
        match estimate_library_complexity(4, 1) {
            Value::Float(c) => assert!((c - 6.602185563963388).abs() < 1e-9, "{c}"),
            other => panic!("{other:?}"),
        }
        assert_eq!(estimate_library_complexity(0, 0), Value::Int(0));
        assert_eq!(estimate_library_complexity(5, 5), Value::Int(0));
        match estimate_library_complexity(10, 0) {
            Value::Float(c) => assert!(c.is_nan()),
            other => panic!("{other:?}"),
        }
    }
}
