//! Sort entries and the pairtools-compatible comparison function.

use std::cmp::Ordering;
use std::sync::Arc;

use crate::chroms::ChromRanks;
use crate::pairs::record::{PairKey, split_fields};

/// One record inside a [`ParsedChunk`]: its key plus the location of the
/// original line in the chunk arena.
#[derive(Debug, Clone, Copy)]
#[repr(C)]
pub struct SortEntry {
    /// Parsed key.
    pub key: PairKey,
    /// Index of the chunk holding the line (assigned when a run is built).
    pub chunk: u32,
    /// Byte offset of the line inside the chunk data.
    pub off: u32,
    /// Line length (without newline).
    pub len: u32,
}

/// A block of input lines with their parsed entries.
#[derive(Debug, Default)]
pub struct ParsedChunk {
    /// Line bytes (lines separated by `\n`; entries point into this).
    pub data: Vec<u8>,
    /// Parsed entries in input order.
    pub entries: Vec<SortEntry>,
}

impl ParsedChunk {
    /// Approximate resident size in bytes while the chunk is part of a run
    /// being sorted: the line data plus the entries, which exist twice (in
    /// the chunk and in the concatenated sorted vector).
    pub fn memory_size(&self) -> usize {
        self.data.capacity() + 2 * self.entries.capacity() * std::mem::size_of::<SortEntry>()
    }

    /// Line bytes of an entry that belongs to this chunk.
    #[inline]
    pub fn line(&self, e: &SortEntry) -> &[u8] {
        &self.data[e.off as usize..(e.off + e.len) as usize]
    }
}

/// Everything needed to compare two entries: chromosome ranks and, for the
/// rare ties on inline fields, access to the lines themselves.
#[derive(Debug, Clone)]
pub struct SortKeyContext {
    /// Lexicographic chromosome ranks.
    pub ranks: Arc<ChromRanks>,
    /// Extra columns (0-based indices) compared lexicographically after
    /// `pair_type` (pairtools `--extra-col`).
    pub extra_cols: Arc<Vec<usize>>,
    /// Index of the `pair_type` column (for truncated pair types).
    pub pair_type_col: Option<usize>,
}

impl SortKeyContext {
    /// Build a context.
    pub fn new(ranks: Arc<ChromRanks>, extra_cols: Arc<Vec<usize>>) -> Self {
        Self {
            ranks,
            extra_cols,
            pair_type_col: Some(7),
        }
    }

    /// Set the `pair_type` column index used for the slow path.
    #[must_use]
    pub fn with_pair_type_col(mut self, col: Option<usize>) -> Self {
        self.pair_type_col = col;
        self
    }

    /// Compare two keys on the inline fields only (chromosomes, positions,
    /// inline pair type). Returns `None` when the inline data cannot decide
    /// (truncated pair types or extra columns requested).
    #[inline]
    pub fn cmp_inline(&self, a: &PairKey, b: &PairKey) -> Option<Ordering> {
        let o = self
            .ranks
            .rank(a.chrom1)
            .cmp(&self.ranks.rank(b.chrom1))
            .then_with(|| self.ranks.rank(a.chrom2).cmp(&self.ranks.rank(b.chrom2)))
            .then_with(|| a.pos1.cmp(&b.pos1))
            .then_with(|| a.pos2.cmp(&b.pos2));
        if o != Ordering::Equal {
            return Some(o);
        }
        let pa = a.pair_type_bytes();
        let pb = b.pair_type_bytes();
        let o = pa.cmp(pb);
        if o != Ordering::Equal {
            // Prefix comparison decides unless one is a strict prefix of the
            // other and truncated.
            let prefix_tie = (a.pair_type_truncated() || b.pair_type_truncated())
                && (pa.starts_with(pb) || pb.starts_with(pa));
            if !prefix_tie {
                return Some(o);
            }
            return None;
        }
        if a.pair_type_truncated() || b.pair_type_truncated() || !self.extra_cols.is_empty() {
            return None;
        }
        Some(a.seq.cmp(&b.seq))
    }

    /// Full comparison with access to the lines for the slow path.
    #[inline]
    pub fn cmp_full(&self, a: &PairKey, la: &[u8], b: &PairKey, lb: &[u8]) -> Ordering {
        match self.cmp_inline(a, b) {
            Some(o) => o,
            None => self.cmp_slow(a, la, b, lb),
        }
    }

    fn cmp_slow(&self, a: &PairKey, la: &[u8], b: &PairKey, lb: &[u8]) -> Ordering {
        let mut ea = Vec::with_capacity(16);
        let mut eb = Vec::with_capacity(16);
        split_fields(la, &mut ea);
        split_fields(lb, &mut eb);
        let fa = |i: usize| field(la, &ea, i);
        let fb = |i: usize| field(lb, &eb, i);
        if (a.pair_type_truncated() || b.pair_type_truncated())
            && let Some(pt_col) = self.pair_type_col
        {
            let o = fa(pt_col).cmp(fb(pt_col));
            if o != Ordering::Equal {
                return o;
            }
        }
        for &c in self.extra_cols.iter() {
            let o = fa(c).cmp(fb(c));
            if o != Ordering::Equal {
                return o;
            }
        }
        a.seq.cmp(&b.seq)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chroms::ChromDict;

    fn key(c1: u32, c2: u32, p1: u64, p2: u64, pt: &[u8], seq: u64) -> PairKey {
        let mut k = PairKey {
            seq,
            pos1: p1,
            pos2: p2,
            chrom1: c1,
            chrom2: c2,
            pair_type: [0; 8],
            strand1: b'+',
            strand2: b'+',
            flags: 0,
            pair_type_len: pt.len().min(8) as u8,
        };
        k.pair_type[..pt.len().min(8)].copy_from_slice(&pt[..pt.len().min(8)]);
        if pt.len() > 8 {
            k.flags |= PairKey::PT_TRUNCATED;
            k.pair_type_len = pt.len() as u8;
        }
        k
    }

    #[test]
    fn lexicographic_chromosomes() {
        let d = ChromDict::with_names(["chr2", "chr10", "chr1"]);
        let ctx = SortKeyContext::new(Arc::new(d.ranks()), Arc::new(Vec::new()));
        let a = key(2, 2, 5, 5, b"UU", 0); // chr1
        let b = key(1, 1, 1, 1, b"UU", 1); // chr10
        let c = key(0, 0, 1, 1, b"UU", 2); // chr2
        assert_eq!(ctx.cmp_inline(&a, &b), Some(Ordering::Less));
        assert_eq!(ctx.cmp_inline(&b, &c), Some(Ordering::Less));
        // positions numeric, pair type lexicographic, seq for ties
        let x = key(0, 0, 9, 1, b"UU", 0);
        let y = key(0, 0, 10, 1, b"UU", 1);
        assert_eq!(ctx.cmp_inline(&x, &y), Some(Ordering::Less));
        let x = key(0, 0, 1, 1, b"UR", 5);
        let y = key(0, 0, 1, 1, b"UU", 1);
        assert_eq!(ctx.cmp_inline(&x, &y), Some(Ordering::Less));
        let x = key(0, 0, 1, 1, b"UU", 5);
        let y = key(0, 0, 1, 1, b"UU", 6);
        assert_eq!(ctx.cmp_inline(&x, &y), Some(Ordering::Less));
        assert_eq!(ctx.cmp_inline(&y, &x), Some(Ordering::Greater));
    }

    #[test]
    fn truncated_pair_types_use_lines() {
        let d = ChromDict::with_names(["chr1"]);
        let ctx = SortKeyContext::new(Arc::new(d.ranks()), Arc::new(Vec::new()));
        let a = key(0, 0, 1, 1, b"abcdefghi", 0);
        let b = key(0, 0, 1, 1, b"abcdefghj", 1);
        assert!(ctx.cmp_inline(&a, &b).is_none());
        let la = b"r\tchr1\t1\tchr1\t1\t+\t+\tabcdefghi";
        let lb = b"r\tchr1\t1\tchr1\t1\t+\t+\tabcdefghj";
        assert_eq!(ctx.cmp_full(&a, la, &b, lb), Ordering::Less);
        assert_eq!(ctx.cmp_full(&b, lb, &a, la), Ordering::Greater);
    }

    #[test]
    fn extra_columns() {
        let d = ChromDict::with_names(["chr1"]);
        let ctx = SortKeyContext::new(Arc::new(d.ranks()), Arc::new(vec![8]));
        let a = key(0, 0, 1, 1, b"UU", 0);
        let b = key(0, 0, 1, 1, b"UU", 1);
        let la = b"r\tchr1\t1\tchr1\t1\t+\t+\tUU\tz";
        let lb = b"r\tchr1\t1\tchr1\t1\t+\t+\tUU\ta";
        assert_eq!(ctx.cmp_full(&a, la, &b, lb), Ordering::Greater);
    }
}
