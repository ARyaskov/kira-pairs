//! Upper-triangular flipping with pairtools `flip` semantics.
//!
//! Given a chromosome order (`!` first, then the chromosomes file order),
//! a pair is in the correct orientation when
//! `(order(chrom1), pos1) <= (order(chrom2), pos2)`. Pairs with exactly one
//! annotated chromosome put the annotated side first; pairs with no
//! annotated chromosome compare the names lexicographically (`chrom1 <
//! chrom2` is correct, so equal unannotated names are flipped, exactly as
//! pairtools does). Flipping swaps every `xxx1`/`xxx2` column pair and
//! reverses the two characters of `pair_type`.

use std::sync::Arc;

use crate::chroms::{ChromDict, ChromOrder};
use crate::pairs::columns::ColumnMap;
use crate::pairs::record::{PairKey, PairRecordRef};

/// Decides and applies flips.
pub struct Flipper {
    order: ChromOrder,
    dict: Arc<ChromDict>,
    enum_by_id: Vec<Option<u32>>,
    /// Column permutation applied when flipping (`out[i] = in[perm[i]]`).
    perm: Vec<usize>,
    pair_type_col: Option<usize>,
    unannotated_seen: bool,
}

impl Flipper {
    /// Build a flipper for the given chromosome order and columns.
    pub fn new(order: ChromOrder, cols: &ColumnMap, dict: Arc<ChromDict>) -> Self {
        let n = cols.len();
        let mut perm: Vec<usize> = (0..n).collect();
        for (a, b) in cols.side_pairs() {
            perm[a] = b;
            perm[b] = a;
        }
        Self {
            order,
            dict,
            enum_by_id: Vec::new(),
            perm,
            pair_type_col: cols.pair_type,
            unannotated_seen: false,
        }
    }

    /// Enumeration value of a chromosome id (`None` when not annotated).
    #[inline]
    pub fn enum_of(&mut self, id: u32) -> Option<u32> {
        let i = id as usize;
        if i >= self.enum_by_id.len() {
            let n = self.dict.len().max(i + 1);
            for j in self.enum_by_id.len()..n {
                let name = self.dict.name(j as u32);
                self.enum_by_id.push(self.order.get(&name));
            }
        }
        self.enum_by_id[i]
    }

    /// True when the pair must be flipped (pairtools `flip` rule).
    pub fn needs_flip(&mut self, key: &PairKey) -> bool {
        let e1 = self.enum_of(key.chrom1);
        let e2 = self.enum_of(key.chrom2);
        let correct = match (e1, e2) {
            (Some(a), Some(b)) => (a, key.pos1) <= (b, key.pos2),
            (Some(_), None) => true,
            (None, Some(_)) => false,
            (None, None) => {
                if !self.unannotated_seen {
                    log::warn!("unannotated chromosomes in the pairs file");
                    self.unannotated_seen = true;
                }
                self.dict.name(key.chrom1) < self.dict.name(key.chrom2)
            }
        };
        if (e1.is_none() || e2.is_none()) && !self.unannotated_seen {
            log::warn!("unannotated chromosomes in the pairs file");
            self.unannotated_seen = true;
        }
        !correct
    }

    /// Key with sides swapped.
    pub fn flip_key(key: &PairKey) -> PairKey {
        let mut k = *key;
        std::mem::swap(&mut k.chrom1, &mut k.chrom2);
        std::mem::swap(&mut k.pos1, &mut k.pos2);
        std::mem::swap(&mut k.strand1, &mut k.strand2);
        if k.pair_type_len >= 2 {
            k.pair_type.swap(0, 1);
        }
        k
    }

    /// Write the flipped line into `out`.
    pub fn flip_line(&self, rec: &PairRecordRef<'_>, out: &mut Vec<u8>) {
        out.clear();
        let n = rec.n_fields();
        for i in 0..n {
            if i > 0 {
                out.push(b'\t');
            }
            let src = self.perm.get(i).copied().unwrap_or(i);
            let src = if src < n { src } else { i };
            let f = rec.field(src).unwrap_or(b"");
            if Some(i) == self.pair_type_col && f.len() >= 2 {
                out.push(f[1]);
                out.push(f[0]);
                out.extend_from_slice(&f[2..]);
            } else {
                out.extend_from_slice(f);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chroms::ChromSizes;
    use crate::pairs::record::{parse_line, split_fields};

    fn setup() -> (Flipper, ColumnMap, Arc<ChromDict>) {
        let cs = ChromSizes::parse("chr1\t100\nchr2\t100\n", "t").unwrap();
        let names: Vec<String> = [
            "readID",
            "chrom1",
            "pos1",
            "chrom2",
            "pos2",
            "strand1",
            "strand2",
            "pair_type",
            "mapq1",
            "mapq2",
            "extra",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let cols = ColumnMap::from_names(names).unwrap();
        let dict = Arc::new(ChromDict::new());
        (
            Flipper::new(ChromOrder::from_chromsizes(&cs), &cols, Arc::clone(&dict)),
            cols,
            dict,
        )
    }

    fn check(line: &str, expect_flip: bool, expect: &str) {
        let (mut f, cols, dict) = setup();
        let mut ends = Vec::new();
        let key = parse_line(
            line.as_bytes(),
            &cols,
            &dict,
            0,
            &mut ends,
            Default::default,
        )
        .unwrap();
        assert_eq!(f.needs_flip(&key), expect_flip, "{line}");
        if expect_flip {
            split_fields(line.as_bytes(), &mut ends);
            let rec = PairRecordRef::new(line.as_bytes(), &ends);
            let mut out = Vec::new();
            f.flip_line(&rec, &mut out);
            assert_eq!(String::from_utf8(out).unwrap(), expect);
        }
    }

    #[test]
    fn flip_rules() {
        check("r\tchr1\t10\tchr1\t20\t+\t-\tUU\t30\t40\tx", false, "");
        check(
            "r\tchr1\t20\tchr1\t10\t+\t-\tUR\t30\t40\tx",
            true,
            "r\tchr1\t10\tchr1\t20\t-\t+\tRU\t40\t30\tx",
        );
        check(
            "r\tchr2\t1\tchr1\t50\t+\t-\tUU\t30\t40\tx",
            true,
            "r\tchr1\t50\tchr2\t1\t-\t+\tUU\t40\t30\tx",
        );
        check("r\t!\t0\tchr1\t50\t-\t+\tNU\t0\t40\tx", false, "");
        check(
            "r\tchr1\t50\t!\t0\t+\t-\tUN\t40\t0\tx",
            true,
            "r\t!\t0\tchr1\t50\t-\t+\tNU\t0\t40\tx",
        );
        // annotated side first
        check(
            "r\tchrZ\t5\tchr1\t50\t+\t-\tUU\t1\t2\tx",
            true,
            "r\tchr1\t50\tchrZ\t5\t-\t+\tUU\t2\t1\tx",
        );
        // both unannotated: lexicographic; equal names flip (pairtools quirk)
        check("r\tchrA\t5\tchrB\t50\t+\t-\tUU\t1\t2\tx", false, "");
        check(
            "r\tchrB\t5\tchrA\t50\t+\t-\tUU\t1\t2\tx",
            true,
            "r\tchrA\t50\tchrB\t5\t-\t+\tUU\t2\t1\tx",
        );
        check(
            "r\tchrA\t5\tchrA\t50\t+\t-\tUU\t1\t2\tx",
            true,
            "r\tchrA\t50\tchrA\t5\t-\t+\tUU\t2\t1\tx",
        );
    }

    #[test]
    fn flip_key_is_involution() {
        let (_, cols, dict) = setup();
        let mut ends = Vec::new();
        let key = parse_line(
            b"r\tchr2\t1\tchr1\t50\t+\t-\tUR\t30\t40\tx",
            &cols,
            &dict,
            0,
            &mut ends,
            Default::default,
        )
        .unwrap();
        let f = Flipper::flip_key(&key);
        assert_eq!(f.chrom1, key.chrom2);
        assert_eq!(f.pair_type_bytes(), b"RU");
        assert_eq!(Flipper::flip_key(&f), key);
    }
}
