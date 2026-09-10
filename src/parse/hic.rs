//! Pair classification: walks, rescue of single ligations, flipping.

use crate::parse::alignment::Alignment;

/// Walks policy (pairtools `--walks-policy`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WalksPolicy {
    /// Mask unrescuable walks.
    Mask,
    /// 5'-most alignment on each side.
    FiveAny,
    /// 5'-most unique alignment on each side.
    FiveUnique,
    /// 3'-most alignment on each side.
    ThreeAny,
    /// 3'-most unique alignment on each side.
    ThreeUnique,
    /// Report all alignments (not implemented).
    All,
}

/// pairtools `pair_index`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairIndex {
    /// `(1, "R1-2")`
    R12,
    /// `(1, "R1")`
    R1,
    /// `(1, "R2")`
    R2,
}

impl PairIndex {
    /// The `walk_pair_type` label.
    pub fn label(self) -> &'static str {
        match self {
            PairIndex::R12 => "R1-2",
            PairIndex::R1 => "R1",
            PairIndex::R2 => "R2",
        }
    }
}

/// pairtools `_convert_gaps_into_alignments`, including its iteration
/// quirk (the loop bound is fixed before insertions, so alignments shifted
/// past the original length are not revisited).
fn convert_gaps(algns: &mut Vec<Alignment>, max_gap: u64) {
    if algns.len() == 1 && !algns[0].is_mapped {
        return;
    }
    let mut last_5_pos = 0u64;
    let n = algns.len();
    for i in 0..n {
        let d5 = algns[i].dist_to_5;
        let span = algns[i].algn_read_span;
        let read_len = algns[i].read_len;
        if d5 > last_5_pos && d5 - last_5_pos > max_gap {
            let mut gap = Alignment::empty();
            gap.dist_to_5 = last_5_pos;
            gap.algn_read_span = d5 - last_5_pos;
            gap.read_len = read_len;
            gap.dist_to_3 = read_len.saturating_sub(d5);
            last_5_pos = d5 + span;
            algns.insert(i, gap);
        } else {
            last_5_pos = last_5_pos.max(d5 + span);
        }
    }
}

/// pairtools `normalize_alignment_list`.
pub fn normalize(algns: &mut Vec<Alignment>, side: u64, max_gap: Option<u64>) {
    if algns.is_empty() {
        algns.push(Alignment::empty());
    }
    algns.sort_by_key(|a| a.dist_to_5);
    if let Some(g) = max_gap {
        convert_gaps(algns, g);
    }
    let n = algns.len() as u64;
    for (i, a) in algns.iter_mut().enumerate() {
        a.read_side = Some(side);
        a.algn_idx = Some(i as u64);
        a.same_side_count = Some(n);
    }
}

/// pairtools `rescue_walk`. Returns the linear side (1 or 2) on success.
pub fn rescue_walk(
    algns1: &mut [Alignment],
    algns2: &mut [Alignment],
    max_molecule_size: u64,
) -> Option<u8> {
    let (n1, n2) = (algns1.len(), algns2.len());
    if n1 <= 1 && n2 <= 1 {
        return None;
    }
    if !((n1 == 1 && n2 == 2) || (n1 == 2 && n2 == 1)) {
        return None;
    }
    let first_is_chimeric = n1 > 1;
    let (chim5, chim3, linear) = if first_is_chimeric {
        (&algns1[0], &algns1[1], &algns2[0])
    } else {
        (&algns2[0], &algns2[1], &algns1[0])
    };
    if !(linear.is_mapped && linear.is_unique) {
        return None;
    }
    let mut can_rescue = true;
    if chim3.is_mapped && chim5.is_unique {
        can_rescue &= chim3.chrom == linear.chrom;
        can_rescue &= chim3.strand != linear.strand;
        if linear.strand == b'+' {
            can_rescue &= linear.pos5 < chim3.pos5;
        } else {
            can_rescue &= linear.pos5 > chim3.pos5;
        }
        let molecule_size: i128 = if linear.strand == b'+' {
            i128::from(chim3.pos5) - i128::from(linear.pos5)
                + i128::from(chim3.dist_to_5)
                + i128::from(linear.dist_to_5)
        } else {
            i128::from(linear.pos5) - i128::from(chim3.pos5)
                + i128::from(chim3.dist_to_5)
                + i128::from(linear.dist_to_5)
        };
        can_rescue &= molecule_size <= i128::from(max_molecule_size);
    }
    if can_rescue {
        if first_is_chimeric {
            algns1[1].kind = b'X';
            algns2[0].kind = b'R';
            Some(1)
        } else {
            algns1[0].kind = b'R';
            algns2[1].kind = b'X';
            Some(2)
        }
    } else {
        None
    }
}

/// pairtools `parse_read` for the standard policies. Returns the two
/// reported alignments and the pair index; `algns1`/`algns2` are
/// normalised in place.
pub fn parse_read(
    algns1: &mut Vec<Alignment>,
    algns2: &mut Vec<Alignment>,
    max_molecule_size: u64,
    max_inter_align_gap: Option<u64>,
    policy: WalksPolicy,
) -> (Alignment, Alignment, PairIndex) {
    normalize(algns1, 1, max_inter_align_gap);
    normalize(algns2, 2, max_inter_align_gap);
    let mut pair_index = PairIndex::R12;
    let is_chimeric_1 = algns1.len() > 1;
    let is_chimeric_2 = algns2.len() > 1;
    if !(is_chimeric_1 || is_chimeric_2) {
        return (algns1[0].clone(), algns2[0].clone(), pair_index);
    }
    let rescued = rescue_walk(algns1, algns2, max_molecule_size);
    if let Some(side) = rescued {
        pair_index = if side == 1 {
            PairIndex::R1
        } else {
            PairIndex::R2
        };
        return (algns1[0].clone(), algns2[0].clone(), pair_index);
    }
    let pick = |algns: &[Alignment], policy: WalksPolicy| -> Alignment {
        match policy {
            WalksPolicy::FiveAny => algns[0].clone(),
            WalksPolicy::ThreeAny => algns[algns.len() - 1].clone(),
            WalksPolicy::FiveUnique => algns
                .iter()
                .find(|a| a.is_mapped && a.is_unique)
                .cloned()
                .unwrap_or_else(|| algns[0].clone()),
            WalksPolicy::ThreeUnique => algns
                .iter()
                .rev()
                .find(|a| a.is_mapped && a.is_unique)
                .cloned()
                .unwrap_or_else(|| algns[algns.len() - 1].clone()),
            WalksPolicy::Mask | WalksPolicy::All => {
                let mut a = algns[0].clone();
                a.mask();
                a.kind = b'W';
                a
            }
        }
    };
    (pick(algns1, policy), pick(algns2, policy), pair_index)
}

/// pairtools `check_pair_order`: true when the pair is in upper-triangular
/// order. `ref_enum[i]` is the enumeration value of reference `i`.
pub fn check_pair_order(a1: &Alignment, a2: &Alignment, ref_enum: &[u32]) -> bool {
    let mut correct = (a1.is_mapped, a1.is_unique) <= (a2.is_mapped, a2.is_unique);
    if let (Some(c1), Some(c2)) = (a1.chrom, a2.chrom) {
        let e1 = ref_enum.get(c1).copied().unwrap_or(u32::MAX);
        let e2 = ref_enum.get(c2).copied().unwrap_or(u32::MAX);
        correct = (e1, a1.pos) <= (e2, a2.pos);
    }
    correct
}

#[cfg(test)]
mod tests {
    use super::*;

    fn mapped(chrom: usize, pos5: u64, strand: u8, d5: u64, span: u64, unique: bool) -> Alignment {
        let mut a = Alignment::empty();
        a.chrom = Some(chrom);
        a.pos5 = pos5;
        a.pos3 = if strand == b'+' {
            pos5 + span - 1
        } else {
            pos5 + 1 - span
        };
        a.pos = pos5;
        a.strand = strand;
        a.is_mapped = true;
        a.is_unique = unique;
        a.mapq = if unique { 60 } else { 0 };
        a.dist_to_5 = d5;
        a.algn_read_span = span;
        a.algn_ref_span = span;
        a.read_len = 100;
        a.kind = if unique { b'U' } else { b'M' };
        a
    }

    #[test]
    fn simple_pair() {
        let mut a1 = vec![mapped(0, 100, b'+', 0, 100, true)];
        let mut a2 = vec![mapped(0, 500, b'-', 0, 100, true)];
        let (h1, h2, pi) = parse_read(&mut a1, &mut a2, 750, Some(20), WalksPolicy::FiveUnique);
        assert_eq!((h1.kind, h2.kind), (b'U', b'U'));
        assert_eq!(pi, PairIndex::R12);
        assert!(check_pair_order(&h1, &h2, &[1]));
        assert!(!check_pair_order(&h2, &h1, &[1]));
    }

    #[test]
    fn rescue_single_ligation() {
        // Read 1 chimeric: 5' part on chr0 at 1000 (+), 3' part at 1200 (-);
        // read 2 linear on chr0 at 1150 (-) pointing towards the 3' part... pairtools:
        // linear must be on the same chrom as chim3, opposite strand, and ordered.
        let mut a1 = vec![
            mapped(0, 1000, b'+', 0, 60, true),
            mapped(0, 1300, b'-', 60, 40, true),
        ];
        let mut a2 = vec![mapped(0, 1250, b'+', 0, 100, true)];
        let (h1, h2, pi) = parse_read(&mut a1, &mut a2, 750, Some(20), WalksPolicy::FiveUnique);
        assert_eq!((h1.kind, h2.kind), (b'U', b'R'));
        assert_eq!(pi, PairIndex::R1);
        assert_eq!(h1.pos5, 1000);
    }

    #[test]
    fn non_unique_5prime_part_auto_rescues() {
        // pairtools checks the 5' part's uniqueness: a multi-mapped 5' part
        // with a unique 3' part is rescued without geometry checks.
        let mut a1 = vec![
            mapped(0, 1000, b'+', 0, 50, false),
            mapped(1, 5000, b'+', 50, 50, true),
        ];
        let mut a2 = vec![mapped(2, 100, b'-', 0, 100, true)];
        let (h1, h2, pi) = parse_read(&mut a1, &mut a2, 750, Some(20), WalksPolicy::FiveUnique);
        assert_eq!((h1.kind, h2.kind), (b'M', b'R'));
        assert_eq!(pi, PairIndex::R1);
    }

    #[test]
    fn unrescuable_walk_policies() {
        // Three alignments on side 1: never rescuable.
        let mut a1 = vec![
            mapped(0, 1000, b'+', 0, 30, false),
            mapped(1, 5000, b'+', 30, 30, true),
            mapped(1, 6000, b'+', 60, 40, true),
        ];
        let mut a2 = vec![mapped(2, 100, b'-', 0, 100, true)];
        let (h1, h2, _) = parse_read(
            &mut a1.clone(),
            &mut a2.clone(),
            750,
            Some(20),
            WalksPolicy::FiveUnique,
        );
        assert_eq!((h1.kind, h2.kind), (b'U', b'U'));
        assert_eq!((h1.chrom, h1.pos5), (Some(1), 5000));
        let (h1, _, _) = parse_read(
            &mut a1.clone(),
            &mut a2.clone(),
            750,
            Some(20),
            WalksPolicy::FiveAny,
        );
        assert_eq!(h1.kind, b'M');
        let (h1, _, _) = parse_read(
            &mut a1.clone(),
            &mut a2.clone(),
            750,
            Some(20),
            WalksPolicy::ThreeAny,
        );
        assert_eq!((h1.chrom, h1.pos5), (Some(1), 6000));
        let (h1, h2, _) = parse_read(&mut a1, &mut a2, 750, Some(20), WalksPolicy::Mask);
        assert_eq!((h1.kind, h2.kind), (b'W', b'W'));
        assert_eq!(h1.chrom, None);
    }

    #[test]
    fn gaps_become_null_alignments() {
        let mut a1 = vec![mapped(0, 1000, b'+', 30, 70, true)];
        let mut a2 = vec![mapped(0, 5000, b'-', 0, 100, true)];
        let (h1, h2, _) = parse_read(&mut a1, &mut a2, 750, Some(20), WalksPolicy::FiveUnique);
        // The leading 30 bp gap becomes a null alignment; the null 5' part
        // auto-rescues the walk, so pairtools reports the null side (NR).
        assert_eq!(a1.len(), 2);
        assert_eq!(a1[0].kind, b'N');
        assert_eq!((h1.kind, h2.kind), (b'N', b'R'));
        // A small gap is ignored.
        let mut a1 = vec![mapped(0, 1000, b'+', 10, 90, true)];
        let mut a2 = vec![mapped(0, 5000, b'-', 0, 100, true)];
        let (h1, h2, _) = parse_read(&mut a1, &mut a2, 750, Some(20), WalksPolicy::FiveUnique);
        assert_eq!(a1.len(), 1);
        assert_eq!((h1.kind, h2.kind), (b'U', b'U'));
    }
}
