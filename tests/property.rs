//! Property-based tests for the core invariants.

use std::sync::Arc;

use kira_pairs::chroms::{ChromDict, ChromOrder, ChromSizes};
use kira_pairs::dedup::{Clustering, DedupConfig, Deduper, Emitted, Method, Outcome};
use kira_pairs::flip::Flipper;
use kira_pairs::pairs::columns::ColumnMap;
use kira_pairs::pairs::record::{PairRecordRef, parse_line, split_fields};
use kira_pairs::sort::key::SortKeyContext;
use kira_pairs::stats::{DistBins, StatsAccumulator};
use proptest::prelude::*;

fn chrom_name() -> impl Strategy<Value = String> {
    prop_oneof![
        Just("chr1".to_string()),
        Just("chr2".to_string()),
        Just("chr10".to_string()),
        Just("chrX".to_string()),
        Just("scaffold_7".to_string()),
    ]
}

fn strand() -> impl Strategy<Value = char> {
    prop_oneof![Just('+'), Just('-')]
}

fn pair_line() -> impl Strategy<Value = String> {
    (
        "[A-Z0-9]{4,12}",
        chrom_name(),
        1u64..1_000_000,
        chrom_name(),
        1u64..1_000_000,
        strand(),
        strand(),
        prop_oneof![Just("UU"), Just("UR"), Just("RU"), Just("DD")],
        0u32..61,
        0u32..61,
    )
        .prop_map(|(id, c1, p1, c2, p2, s1, s2, pt, m1, m2)| {
            format!("{id}\t{c1}\t{p1}\t{c2}\t{p2}\t{s1}\t{s2}\t{pt}\t{m1}\t{m2}")
        })
}

fn columns() -> ColumnMap {
    ColumnMap::from_names(
        [
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
        ]
        .iter()
        .map(|s| s.to_string())
        .collect(),
    )
    .unwrap()
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn flip_is_an_involution(line in pair_line()) {
        let cs = ChromSizes::parse("chr1\t10\nchr2\t10\nchr10\t10\nchrX\t10\n", "t").unwrap();
        let cols = columns();
        let dict = Arc::new(ChromDict::new());
        let flipper = Flipper::new(ChromOrder::from_chromsizes(&cs), &cols, Arc::clone(&dict));
        let mut ends = Vec::new();
        split_fields(line.as_bytes(), &mut ends);
        let mut once = Vec::new();
        flipper.flip_line(&PairRecordRef::new(line.as_bytes(), &ends), &mut once);
        let mut ends2 = Vec::new();
        split_fields(&once, &mut ends2);
        let mut twice = Vec::new();
        flipper.flip_line(&PairRecordRef::new(&once, &ends2), &mut twice);
        prop_assert_eq!(twice, line.as_bytes().to_vec());
        // Keys agree with lines.
        let k = parse_line(line.as_bytes(), &cols, &dict, 0, &mut ends, Default::default).unwrap();
        let fk = Flipper::flip_key(&k);
        let k_once = parse_line(&once, &cols, &dict, 0, &mut ends, Default::default).unwrap();
        prop_assert_eq!(fk, k_once);
    }

    #[test]
    fn a_flipped_record_is_in_order(line in pair_line()) {
        let cs = ChromSizes::parse("chr1\t10\nchr2\t10\nchr10\t10\nchrX\t10\n", "t").unwrap();
        let cols = columns();
        let dict = Arc::new(ChromDict::new());
        let mut flipper = Flipper::new(ChromOrder::from_chromsizes(&cs), &cols, Arc::clone(&dict));
        let mut ends = Vec::new();
        let k = parse_line(line.as_bytes(), &cols, &dict, 0, &mut ends, Default::default).unwrap();
        if flipper.needs_flip(&k) {
            let fk = Flipper::flip_key(&k);
            // After flipping, annotated pairs never need flipping again.
            let both = flipper.enum_of(fk.chrom1).is_some() && flipper.enum_of(fk.chrom2).is_some();
            if both {
                prop_assert!(!flipper.needs_flip(&fk));
            }
        }
    }

    #[test]
    fn sort_order_is_total_and_lexicographic(lines in prop::collection::vec(pair_line(), 1..200)) {
        let cols = columns();
        let dict = ChromDict::new();
        let mut ends = Vec::new();
        let mut keyed: Vec<_> = lines
            .iter()
            .enumerate()
            .map(|(i, l)| (parse_line(l.as_bytes(), &cols, &dict, i as u64, &mut ends, Default::default).unwrap(), l.clone()))
            .collect();
        let ctx = SortKeyContext::new(Arc::new(dict.ranks()), Arc::new(Vec::new()));
        keyed.sort_by(|a, b| ctx.cmp_full(&a.0, a.1.as_bytes(), &b.0, b.1.as_bytes()));
        for w in keyed.windows(2) {
            let fa: Vec<&str> = w[0].1.split('\t').collect();
            let fb: Vec<&str> = w[1].1.split('\t').collect();
            let ka = (fa[1].as_bytes(), fa[3].as_bytes(), fa[2].parse::<u64>().unwrap(), fa[4].parse::<u64>().unwrap(), fa[7].as_bytes());
            let kb = (fb[1].as_bytes(), fb[3].as_bytes(), fb[2].parse::<u64>().unwrap(), fb[4].parse::<u64>().unwrap(), fb[7].as_bytes());
            prop_assert!(ka <= kb, "{:?} > {:?}", w[0].1, w[1].1);
            if ka == kb {
                prop_assert!(w[0].0.seq < w[1].0.seq);
            }
        }
    }

    #[test]
    fn dedup_output_has_no_close_unique_pairs(
        pos in prop::collection::vec((0u64..200, 0u64..200, 0u8..4), 1..300),
        max_mismatch in 0u64..6,
        sum in any::<bool>(),
    ) {
        // Build a sorted block on chr1/chr1.
        let mut recs: Vec<(u64, u64, u8)> = pos;
        recs.sort();
        let lines: Vec<String> = recs
            .iter()
            .enumerate()
            .map(|(i, (p1, p2, s))| format!("r{i}\tchr1\t{p1}\tchr1\t{p2}\t{}\t{}\tUU", if s & 1 == 0 { '+' } else { '-' }, if s & 2 == 0 { '+' } else { '-' }))
            .collect();
        let dict = ChromDict::with_names(["!"]);
        let cols = ColumnMap::standard();
        let cfg = DedupConfig {
            max_mismatch,
            method: if sum { Method::Sum } else { Method::Max },
            clustering: Clustering::Transitive,
            ..Default::default()
        };
        let mut d = Deduper::new(cfg, &dict);
        let mut kept: Vec<(u64, u64, u8, u8)> = Vec::new();
        let mut n_out = 0usize;
        let mut sink = |e: Emitted<'_>| {
            n_out += 1;
            if e.outcome == Outcome::Unique {
                kept.push((e.key.pos1, e.key.pos2, e.key.strand1, e.key.strand2));
            }
            Ok(())
        };
        let mut ends = Vec::new();
        for (i, l) in lines.iter().enumerate() {
            let k = parse_line(l.as_bytes(), &cols, &dict, i as u64, &mut ends, Default::default).unwrap();
            d.push(&k, l.as_bytes(), i as u64, &mut sink).unwrap();
        }
        d.finish(&mut sink).unwrap();
        prop_assert_eq!(n_out, lines.len());
        // No two kept records with equal strands are within reach of each other.
        for i in 0..kept.len() {
            for j in i + 1..kept.len() {
                let (a, b) = (kept[i], kept[j]);
                if a.2 != b.2 || a.3 != b.3 {
                    continue;
                }
                let d1 = a.0.abs_diff(b.0);
                let d2 = a.1.abs_diff(b.1);
                let close = if sum { d1 + d2 <= max_mismatch } else { d1.max(d2) <= max_mismatch };
                prop_assert!(!close, "kept records {:?} and {:?} are within {} ({})", a, b, max_mismatch, if sum { "sum" } else { "max" });
            }
        }
    }

    #[test]
    fn stats_counts_are_conserved(lines in prop::collection::vec(pair_line(), 0..100)) {
        let cols = columns();
        let dict = ChromDict::with_names(["!"]);
        let mut acc = StatsAccumulator::new(DistBins::default(), dict.get(b"!"));
        let mut ends = Vec::new();
        for (i, l) in lines.iter().enumerate() {
            let k = parse_line(l.as_bytes(), &cols, &dict, i as u64, &mut ends, Default::default).unwrap();
            acc.observe_plain(&k);
        }
        prop_assert_eq!(acc.total, lines.len() as u64);
        prop_assert_eq!(acc.total_mapped + acc.total_unmapped + acc.total_single_sided_mapped, acc.total);
        prop_assert_eq!(acc.total_dups + acc.total_nodups, acc.total_mapped);
        prop_assert_eq!(acc.cis + acc.trans, acc.total_nodups);
        let snap = acc.snapshot(&dict);
        let dist_total: u64 = snap.dist_freq.iter().map(|v| v.iter().sum::<u64>()).sum();
        prop_assert_eq!(dist_total, acc.cis);
        let pt_total: u64 = snap.pair_types.iter().map(|(_, n)| n).sum();
        prop_assert_eq!(pt_total, acc.total);
    }
}
