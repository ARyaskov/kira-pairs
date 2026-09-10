//! Microbenchmarks: dedup neighbourhood lookup on sorted synthetic data and
//! bin-key construction.

use std::hint::black_box;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use kira_pairs::binning::BinLayout;
use kira_pairs::chroms::{ChromDict, ChromSizes};
use kira_pairs::dedup::{DedupConfig, Deduper};
use kira_pairs::pairs::record::PairKey;

fn sorted_keys(n: usize, dup_rate: f64) -> Vec<(PairKey, Vec<u8>)> {
    let mut x = 0x2545_f491_4f6c_dd1du64;
    let mut out = Vec::with_capacity(n);
    let mut pos1 = 1u64;
    for i in 0..n {
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        let is_dup = !out.is_empty() && (x % 1000) as f64 / 1000.0 < dup_rate;
        // Keep pos1 non-decreasing (block-sorted input): duplicates reuse the
        // current pos1 with a small forward jitter, unique records advance.
        let (p1, p2) = if is_dup {
            let prev: &(PairKey, Vec<u8>) = &out[out.len() - 1];
            let p1 = prev.0.pos1 + (x % 3);
            pos1 = pos1.max(p1);
            (p1, prev.0.pos2 + ((x >> 3) % 3))
        } else {
            pos1 += 3 + (x >> 5) % 40;
            (pos1, 1_000_000 + (x >> 9) % 100_000_000)
        };
        let k = PairKey {
            seq: i as u64,
            pos1: p1,
            pos2: p2,
            chrom1: 1,
            chrom2: 1,
            pair_type: *b"UU\0\0\0\0\0\0",
            strand1: b'+',
            strand2: b'-',
            flags: 0,
            pair_type_len: 2,
        };
        let line = format!("r{i}\tchr1\t{p1}\tchr1\t{p2}\t+\t-\tUU").into_bytes();
        out.push((k, line));
    }
    out
}

fn bench_dedup(c: &mut Criterion) {
    let dict = ChromDict::with_names(["!", "chr1"]);
    let data = sorted_keys(200_000, 0.15);
    let mut g = c.benchmark_group("dedup");
    g.throughput(Throughput::Elements(data.len() as u64));
    g.bench_function("sweep_max3", |b| {
        b.iter(|| {
            let mut d = Deduper::new(DedupConfig::default(), &dict);
            let mut n = 0u64;
            let mut sink = |_e: kira_pairs::dedup::Emitted<'_>| {
                n += 1;
                Ok(())
            };
            for (i, (k, l)) in data.iter().enumerate() {
                d.push(k, l, i as u64, &mut sink).unwrap();
            }
            d.finish(&mut sink).unwrap();
            black_box(n)
        })
    });
    g.finish();
}

fn bench_bin_key(c: &mut Criterion) {
    let cs = ChromSizes::parse(
        &(1..=24)
            .map(|i| format!("chr{i}\t{}\n", 200_000_000 - i * 1000))
            .collect::<String>(),
        "t",
    )
    .unwrap();
    let layout = BinLayout::new(&cs, 10_000);
    let data = sorted_keys(100_000, 0.0);
    let mut g = c.benchmark_group("bin_key");
    g.throughput(Throughput::Elements(data.len() as u64));
    g.bench_function("global_bin_ids", |b| {
        b.iter(|| {
            let mut acc = 0u64;
            for (k, _) in &data {
                let c1 = (k.chrom1 as usize) % 24;
                let c2 = (k.chrom2 as usize + 3) % 24;
                let b1 =
                    layout.offsets[c1] + ((k.pos1 - 1) / 10_000).min(layout.chrom_bins(c1) - 1);
                let b2 =
                    layout.offsets[c2] + ((k.pos2 - 1) / 10_000).min(layout.chrom_bins(c2) - 1);
                acc = acc.wrapping_add(b1.min(b2) * 31 + b1.max(b2));
            }
            black_box(acc)
        })
    });
    g.finish();
}

criterion_group!(benches, bench_dedup, bench_bin_key);
criterion_main!(benches);
