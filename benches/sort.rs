//! Microbenchmarks: sort key comparison and in-memory sorting.

use std::hint::black_box;
use std::sync::Arc;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use kira_pairs::chroms::ChromDict;
use kira_pairs::pairs::record::PairKey;
use kira_pairs::sort::key::SortKeyContext;

fn keys(n: usize) -> Vec<PairKey> {
    let mut x = 0x9e37_79b9u64;
    (0..n)
        .map(|i| {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            PairKey {
                seq: i as u64,
                pos1: x % 250_000_000,
                pos2: (x >> 20) % 250_000_000,
                chrom1: (x % 24) as u32,
                chrom2: ((x >> 8) % 24) as u32,
                pair_type: *b"UU\0\0\0\0\0\0",
                strand1: b'+',
                strand2: b'-',
                flags: 0,
                pair_type_len: 2,
            }
        })
        .collect()
}

fn bench_compare(c: &mut Criterion) {
    let dict = ChromDict::with_names((1..=24).map(|i| format!("chr{i}")));
    let ctx = SortKeyContext::new(Arc::new(dict.ranks()), Arc::new(Vec::new()));
    let ks = keys(10_000);
    let mut g = c.benchmark_group("sort_key");
    g.throughput(Throughput::Elements(ks.len() as u64));
    g.bench_function("cmp_inline", |b| {
        b.iter(|| {
            let mut acc = 0usize;
            for w in ks.windows(2) {
                if ctx
                    .cmp_inline(black_box(&w[0]), black_box(&w[1]))
                    .is_some_and(|o| o.is_lt())
                {
                    acc += 1;
                }
            }
            black_box(acc)
        })
    });
    g.finish();
    let mut g = c.benchmark_group("in_memory_sort");
    for &n in &[100_000usize, 1_000_000] {
        let ks = keys(n);
        g.throughput(Throughput::Elements(n as u64));
        g.bench_function(format!("sort_unstable_{n}"), |b| {
            b.iter_batched(
                || ks.clone(),
                |mut v| {
                    v.sort_unstable_by(|a, b| {
                        ctx.cmp_inline(a, b).unwrap_or(std::cmp::Ordering::Equal)
                    });
                    black_box(v.len())
                },
                criterion::BatchSize::LargeInput,
            )
        });
    }
    g.finish();
}

criterion_group!(benches, bench_compare);
criterion_main!(benches);
