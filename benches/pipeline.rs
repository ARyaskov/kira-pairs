//! End-to-end library benchmark on a small in-memory dataset: generate ->
//! sort -> dedup -> stats. Large-file numbers come from
//! `scripts/bench_end_to_end.sh`, not from here.

use std::hint::black_box;
use std::sync::Arc;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use kira_pairs::chroms::ChromDict;
use kira_pairs::dedup::{DedupConfig, Deduper};
use kira_pairs::generate::{GenerateConfig, generate};
use kira_pairs::io::buffered::LineBlock;
use kira_pairs::memory::MemoryBudget;
use kira_pairs::pairs::columns::ColumnMap;
use kira_pairs::sort::external::{ExternalSorter, SortConfig};
use kira_pairs::stats::StatsAccumulator;

fn dataset(n: u64) -> Vec<u8> {
    let cfg = GenerateConfig {
        records: n,
        ..Default::default()
    };
    let mut v = Vec::new();
    generate(&cfg, &mut v).unwrap();
    // strip header
    let body_start = v
        .iter()
        .enumerate()
        .find(|(_, _)| false)
        .map(|(i, _)| i)
        .unwrap_or_else(|| {
            let mut pos = 0;
            while pos < v.len() {
                let end = v[pos..]
                    .iter()
                    .position(|b| *b == b'\n')
                    .map(|p| pos + p + 1)
                    .unwrap_or(v.len());
                if v[pos] != b'#' {
                    return pos;
                }
                pos = end;
            }
            v.len()
        });
    v[body_start..].to_vec()
}

fn bench_pipeline(c: &mut Criterion) {
    let n = 200_000u64;
    let body = dataset(n);
    let mut g = c.benchmark_group("pipeline_in_memory");
    g.sample_size(10);
    g.throughput(Throughput::Elements(n));
    g.bench_function("sort_dedup_stats_200k", |b| {
        b.iter(|| {
            let dict = Arc::new(ChromDict::with_names(["!"]));
            let budget = MemoryBudget::new(1 << 30, 4, 2).unwrap();
            let cfg = SortConfig::new(4, budget);
            let mut sorter =
                ExternalSorter::new(cfg, ColumnMap::standard(), Arc::clone(&dict)).unwrap();
            let n_lines = body.iter().filter(|b| **b == b'\n').count() as u64;
            sorter
                .push_block(LineBlock {
                    data: body.clone(),
                    first_line: 1,
                    n_lines,
                    index: 0,
                })
                .unwrap();
            let (mut stream, _) = sorter.finish().unwrap();
            let mut dedup = Deduper::new(DedupConfig::default(), &dict);
            let mut stats = StatsAccumulator::with_dict(&dict);
            let mut sink = |e: kira_pairs::dedup::Emitted<'_>| {
                stats.observe(
                    e.key,
                    None,
                    matches!(e.outcome, kira_pairs::dedup::Outcome::Duplicate { .. }),
                );
                Ok(())
            };
            while let Some((k, l)) = stream.next_record().unwrap() {
                dedup.push(k, l, k.seq, &mut sink).unwrap();
            }
            dedup.finish(&mut sink).unwrap();
            black_box(stats.total)
        })
    });
    g.finish();
}

criterion_group!(benches, bench_pipeline);
criterion_main!(benches);
