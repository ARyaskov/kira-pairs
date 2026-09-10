//! Microbenchmarks: TSV row parsing, integer parsing, dictionary lookup.

use std::hint::black_box;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use kira_pairs::chroms::ChromDict;
use kira_pairs::pairs::columns::ColumnMap;
use kira_pairs::pairs::record::{parse_line, split_fields};
use kira_pairs::util::int::parse_u64;

fn sample_lines(n: usize) -> Vec<Vec<u8>> {
    (0..n)
        .map(|i| {
            format!(
                "READ{i:012}\tchr{}\t{}\tchr{}\t{}\t+\t-\tUU\t60\t60",
                i % 23 + 1,
                1_000_000 + i * 37,
                i % 19 + 1,
                2_000_000 + i * 11
            )
            .into_bytes()
        })
        .collect()
}

fn bench_row_parsing(c: &mut Criterion) {
    let lines = sample_lines(10_000);
    let bytes: usize = lines.iter().map(|l| l.len() + 1).sum();
    let dict = ChromDict::new();
    let cols = ColumnMap::standard();
    let mut ends = Vec::with_capacity(16);
    let mut g = c.benchmark_group("row_parsing");
    g.throughput(Throughput::Bytes(bytes as u64));
    g.bench_function("split_fields", |b| {
        b.iter(|| {
            for l in &lines {
                split_fields(black_box(l), &mut ends);
                black_box(ends.len());
            }
        })
    });
    g.bench_function("parse_line_to_key", |b| {
        b.iter(|| {
            for (i, l) in lines.iter().enumerate() {
                let k = parse_line(
                    black_box(l),
                    &cols,
                    &dict,
                    i as u64,
                    &mut ends,
                    Default::default,
                )
                .unwrap();
                black_box(k);
            }
        })
    });
    g.bench_function("naive_split_collect", |b| {
        b.iter(|| {
            for l in &lines {
                let v: Vec<&[u8]> = l.split(|c| *c == b'\t').collect();
                black_box(v.len());
            }
        })
    });
    g.finish();
}

fn bench_int_parsing(c: &mut Criterion) {
    let nums: Vec<Vec<u8>> = (0..10_000u64)
        .map(|i| (i * 7_919 + 123_456_789).to_string().into_bytes())
        .collect();
    let mut g = c.benchmark_group("int_parsing");
    g.throughput(Throughput::Elements(nums.len() as u64));
    g.bench_function("parse_u64", |b| {
        b.iter(|| {
            for n in &nums {
                black_box(parse_u64(black_box(n)));
            }
        })
    });
    g.bench_function("std_str_parse", |b| {
        b.iter(|| {
            for n in &nums {
                let s = std::str::from_utf8(n).unwrap();
                black_box(s.parse::<u64>().ok());
            }
        })
    });
    g.finish();
}

fn bench_dict(c: &mut Criterion) {
    let dict = ChromDict::with_names((1..=24).map(|i| format!("chr{i}")));
    let names: Vec<Vec<u8>> = (0..10_000)
        .map(|i| format!("chr{}", i % 24 + 1).into_bytes())
        .collect();
    let mut g = c.benchmark_group("chrom_dict");
    g.throughput(Throughput::Elements(names.len() as u64));
    g.bench_function("intern_hit", |b| {
        b.iter(|| {
            for n in &names {
                black_box(dict.intern(black_box(n)));
            }
        })
    });
    g.finish();
}

criterion_group!(benches, bench_row_parsing, bench_int_parsing, bench_dict);
criterion_main!(benches);
