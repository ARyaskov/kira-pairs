# Benchmarks

Results are filled in by `scripts/bench_end_to_end.sh`; the table below is
the record of the runs made for the 0.1.0 release. Do not extrapolate from
the small fixtures in `tests/`: the numbers here use datasets larger than
the CPU caches and, where noted, larger than the memory budget.

## Methodology

* Synthetic data from `kira-pairs generate --seed 42 --duplicate-rate 0.15`
  (24 chromosomes, 100 Mb each, 75% cis with a heavy-tailed distance
  distribution, 2% unmapped, 24-character read IDs), written plain and as
  BGZF.
* Both tools run with the same `--threads`/`--nproc`, `--memory` and
  `--tmpdir`; page cache dropped between runs when permitted (otherwise
  warm-cache numbers are reported and labelled).
* `/usr/bin/time -f "%e %U %S %M"` for wall seconds, user/system CPU seconds
  and peak RSS (KiB). Temporary disk usage from `--metrics`
  (`temporary_bytes_written`).
* Outputs of the two tools are compared for identity after each pair of runs.
* Microbenchmarks: `cargo bench` (Criterion) in `benches/`.

## Environment

| Item | Value |
| --- | --- |
| CPU | Intel Core i7-12700 (8 P-cores + 4 E-cores, 20 hardware threads) |
| RAM | 30 GiB |
| Storage | NVMe SSD (ext4 on LVM); temporary runs written to the same volume |
| OS | Linux 7.0.0-31-generic (Ubuntu), x86_64 |
| Rust | rustc 1.98.0 (edition 2024), default portable build (`cargo build --release`, no `target-cpu=native`) |
| pairtools | 1.1.3 (Python 3.14, pandas 2.3, pysam 0.24; `lz4c` absent so pairtools sort compresses runs with gzip) |
| Cache state | warm (page cache could not be dropped without root); every input was read once by the generator before timing |

## Microbenchmarks (Criterion, single thread)

`cargo bench`, mean of 100 samples; throughput in millions of elements or MB of input text per second.

| Benchmark | Mean | Throughput |
| --- | ---: | ---: |
| `row_parsing/split_fields` (memchr tab scan, 10k rows) | 314 µs | 1817 MB/s |
| `row_parsing/parse_line_to_key` (full key parse, 10k rows) | 880 µs | 649 MB/s (11.4 M rows/s) |
| `row_parsing/naive_split_collect` (`split().collect::<Vec>()`, for comparison) | 999 µs | 572 MB/s |
| `int_parsing/parse_u64` (10k numbers) | 47 µs | 213 M/s |
| `int_parsing/std_str_parse` (for comparison) | 91 µs | 110 M/s |
| `chrom_dict/intern_hit` (10k lookups, RwLock + HashMap) | 191 µs | 52 M/s |
| `sort_key/cmp_inline` (10k comparisons) | 13.5 µs | 743 M/s |
| `in_memory_sort/sort_unstable_100000` | 8.6 ms | 11.6 M keys/s |
| `in_memory_sort/sort_unstable_1000000` | 95 ms | 10.5 M keys/s |
| `dedup/sweep_max3` (200k sorted records, 15% near-duplicates) | 11.6 ms | 17.3 M records/s |
| `bin_key/global_bin_ids` (100k pairs) | 248 µs | 403 M/s |
| `pipeline_in_memory/sort_dedup_stats_200k` (4 threads, incl. thread start-up) | 67 ms | 3.0 M records/s |

Per record, the full key parse costs roughly 90 ns, a sort comparison 1.3 ns,
and a dedup step 58 ns; the parser therefore sustains well over 10 M
records/s per core, which is why decompression and compression, not parsing,
dominate the end-to-end numbers below.

## End-to-end results

(filled in below by the release run)
