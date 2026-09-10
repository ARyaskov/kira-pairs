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

Dataset: `kira-pairs generate --records 10000000 --seed 42 --duplicate-rate 0.15`
(10 M records, 607 MB plain, 309 MB BGZF, 24 chromosomes), `--threads 16`
(`--nproc 16` for pairtools), `--memory 4G`, temporary files on the same NVMe
volume, warm page cache. `scripts/bench_end_to_end.sh --records 10000000
--threads 16 --memory 4G`. Sort outputs of both tools were verified identical
after the run; the dedup outputs differ by the chunk-boundary effect
described in [PAIRTOOLS_COMPATIBILITY.md](PAIRTOOLS_COMPATIBILITY.md) (see
"Dedup difference at scale" below).

| Command | Wall (s) | CPU (s) | Peak RSS (MB) | Records/s | Input MB/s |
| --- | ---: | ---: | ---: | ---: | ---: |
| kira sort (plain in, plain out) | 2.66 | 30.1 | 1811 | 3.76 M | 228 |
| kira sort (gz in, gz out) | 2.79 | 39.7 | 1822 | 3.58 M | 111 |
| pairtools sort (plain in, plain out) | 7.94 | 76.0 | 2113 | 1.26 M | 76 |
| pairtools sort (gz in, gz out) | 12.74 | 86.6 | 2113 | 0.78 M | 24 |
| kira dedup (plain) | 1.73 | 21.4 | 770 | 5.78 M | 351 |
| pairtools dedup (plain) | 46.61 | 47.9 | 176 | 0.21 M | 13 |
| kira stats | 1.55 | 22.5 | 629 | 6.45 M | 392 |
| pairtools stats | 10.06 | 14.8 | 235 | 0.99 M | 60 |
| kira sort \| dedup (gz in, gz out) | 3.56 | 46.6 | 1808 | 2.81 M | 87 |
| pairtools sort \| dedup (gz in, gz out) | 52.56 | 129.0 | 2113 | 0.19 M | 6 |

Peak RSS for kira-pairs is dominated by the sort run buffers, which are sized
from `--memory` (a 4G budget yields about 1.3 GB of run buffers, two runs
resident); the 10 M-record input fits in a single run, so the sorter never
touched disk here. pairtools sort delegates to GNU `sort -S 4G --parallel 16`
plus Python I/O and gzip-compressed chunks (`lz4c` was not installed).

| Workflow | pairtools 1.1.3 | kira-pairs 0.1.0 | Speed-up |
| --- | ---: | ---: | ---: |
| sort (plain in, plain out) | 7.9 s | 2.66 s | 3.0× |
| sort (gz in, gz out) | 12.7 s | 2.79 s | 4.6× |
| dedup (plain) | 46.6 s | 1.73 s | 26.9× |
| stats | 10.1 s | 1.55 s | 6.5× |
| sort \| dedup (gz in, gz out) | 52.6 s | 3.56 s | 14.8× |

Speed-ups are wall-clock ratios on this machine and dataset; pairtools
dedup and stats are single-threaded Python/pandas/scipy code, so the gap
there is mostly CPU work per record rather than parallelism. The fused
`process` pipeline was not benchmarked against a large BAM in this release
(no large aligner output was available offline); on the committed 400-read
fixture it produces identical output to `pairtools parse | sort | dedup`.

## Dedup difference at scale

On the 10 M-record dataset kira-pairs marks 1 470 982 duplicates and
pairtools 1 470 976 (0.0004% fewer). The 6 records found only by kira-pairs
were examined individually: each is within `--max-mismatch 3` of a record
that pairtools itself marks as a duplicate, 5 of the 6 are the first row of
a 10 000-row pairtools chunk (offsets 0-19 of the chunk), and their cluster
link is a duplicate in the previous chunk, which pairtools' 100-row
carry-over (non-duplicates only) does not retain. No record is marked by
pairtools and not by kira-pairs. This is pairtools' documented chunking
approximation, not a metric difference; on the 300 k-record differential
dataset and all fixtures the two tools agree exactly.

## Memory-bounded sort at scale

`scripts/verify_large_sort.sh --records 30000000 --memory 256M --threads 16
--max-fan-in 4`: 30 M records (1.83 GB of text) sorted with a 256 M budget and a
deliberately tiny fan-in to force intermediate merge passes.

(filled in below)
