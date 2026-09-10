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

(filled in below by the release run)

## Results

(filled in below by the release run)
