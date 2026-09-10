# Changelog

## 0.1.0 (2026-09-10)

First release.

* `.pairs` reader/writer with header preservation, `#columns:`-driven column
  maps, plain/gzip/BGZF/LZ4 input by magic bytes, parallel BGZF output.
* `sort`: parallel external merge sort with pairtools block-sort semantics
  and stable tie-breaking; LZ4-compressed private run files; loser-tree
  k-way merge with parallel intermediate passes.
* `dedup`: pairtools-compatible duplicate detection (`max`/`sum` metrics,
  transitive or greedy clustering, `--extra-col-pair`, `--keep-parent-id`,
  output routing, embedded stats).
* `stats`: streaming pairtools-compatible statistics (TSV, YAML, JSON),
  merging of stats files.
* `flip`, `select` (safe expression engine), `bin` (COO/BG2, multi-resolution,
  bounded memory).
* `parse`: standard paired-end Hi-C SAM/BAM parsing (walk rescue, walks
  policies, extra columns, pairsam output).
* `process`: fused parse → flip → sort → dedup → stats/bin pipeline.
* `generate`: deterministic synthetic datasets for benchmarks.
