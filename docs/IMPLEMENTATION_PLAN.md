# kira-pairs v0.1 implementation plan

Reference behaviour: pairtools 1.1.3 (source inspected, installed as an oracle
for differential tests). Target toolchain: Rust 1.98.0, edition 2024, stable only.

## Guiding decisions

* **One internal record shape for every stage.** A `PairKey` (chrom ids,
  positions, strands, pair type prefix, monotonic sequence number) plus the
  original TAB-separated line bytes. Every hot stage (sort, dedup, stats, bin)
  works from the key; the line is carried untouched so unknown extra columns
  round-trip. Text is re-parsed only for columns a stage explicitly asks for
  (e.g. `--extra-col-pair`, `mapq1`), and only for that column.
* **Chromosome names are dictionary-encoded** once per process. Sort order is
  pairtools' lexicographic byte order, computed from a rank table over the
  dictionary (never natural/"chr2 < chr10" order).
* **Bounded memory everywhere.** `--memory` is a real budget: the external
  sorter sizes runs from it, the binner spills sorted runs from it, channels
  are bounded and small, and dedup keeps only the open positional window.
* **Determinism.** Sequence numbers are assigned in input order regardless of
  parser thread count; sort ties are broken by sequence number (pairtools uses
  `sort --stable`), merges are exact, and per-thread stats accumulators are
  merged in a fixed order.

## Milestones

1. Foundation: CLI, errors, chromsizes, header, reader/writer, compression
   (gzip/BGZF/lz4 auto-detect, parallel BGZF writer), temp dirs, logging.
2. Sort: run generation (parallel parse -> bounded chunk -> parallel sort ->
   lz4 block-compressed private run format) and loser-tree k-way merge with
   multi-pass parallel merging when the fan-in is large.
3. Dedup: positional sweep with pos2 bucket index, union-find transitive
   clustering (pairtools default `scipy` backend), greedy mode (`cython`
   backend), max/sum metrics, extra-col pairs, parent ids, output routing,
   embedded stats.
4. Stats: incremental `StatsAccumulator`, pairtools TSV/YAML output, JSON
   extension, merging, log-spaced distance bins.
5. Flip and select (expression compiler + evaluator).
6. Bin: streaming aggregation with spill-to-disk merge-reduce, COO and BG2
   outputs, multi-resolution in one pass.
7. Parse: BAM/SAM (noodles) -> pairs for standard paired-end Hi-C
   (walks policies, rescue of simple walks, extra columns, flipping).
8. Process: parse -> flip -> sort -> dedup -> stats/bin without text
   intermediates on disk.
9. Differential tests vs pairtools, benchmarks, docs.

## Private run format (unstable)

`KPRUN` magic, version, flags; then lz4 block-compressed blocks of
`[key][line_len][line]` records. Not a public format; may change any time.
