# kira-pairs

**kira-pairs** is a fast, memory-bounded Rust implementation of the
performance-critical parts of the Hi-C `.pairs` workflow handled by
[pairtools](https://github.com/open2c/pairtools): sorting, duplicate removal,
statistics, flipping, filtering, contact binning, and standard paired-end
BAM parsing, plus a fused `process` pipeline that runs all of them without
text intermediates.

It reads and writes the 4DN/pairtools `.pairs` format (plain, gzip/BGZF, LZ4),
preserves headers and unknown extra columns, and reproduces pairtools 1.1.3
semantics wherever a command is implemented, so that
`pairtools` and `kira-pairs` commands can be mixed freely in a pipeline.

kira-pairs is an independent project. It is not an official Open2C or 4DN tool.

## Status

Version 0.1.0. See [docs/PAIRTOOLS_COMPATIBILITY.md](docs/PAIRTOOLS_COMPATIBILITY.md)
for the feature-by-feature compatibility matrix (every "exact" claim is backed
by an automated differential test against pairtools 1.1.3) and
[docs/BENCHMARKS.md](docs/BENCHMARKS.md) for measured performance.

Implemented in v0.1: `sort`, `flip`, `dedup`, `stats`, `select`, `bin`,
`parse` (standard paired-end Hi-C), `process` (fused pipeline), `generate`
(synthetic data). Not implemented: `parse2`, `phase`, `restrict`,
`filterbycov`, `scaling`, `sample`, `split`, `merge`, `header`, `markasdup`,
`--walks-policy all`, stats `--filter`, `.cool` output.

## Installation

Requires Rust 1.98.0 or newer (edition 2024, stable toolchain).

```bash
git clone https://github.com/riaskov/kira-pairs
cd kira-pairs
cargo build --release
./target/release/kira-pairs --help
```

The default build is portable. For a host-optimised binary:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

An optional `libdeflate` Cargo feature is reserved for a native BGZF codec
backend; the default build uses the pure-Rust `zlib-rs` backend and has no C
dependencies.

## Supported formats

| Input                    | Detection                | Output                      |
| ------------------------ | ------------------------ | --------------------------- |
| `.pairs` (plain text)    | magic bytes              | by extension                |
| `.pairs.gz` (BGZF/bgzip) | magic bytes, parallel    | `.gz`, `.bgz`, `.bgzf` → BGZF (parallel) |
| `.pairs.gz` (plain gzip) | magic bytes              | —                           |
| `.pairs.lz4` (LZ4 frame) | magic bytes              | `.lz4`                      |
| SAM text                 | magic bytes (`parse`)    | —                           |
| BAM                      | magic bytes (`parse`)    | —                           |

`-` (or an omitted path) means stdin/stdout; compression of stdin is detected
from its magic bytes. stdout is data only; diagnostics go to stderr.

## Quick start

```bash
# Sort, dedup and collect stats from an existing .pairs file
kira-pairs sort input.pairs.gz --threads 16 --memory 8G --tmpdir /nvme/tmp \
  | kira-pairs dedup --max-mismatch 3 --output-stats stats.txt -o nodups.pairs.gz

# Everything from a name-grouped BAM in one pass, no text intermediates
kira-pairs process input.bam --chroms-path hg38.chrom.sizes \
    --threads 24 --memory 12G --tmpdir /nvme/tmp \
    --drop-sam --max-mismatch 3 --resolution 10000 \
    --output-pairs contacts.nodups.pairs.gz \
    --output-stats contacts.stats \
    --output-bins contacts.10kb.tsv.gz

# Load the bins into cooler
cooler load -f coo hg38.chrom.sizes:10000 contacts.10kb.tsv.gz contacts.10kb.cool
```

## Commands

```
kira-pairs sort      [INPUT] -o OUT [--extra-col COL]... [--max-fan-in N]
kira-pairs flip      [INPUT] -o OUT -c chrom.sizes
kira-pairs dedup     [INPUT] -o OUT [--output-dups F] [--output-unmapped F] [--output-stats F]
                     [--max-mismatch N] [--method max|sum] [--clustering transitive|greedy]
                     [--mark-dups|--no-mark-dups] [--keep-parent-id] [--extra-col-pair C1 C2]...
kira-pairs stats     [INPUT] -o OUT [--yaml|--json] [--n-dist-bins-decade N]
kira-pairs stats     --merge STATS... -o OUT [--yaml]
kira-pairs select    CONDITION [INPUT] -o OUT [--output-rest F] [-t COL TYPE]... [--chrom-subset F]
kira-pairs bin       [INPUT] -o OUT -c chrom.sizes -r RES [-r RES]... [--format coo|bg2]
                     [--min-mapq N] [--pair-types LIST] [--bins-out F]
kira-pairs parse     [INPUT.bam|.sam] -o OUT -c chrom.sizes [--min-mapq N] [--walks-policy ...]
                     [--drop-sam] [--add-columns LIST] [--output-stats F] ...
kira-pairs process   INPUT.bam -c chrom.sizes [--output-pairs F] [--output-dups F]
                     [--output-stats F] [--output-bins F -r RES]... [dedup/parse options]
kira-pairs generate  --records N --seed S [--duplicate-rate X] [--sorted] -o OUT
```

Common options on every command:

```
--threads N               worker threads (default: all cores)
--memory SIZE             memory budget, e.g. 512M, 8G (default 2G)
--tmpdir PATH             scratch directory for temporary runs
--io-threads N            input decompression threads (default min(threads,4))
--compression-threads N   output compression threads (default threads)
--compression-level 1-9   gzip/BGZF level (default 6)
-v / -vv                  diagnostics; --quiet errors only
--metrics                 runtime counters on stderr at exit
--progress                periodic progress on stderr
```

### Migration from pairtools

| pairtools command        | kira-pairs equivalent                                   |
| ------------------------ | ------------------------------------------------------- |
| `pairtools sort`         | `kira-pairs sort`                                       |
| `pairtools flip`         | `kira-pairs flip`                                       |
| `pairtools dedup`        | `kira-pairs dedup`                                      |
| `pairtools stats`        | `kira-pairs stats`                                      |
| `pairtools select`       | `kira-pairs select`                                     |
| `pairtools parse`        | `kira-pairs parse`                                      |
| `pairtools parse \| sort \| dedup` | `kira-pairs process`                          |
| `cooler cload pairs`     | `kira-pairs bin` + `cooler load -f coo` (or `-f bg2`)   |

Option names follow pairtools where the option exists (`--max-mismatch`,
`--method`, `--keep-parent-id`, `--extra-col-pair`, `--walks-policy`,
`--add-columns`, `--drop-sam`, `-c/--chroms-path`, `--output-dups`, ...).
pairtools' `--backend scipy|sklearn|cython` is accepted as an alias for
`--clustering transitive|transitive|greedy`. pairtools' `--nproc` becomes
`--threads`.

Mixed pipelines work in both directions:

```bash
pairtools parse -c hg38.chrom.sizes input.bam | kira-pairs sort --memory 8G | kira-pairs dedup -o out.pairs.gz
kira-pairs parse -c hg38.chrom.sizes input.bam | kira-pairs sort | pairtools dedup -o out.pairs.gz
```

## Performance architecture

* **One decode, one representation.** Every stage works on a 48-byte
  `PairKey` (dictionary-encoded chromosome ids, `u64` positions, strands,
  inline pair type, input sequence number) plus the original line bytes.
  Extra columns are never parsed unless a stage asks for them, and the line is
  carried untouched so unknown columns round-trip.
* **Streaming, bounded pipelines.** Reader → parser workers → consumer stages
  are connected with bounded channels; the reader is throttled by
  back-pressure. Nothing grows with input size.
* **External merge sort.** Input is parsed in parallel into a bounded chunk,
  sorted with a parallel unstable sort on keys (ties broken by the sequence
  number, i.e. pairtools' `sort --stable`), and written as an LZ4
  block-compressed private run file with parallel compression. Runs are merged
  with a loser tree; if more runs exist than the fan-in allows, intermediate
  merge passes run in parallel. Files that fit in memory never touch disk.
* **Sweep-line dedup.** Sorted input is deduplicated with a positional window
  of `max_mismatch` on `pos1`, a `pos2` bucket index for dense loci, and a
  union-find for transitive clusters. Work per record is proportional to the
  local density, not to the block size.
* **Incremental stats and binning.** `StatsAccumulator::observe` is a
  lock-free counter update that runs inside dedup/process; the binner
  aggregates in a hash table that spills sorted runs when it exceeds its budget
  and merge-reduces at the end. Several resolutions are produced in one pass.
* **Parallel BGZF.** Output compression and BGZF input decompression use a
  worker pool; plain gzip input is decompressed on a dedicated thread that
  overlaps with parsing.

Internal run files (`KPRUN` format) are private, versioned and validated on
read; they are not a public storage format and may change between versions.

## Memory behaviour

`--memory` is a budget, not a hint. It is split into a record budget (sort
runs, bin tables), channel buffers, compression buffers and merge read-ahead;
the sorter sizes its runs so that at most two runs are resident (one being
filled, one being sorted/written). Peak RSS is predictable: a run of
`--memory 64M` over a 100 MB file stays far below 512 MB in the test suite,
and the number of runs and bytes written to `--tmpdir` are reported by
`--metrics`. The minimum accepted budget is 64M.

## Reproducibility

Output is deterministic: the same input and arguments give byte-identical
data regardless of `--threads`, chunk boundaries or scheduling. Sort ties keep
input order, merges are exact, and stats are accumulated in input order. BGZF
output carries no timestamps.

## Limitations (v0.1)

* `parse` covers the standard paired-end case only: `--walks-policy all`
  (complex walks), `parse2`, `--add-columns mismatches`,
  `--readid-transform` and `--output-parsed-alignments` are not implemented
  and are rejected with an error. CRAM input is not supported.
* `stats --filter`, by-tile duplicate stats and pairsam-specific SAM-flag
  marking of duplicates are not implemented.
* `sort --extra-col` compares extra columns lexicographically only.
* `dedup` is exact (no chunking), so on very large inputs it finds the few
  duplicate chains that pairtools' 10 000-row chunking misses (6 of 1.47 M on
  the 10 M-record benchmark); outputs are otherwise identical.
* `#chromosomes:` is rewritten with the sorted names only; pairtools 1.1.3
  additionally emits a stray `:` token there (a pairtools bug).
* No `.cool` writer; use `cooler load` on `bin` output.

## Benchmark methodology

See [docs/BENCHMARKS.md](docs/BENCHMARKS.md). In short: synthetic datasets
from `kira-pairs generate` (seeded, tens of millions of records, larger than
CPU caches), `/usr/bin/time -v` for wall/CPU/peak RSS, identical thread and
memory settings for both tools, outputs cross-checked for identity, and
Criterion microbenchmarks in `benches/`.

```bash
cargo bench                                 # microbenchmarks
scripts/bench_end_to_end.sh --records 10000000 --pairtools /path/to/pairtools
```

## Development

```bash
cargo fmt --check
cargo check --all-targets --all-features
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features

# pairtools 1.1.3 differential tests (needs Python; see scripts/setup_pairtools_oracle.sh)
scripts/setup_pairtools_oracle.sh .venv-pairtools
python3 scripts/compare_pairtools.py --kira target/release/kira-pairs --pairtools .venv-pairtools/bin/pairtools
PAIRTOOLS=.venv-pairtools/bin/pairtools PYTHON=.venv-pairtools/bin/python scripts/regenerate_golden.sh
```

## Acknowledgements

kira-pairs reimplements semantics defined by the
[pairtools](https://github.com/open2c/pairtools) project (Open2C) and the
[4DN `.pairs` specification](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md).
BAM/SAM parsing uses the [noodles](https://github.com/zaeleus/noodles) crates.

## License

MIT. See [LICENSE](LICENSE).
