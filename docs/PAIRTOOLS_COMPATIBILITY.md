# pairtools 1.1.3 compatibility

Reference: pairtools 1.1.3 (source inspected; installed as an executable
oracle via `scripts/setup_pairtools_oracle.sh`). Differential tests live in
`scripts/compare_pairtools.py` (runs both tools on fixtures and a 300k-record
synthetic dataset) and in `tests/cli.rs` (compares against committed pairtools
outputs in `tests/golden/`). "exact" means byte-identical output bodies and
headers after removing each tool's own `@PG` provenance record.

| Feature                                   | Status | Evidence |
| ----------------------------------------- | ------ | -------- |
| `.pairs` reading (header, columns, extra columns) | exact | all tests |
| `.pairs` writing                          | exact  | all tests |
| gzip / BGZF input                          | exact  | `stdin_stdout_pipeline_and_gzip` |
| LZ4 input/output                           | exact  | `stdin_stdout_pipeline_and_gzip` |
| stdin/stdout pipelines                     | exact  | compare script, cli tests |
| header preservation (unknown fields, comments) | exact | `edge.pairs` |
| `@PG` provenance record                    | equivalent (own ID/PN/VN) | — |
| `#chromosomes:` after sort                 | deviates: names only, pairtools adds a stray `:` token (bug) | `edge.pairs` |
| sort order (chrom1, chrom2, pos1, pos2, pair_type, stable) | exact | sort tests, 300k synthetic |
| sort `--extra-col`                         | lexicographic only (pairtools: numeric for int dtypes) | unit test |
| `#sorted:` header update                   | exact  | sort tests |
| flip                                       | exact  | `flip_matches_pairtools`, compare script |
| dedup `--method max`                       | exact  | adversarial + synthetic |
| dedup `--method sum`                       | exact  | adversarial + synthetic |
| dedup transitive clustering (scipy/sklearn backends) | exact (see note 1) | adversarial |
| dedup greedy clustering (cython backend)   | exact for classification; parent ids correct in kira-pairs, wrong in pairtools (note 2) | adversarial |
| dedup `--keep-parent-id`                   | exact (transitive)  | goldens |
| dedup `--mark-dups/--no-mark-dups`         | exact  | compare script |
| dedup `--output-dups`, `--output-unmapped` (incl. shared streams) | exact | goldens, compare script |
| dedup `--extra-col-pair`                   | equality on both columns (note 3) | unit test |
| dedup `--output-stats`                     | exact (note 4) | goldens |
| dedup on unsorted input                    | kira-pairs errors; pairtools silently miscounts | `dedup_rejects_unsorted_input` |
| stats core counters, `pair_types`, `chrom_freq`, `dist_freq`, `cis_Nkb+` | exact | goldens |
| stats `summary/*` (fractions, `complexity_naive`, convergence) | exact to 1e-9 relative (note 4) | goldens |
| stats `--yaml`                             | equivalent (same values; pairtools omits zero entries, so do we) | `stats_matches_pairtools_and_merges` |
| stats `--merge`                            | exact  | compare script |
| stats `--n-dist-bins-decade`               | implemented, same bin edges | unit test |
| stats `--filter`, by-tile stats            | not implemented | — |
| select operators (`== != < <= > >= and or not in`), arithmetic, `abs`, `csv_match`, `wildcard_match`, `regex_match`, `COLS[i]` | exact for these forms | goldens, compare script |
| select arbitrary Python                    | not supported (safe subset only) | — |
| select `--output-rest`, `--chrom-subset`, `--remove-columns`, `--type-cast` | implemented | tests |
| bin                                        | kira-pairs extension (cooler bin numbering, COO/BG2) | tests |
| parse `--walks-policy 5unique/5any/3unique/3any/mask` | exact | goldens, compare script |
| parse `--walks-policy all`, `parse2`       | not implemented (error) | — |
| parse `--min-mapq`, `--max-molecule-size`, `--max-inter-align-gap`, `--report-alignment-end`, `--no-flip` | exact | compare script |
| parse pairsam output (`sam1`/`sam2`)       | exact  | `hic.parse.pairsam` |
| parse `--drop-sam`, `--drop-readid`, `--add-pair-index`, `--assembly` | exact | compare script |
| parse `--add-columns` (mapq, pos5, pos3, cigar, read_len, matched_bp, algn_ref_span, algn_read_span, dist_to_5, dist_to_3, read_side, algn_idx, same_side_algn_count, SAM tags) | exact | goldens |
| parse `--add-columns seq`, `--drop-seq`    | implemented; oracle crashes (pysam/Python 3.14), verified by unit tests only | — |
| parse `--add-columns mismatches`, `--readid-transform`, `--output-parsed-alignments` | not implemented | — |
| parse header (`#chromsize` order, `#samheader`) | exact | goldens |
| parse CRAM                                 | not supported | — |
| process (parse→sort→dedup→stats)           | exact vs `pairtools parse \| sort \| dedup` | `process_matches_parse_sort_dedup` |
| phase, restrict, filterbycov, scaling, sample, split, merge, header | not implemented | — |

Notes

1. pairtools' scipy backend processes 10 000-row chunks with a 100-row
   carry-over; a duplicate cluster straddling a chunk boundary can be missed
   by pairtools ("the algorithm might miss a few duplicates", pairtools
   docs). kira-pairs uses an exact sweep with no chunking, so on very large
   inputs pairtools may report slightly *fewer* duplicates than kira-pairs.
   Differential tests use inputs where the effect does not occur.
2. pairtools' cython backend stores parent indices relative to an internal
   buffer that is shifted without rebasing the indices, so `parent_readID`
   values are wrong after the first buffer shrink (observable on a 30-line
   input). kira-pairs reports the true parent. Goldens for the greedy mode
   therefore cover classification only.
3. pairtools' scipy backend has a quirk for column pairs with differing names
   (it intersects arbitrary group labels); the cython backend and kira-pairs
   require both columns to match between duplicates, which is the documented
   intent.
4. Floating point summaries are computed with the same formulas
   (`scipy.special.lambertw` replaced by a Halley iteration) and printed with
   Python's `repr` formatting; values agree to ≥ 1e-12 relative and are
   compared with a tolerance in tests. `pair_types/*` and `chrom_freq/*`
   lines are emitted in sorted order (pairtools: first-seen order).

Known environment note: the oracle was run on Python 3.14 with pandas 2.3
and a `pipes` shim; pairtools' own behaviour is unaffected by these.
