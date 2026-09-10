//! `kira-pairs dedup`.

use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use clap::Args;

use crate::cli::common::{Context, Resources, append_pg, new_dict, stdio_path};
use crate::cli::{BackendArg, ClusteringArg, MethodArg, ResourceOpts};
use crate::dedup::{Clustering, DedupConfig, Deduper, Emitted, Method, Outcome, mark_dd};
use crate::error::{KiraError, Result};
use crate::io::compression::same_output;
use crate::metrics::{Metrics, Progress};
use crate::pairs::parallel::OrderedParser;
use crate::pairs::reader::PairsReader;
use crate::pairs::writer::PairsWriter;
use crate::stats::{DistBins, StatsAccumulator, StatsFormat};

/// Arguments for `dedup`.
#[derive(Debug, Args)]
pub struct DedupArgs {
    /// Input sorted, upper-triangular .pairs file; `-` or omitted = stdin.
    pub input: Option<PathBuf>,
    /// Output file for deduplicated pairs; `-` or omitted = stdout.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// Output file for duplicates. Same path as --output (or `-`) writes
    /// them together with deduplicated pairs. By default duplicates are dropped.
    #[arg(long, value_name = "FILE")]
    pub output_dups: Option<PathBuf>,
    /// Output file for pairs with an unmapped side (may equal --output or
    /// --output-dups). By default unmapped pairs are dropped.
    #[arg(long, value_name = "FILE")]
    pub output_unmapped: Option<PathBuf>,
    /// Output file for statistics.
    #[arg(long, value_name = "FILE")]
    pub output_stats: Option<PathBuf>,
    /// Write statistics as YAML instead of TSV.
    #[arg(long)]
    pub yaml: bool,
    /// Write statistics as JSON (kira-pairs extension).
    #[arg(long, conflicts_with = "yaml")]
    pub json: bool,
    /// Distance bins per decade for the stats distance histogram.
    #[arg(long, default_value_t = 8, value_name = "N")]
    pub n_dist_bins_decade: usize,
    /// Pairs with both sides within this distance (bp) are duplicates.
    #[arg(long, default_value_t = 3, value_name = "N")]
    pub max_mismatch: u64,
    /// Distance metric.
    #[arg(long, value_enum, default_value_t = MethodArg::Max)]
    pub method: MethodArg,
    /// Clustering mode (transitive = pairtools scipy backend).
    #[arg(long, value_enum, default_value_t = ClusteringArg::Transitive)]
    pub clustering: ClusteringArg,
    /// pairtools backend name (alias for --clustering).
    #[arg(long, value_enum, conflicts_with = "clustering")]
    pub backend: Option<BackendArg>,
    /// Mark duplicates as `DD` in pair_type.
    #[arg(long, default_value_t = true, overrides_with = "no_mark_dups")]
    pub mark_dups: bool,
    /// Keep the original pair_type of duplicates.
    #[arg(long)]
    pub no_mark_dups: bool,
    /// Add a `parent_readID` column to duplicates.
    #[arg(long)]
    pub keep_parent_id: bool,
    /// Extra column pair (names or 0-based indices) that must match between
    /// duplicates, e.g. --extra-col-pair phase1 phase2. May be repeated.
    #[arg(long = "extra-col-pair", num_args = 2, value_names = ["COLUMN1", "COLUMN2"])]
    pub extra_col_pair: Vec<String>,
    /// Which outputs receive the header.
    #[arg(long, default_value = "both", value_parser = ["dups", "dedup", "both", "none"])]
    pub send_header_to: String,
    /// Resource options.
    #[command(flatten)]
    pub res: ResourceOpts,
}

/// Output routing for dedup results.
pub struct DedupOutputs {
    /// Deduplicated pairs (`None` = dropped).
    pub main: Option<PairsWriter>,
    /// Duplicates: `Main` means "write into the main output".
    pub dups: Route,
    /// Unmapped pairs.
    pub unmapped: Route,
    /// Whether to mark duplicates as `DD`.
    pub mark_dups: bool,
    /// Whether to append the parent read ID to duplicates.
    pub keep_parent_id: bool,
    /// pairtools cython-backend semantics: only duplicates carry a parent id.
    /// Otherwise (scipy backend) kept pairs carry their own read id when
    /// duplicates share the output stream and unmapped pairs an empty one.
    pub parent_only_on_dups: bool,
    /// Column index of `readID`.
    pub readid_col: usize,
    /// Pair type column (for `DD` marking).
    pub pair_type_col: Option<usize>,
    /// Statistics accumulator.
    pub stats: Option<StatsAccumulator>,
    /// Scratch buffer for rewritten lines.
    pub scratch: Vec<u8>,
}

/// Where a class of records goes.
pub enum Route {
    /// Dropped.
    Drop,
    /// Written to the main output.
    Main,
    /// Written to the duplicates output.
    Dups,
    /// Written to a dedicated writer.
    Own(PairsWriter),
}

impl DedupOutputs {
    /// Route one emitted record.
    pub fn handle(&mut self, e: Emitted<'_>) -> Result<()> {
        match &e.outcome {
            Outcome::Unmapped => {
                if let Some(s) = self.stats.as_mut() {
                    s.observe(e.key, None, false);
                }
                let line: &[u8] = if self.keep_parent_id && !self.parent_only_on_dups {
                    self.scratch.clear();
                    self.scratch.extend_from_slice(e.line);
                    self.scratch.push(b'\t');
                    &self.scratch
                } else {
                    e.line
                };
                match &mut self.unmapped {
                    Route::Drop => {}
                    Route::Main => {
                        if let Some(w) = self.main.as_mut() {
                            w.write_line(line)?;
                        }
                    }
                    Route::Dups => {
                        if let Route::Own(w) = &mut self.dups {
                            w.write_line(line)?;
                        } else if let Some(w) = self.main.as_mut() {
                            w.write_line(line)?;
                        }
                    }
                    Route::Own(w) => w.write_line(line)?,
                }
            }
            Outcome::Unique => {
                if let Some(s) = self.stats.as_mut() {
                    s.observe(e.key, None, false);
                }
                let shared = matches!(self.dups, Route::Main);
                if self.keep_parent_id && !self.parent_only_on_dups && shared {
                    let mut ends = Vec::with_capacity(16);
                    crate::pairs::record::split_fields(e.line, &mut ends);
                    let rec = crate::pairs::record::PairRecordRef::new(e.line, &ends);
                    self.scratch.clear();
                    self.scratch.extend_from_slice(e.line);
                    self.scratch.push(b'\t');
                    self.scratch
                        .extend_from_slice(rec.field(self.readid_col).unwrap_or(b""));
                    if let Some(w) = self.main.as_mut() {
                        w.write_line(&self.scratch)?;
                    }
                } else if let Some(w) = self.main.as_mut() {
                    w.write_line(e.line)?;
                }
            }
            Outcome::Duplicate { parent } => {
                if let Some(s) = self.stats.as_mut() {
                    s.observe(e.key, if self.mark_dups { Some(b"DD") } else { None }, true);
                }
                let writer: Option<&mut PairsWriter> = match &mut self.dups {
                    Route::Drop => None,
                    Route::Main | Route::Dups => self.main.as_mut(),
                    Route::Own(w) => Some(w),
                };
                if let Some(w) = writer {
                    let scratch = &mut self.scratch;
                    if self.mark_dups {
                        mark_dd(e.line, self.pair_type_col, scratch);
                    } else {
                        scratch.clear();
                        scratch.extend_from_slice(e.line);
                    }
                    if self.keep_parent_id {
                        scratch.push(b'\t');
                        scratch.extend_from_slice(parent);
                    }
                    w.write_line(scratch)?;
                }
            }
        }
        Ok(())
    }

    /// Finish all writers.
    pub fn finish(self) -> Result<Option<StatsAccumulator>> {
        if let Some(w) = self.main {
            w.finish()?;
        }
        if let Route::Own(w) = self.dups {
            w.finish()?;
        }
        if let Route::Own(w) = self.unmapped {
            w.finish()?;
        }
        Ok(self.stats)
    }
}

/// Resolve `--extra-col-pair` arguments to column index pairs.
pub fn resolve_extra_col_pairs(
    args: &[String],
    cols: &crate::pairs::columns::ColumnMap,
) -> Result<Vec<(usize, usize)>> {
    let mut out = Vec::new();
    for pair in args.chunks(2) {
        if pair.len() != 2 {
            return Err(KiraError::arg("--extra-col-pair needs two columns"));
        }
        out.push((cols.require(&pair[0])?, cols.require(&pair[1])?));
    }
    Ok(out)
}

/// Run `dedup`.
pub fn run(a: DedupArgs, ctx: &Context) -> Result<()> {
    let res = Resources::from_opts(&a.res)?;
    let mut metrics = Metrics::start();
    let mark_dups = a.mark_dups && !a.no_mark_dups;
    let clustering = match a.backend {
        Some(BackendArg::Cython) => Clustering::Greedy,
        Some(_) => Clustering::Transitive,
        None => match a.clustering {
            ClusteringArg::Transitive => Clustering::Transitive,
            ClusteringArg::Greedy => Clustering::Greedy,
        },
    };
    let input = stdio_path(&a.input);
    let reader = PairsReader::open(input, res.io_threads)?;
    let (mut header, cols, body) = reader.into_parts();
    let name = body.name().to_string();
    if !header.is_sorted() && !header.is_empty() {
        log::warn!(
            "pairs file appears not to be sorted (no #sorted: header); dedup might produce wrong results"
        );
    }
    append_pg(&mut header, "dedup", ctx)?;
    let mut dups_header = header.clone();
    if a.keep_parent_id && !dups_header.is_empty() {
        dups_header.append_columns(&["parent_readID".to_string()]);
    }
    let out_path = stdio_path(&a.output);
    let dups_path = stdio_path(&a.output_dups);
    let unmapped_path = stdio_path(&a.output_unmapped);
    let dups_same_as_main =
        a.output_dups.is_some() && (dups_path.is_none() || same_output(dups_path, out_path));
    if dups_same_as_main {
        header = dups_header.clone();
    }
    let send_main = matches!(a.send_header_to.as_str(), "both" | "dedup");
    let send_dups = matches!(a.send_header_to.as_str(), "both" | "dups");

    let mut main = res.open_writer(out_path)?;
    if send_main {
        main.write_header(&header)?;
    }
    let dups = if a.output_dups.is_none() {
        Route::Drop
    } else if dups_same_as_main {
        Route::Main
    } else {
        let mut w = res.open_writer(dups_path)?;
        if send_dups {
            w.write_header(&dups_header)?;
        }
        Route::Own(w)
    };
    let unmapped = if a.output_unmapped.is_none() {
        Route::Drop
    } else if unmapped_path.is_none() || same_output(unmapped_path, out_path) {
        Route::Main
    } else if a.output_dups.is_some() && same_output(unmapped_path, dups_path) {
        Route::Dups
    } else {
        let mut w = res.open_writer(unmapped_path)?;
        w.write_header(&header)?;
        Route::Own(w)
    };
    let dict = new_dict();
    let stats = a.output_stats.as_ref().map(|_| {
        StatsAccumulator::new(
            DistBins::new(a.n_dist_bins_decade),
            dict.get(crate::chroms::UNMAPPED_CHROM),
        )
    });
    let mut outputs = DedupOutputs {
        main: Some(main),
        dups,
        unmapped,
        mark_dups,
        keep_parent_id: a.keep_parent_id,
        parent_only_on_dups: clustering == Clustering::Greedy,
        readid_col: cols.readid,
        pair_type_col: cols.pair_type,
        stats,
        scratch: Vec::with_capacity(256),
    };
    let cfg = DedupConfig {
        max_mismatch: a.max_mismatch,
        method: match a.method {
            MethodArg::Max => Method::Max,
            MethodArg::Sum => Method::Sum,
        },
        clustering,
        keep_parent_id: a.keep_parent_id,
        extra_col_pairs: resolve_extra_col_pairs(&a.extra_col_pair, &cols)?,
        readid_col: cols.readid,
        input_name: name,
    };
    let mut deduper = Deduper::new(cfg, &dict);
    let cols = Arc::new(cols);
    let mut parser = OrderedParser::with_depth(
        body,
        Arc::clone(&cols),
        Arc::clone(&dict),
        res.threads,
        res.budget.channel_depth(res.threads),
    );
    let mut progress = Progress::new(ctx.progress, std::time::Duration::from_secs(5));
    let t0 = Instant::now();
    let mut records = 0u64;
    let mut sink = |e: Emitted<'_>| outputs.handle(e);
    while let Some(chunk) = parser.next_chunk()? {
        for e in &chunk.entries {
            records += 1;
            deduper.push(&e.key, chunk.line(e), e.key.seq, &mut sink)?;
        }
        progress.tick(records, parser.bytes_read());
    }
    deduper.finish(&mut sink)?;
    let dm = deduper.metrics().clone();
    let bytes_read = parser.bytes_read();
    drop(parser);
    let stats = outputs.finish()?;
    if let (Some(path), Some(s)) = (&a.output_stats, stats) {
        let fmt = if a.yaml {
            StatsFormat::Yaml
        } else if a.json {
            StatsFormat::Json
        } else {
            StatsFormat::Tsv
        };
        let text = crate::stats::format::render(&s.snapshot(&dict), fmt);
        let mut w = crate::io::compression::open_output(
            stdio_path(&Some(path.clone())),
            res.output_options(),
        )?;
        std::io::Write::write_all(&mut w, text.as_bytes())?;
        w.finish()?;
    }
    progress.finish(records, bytes_read);
    if ctx.metrics {
        metrics.set("records_read", dm.records);
        metrics.set("records_unmapped", dm.unmapped);
        metrics.set("duplicates", dm.duplicates);
        metrics.set("bytes_read", bytes_read);
        metrics.set("peak_pending_records", dm.peak_pending);
        metrics.set("peak_window_records", dm.peak_window);
        metrics.set_f("dedup_seconds", t0.elapsed().as_secs_f64());
        metrics.finalize();
        metrics.print();
    }
    Ok(())
}
