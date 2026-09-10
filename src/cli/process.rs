//! `kira-pairs process`: the fused pipeline.

use std::io::Write;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use clap::Args;

use crate::binning::Binner;
use crate::chroms::ChromSizes;
use crate::cli::bin::{write_bins, write_table};
use crate::cli::common::{Context, Resources, append_pg, new_dict, stdio_path};
use crate::cli::dedup::{DedupOutputs, Route};
use crate::cli::parse::ParseOpts;
use crate::cli::{BinFormatArg, ClusteringArg, MethodArg, ResourceOpts};
use crate::dedup::{Clustering, DedupConfig, Method};
use crate::error::{KiraError, Result};
use crate::metrics::Metrics;
use crate::parse::{HicParser, open_alignments};
use crate::pipeline::{bin_config, consume_sorted, parse_and_sort};
use crate::sort::external::SortConfig;
use crate::stats::{DistBins, StatsAccumulator, StatsFormat};

/// Arguments for `process`.
#[derive(Debug, Args)]
pub struct ProcessArgs {
    /// Input SAM or BAM (name-grouped); `-`/omitted = stdin.
    pub input: Option<PathBuf>,
    /// Output deduplicated .pairs file (omit to skip).
    #[arg(long, value_name = "FILE")]
    pub output_pairs: Option<PathBuf>,
    /// Output duplicates (may equal --output-pairs).
    #[arg(long, value_name = "FILE")]
    pub output_dups: Option<PathBuf>,
    /// Output pairs with an unmapped side (may equal --output-pairs).
    #[arg(long, value_name = "FILE")]
    pub output_unmapped: Option<PathBuf>,
    /// Output statistics (TSV, or YAML/JSON with --yaml/--json).
    #[arg(long, value_name = "FILE")]
    pub output_stats: Option<PathBuf>,
    /// Statistics as YAML.
    #[arg(long)]
    pub yaml: bool,
    /// Statistics as JSON.
    #[arg(long, conflicts_with = "yaml")]
    pub json: bool,
    /// Output binned contacts; one path per --resolution or a path with `{res}`.
    #[arg(long, value_name = "FILE")]
    pub output_bins: Vec<PathBuf>,
    /// Bin table output(s) (`chrom start end`).
    #[arg(long, value_name = "FILE")]
    pub output_bin_table: Vec<PathBuf>,
    /// Bin resolution(s) in bp.
    #[arg(short = 'r', long, value_name = "N")]
    pub resolution: Vec<u64>,
    /// Bin output format.
    #[arg(long, value_enum, default_value_t = BinFormatArg::Coo)]
    pub bin_format: BinFormatArg,
    /// Maximal mismatch for duplicates.
    #[arg(long, default_value_t = 3, value_name = "N")]
    pub max_mismatch: u64,
    /// Duplicate metric.
    #[arg(long, value_enum, default_value_t = MethodArg::Max)]
    pub method: MethodArg,
    /// Duplicate clustering.
    #[arg(long, value_enum, default_value_t = ClusteringArg::Transitive)]
    pub clustering: ClusteringArg,
    /// Keep the original pair_type of duplicates instead of `DD`.
    #[arg(long)]
    pub no_mark_dups: bool,
    /// Add parent_readID to duplicates.
    #[arg(long)]
    pub keep_parent_id: bool,
    /// Parse options.
    #[command(flatten)]
    pub parse: ParseOpts,
    /// Resource options.
    #[command(flatten)]
    pub res: ResourceOpts,
}

/// Run `process`.
pub fn run(a: ProcessArgs, ctx: &Context) -> Result<()> {
    let res = Resources::from_opts(&a.res)?;
    let mut metrics = Metrics::start();
    if a.output_pairs.is_none() && a.output_stats.is_none() && a.output_bins.is_empty() {
        return Err(KiraError::arg(
            "nothing to do: give --output-pairs, --output-stats and/or --output-bins",
        ));
    }
    if !a.output_bins.is_empty() && a.resolution.is_empty() {
        return Err(KiraError::arg("--output-bins requires --resolution"));
    }
    if a.parse.no_flip {
        log::warn!(
            "--no-flip: pairs are not upper-triangular, duplicate detection may be incomplete"
        );
    }
    let chroms = ChromSizes::from_path(&a.parse.chroms_path)?;
    let cfg = a.parse.to_config();
    let (source, sam_header, name) = open_alignments(stdio_path(&a.input), res.io_threads)?;
    let dict = new_dict();
    let mut parser = HicParser::new(cfg, sam_header, &chroms, Arc::clone(&dict))?;
    append_pg(parser.header_mut(), "process", ctx)?;
    let mut header = parser.header().clone();
    header.mark_sorted()?;
    let cols = crate::pairs::columns::ColumnMap::from_names(header.columns())?;
    // Sorter configuration.
    let mut sort_cfg = SortConfig::new(res.threads, res.budget);
    sort_cfg.tmpdir = res.tmpdir.clone();
    sort_cfg.input_name = name.clone();
    let t0 = Instant::now();
    let (mut stream, mut pm) = parse_and_sort(source, &mut parser, &name, sort_cfg, None)?;
    // Outputs.
    let mark_dups = !a.no_mark_dups;
    let out_path = stdio_path(&a.output_pairs);
    let dups_path = stdio_path(&a.output_dups);
    let unmapped_path = stdio_path(&a.output_unmapped);
    let mut dups_header = header.clone();
    if a.keep_parent_id {
        dups_header.append_columns(&["parent_readID".to_string()]);
    }
    let dups_same = a.output_dups.is_some()
        && a.output_pairs.is_some()
        && crate::io::compression::same_output(dups_path, out_path);
    let main = match &a.output_pairs {
        Some(_) => {
            let mut w = res.open_writer(out_path)?;
            w.write_header(if dups_same { &dups_header } else { &header })?;
            Some(w)
        }
        None => None,
    };
    let dups = if a.output_dups.is_none() {
        Route::Drop
    } else if dups_same {
        Route::Main
    } else {
        let mut w = res.open_writer(dups_path)?;
        w.write_header(&dups_header)?;
        Route::Own(w)
    };
    let unmapped = if a.output_unmapped.is_none() {
        Route::Drop
    } else if a.output_pairs.is_some()
        && crate::io::compression::same_output(unmapped_path, out_path)
    {
        Route::Main
    } else if a.output_dups.is_some()
        && crate::io::compression::same_output(unmapped_path, dups_path)
    {
        Route::Dups
    } else {
        let mut w = res.open_writer(unmapped_path)?;
        w.write_header(&header)?;
        Route::Own(w)
    };
    let stats = a.output_stats.as_ref().map(|_| {
        StatsAccumulator::new(DistBins::default(), dict.get(crate::chroms::UNMAPPED_CHROM))
    });
    let mut outputs = DedupOutputs {
        main,
        dups,
        unmapped,
        mark_dups,
        keep_parent_id: a.keep_parent_id,
        parent_only_on_dups: a.clustering == ClusteringArg::Greedy,
        readid_col: cols.readid,
        pair_type_col: cols.pair_type,
        stats,
        scratch: Vec::with_capacity(256),
    };
    let mut binner = if a.output_bins.is_empty() {
        None
    } else {
        let bcfg = bin_config(
            a.resolution.clone(),
            &chroms,
            res.budget.records / 4,
            res.tmpdir.as_deref(),
        );
        Some(Binner::new(bcfg, Arc::clone(&dict))?)
    };
    let dedup_cfg = DedupConfig {
        max_mismatch: a.max_mismatch,
        method: match a.method {
            MethodArg::Max => Method::Max,
            MethodArg::Sum => Method::Sum,
        },
        clustering: match a.clustering {
            ClusteringArg::Transitive => Clustering::Transitive,
            ClusteringArg::Greedy => Clustering::Greedy,
        },
        keep_parent_id: a.keep_parent_id,
        extra_col_pairs: Vec::new(),
        readid_col: cols.readid,
        input_name: name.clone(),
    };
    let t1 = Instant::now();
    let dm = consume_sorted(
        &mut stream,
        dedup_cfg,
        &dict,
        |e| outputs.handle(e),
        binner.as_mut(),
    )?;
    drop(stream);
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
        let mut w = crate::io::compression::open_output(Some(path), res.output_options())?;
        w.write_all(text.as_bytes())?;
        w.finish()?;
    }
    if let Some(b) = binner {
        let outputs_paths = resolve_bin_paths(&a.output_bins, &a.resolution)?;
        let table_paths = if a.output_bin_table.is_empty() {
            Vec::new()
        } else {
            resolve_bin_paths(&a.output_bin_table, &a.resolution)?
        };
        let (tables, _bm) = b.finish()?;
        for (i, (layout, mut table)) in tables.into_iter().enumerate() {
            write_table(
                &mut table,
                &layout,
                &chroms,
                a.bin_format,
                outputs_paths[i].as_deref(),
                &res,
            )?;
            if let Some(p) = table_paths.get(i) {
                write_bins(&layout, &chroms, p.as_deref(), &res)?;
            }
        }
    }
    pm.dedup_seconds = t1.elapsed().as_secs_f64();
    pm.duplicates = dm.duplicates;
    if ctx.metrics {
        metrics.set("alignments_read", pm.alignments_read);
        metrics.set("records_read", pm.pairs_parsed);
        metrics.set("duplicates", pm.duplicates);
        metrics.set_f("parse_seconds", pm.parse_seconds);
        metrics.set_f("sort_run_seconds", pm.sort.run_seconds);
        metrics.set("number_of_runs", pm.sort.runs);
        metrics.set("temporary_bytes_written", pm.sort.temp_bytes);
        metrics.set("peak_records_buffered", pm.sort.peak_records_buffered);
        metrics.set_f("merge_dedup_seconds", pm.dedup_seconds);
        metrics.set_f("total_seconds", t0.elapsed().as_secs_f64());
        metrics.finalize();
        metrics.print();
    }
    Ok(())
}

fn resolve_bin_paths(paths: &[PathBuf], resolutions: &[u64]) -> Result<Vec<Option<PathBuf>>> {
    if paths.len() == 1 && paths[0].to_string_lossy().contains("{res}") {
        let t = paths[0].to_string_lossy().to_string();
        return Ok(resolutions
            .iter()
            .map(|r| Some(PathBuf::from(t.replace("{res}", &r.to_string()))))
            .collect());
    }
    if paths.len() != resolutions.len() {
        return Err(KiraError::arg(format!(
            "--output-bins: got {} path(s) for {} resolution(s) (use one path per resolution or a `{{res}}` placeholder)",
            paths.len(),
            resolutions.len()
        )));
    }
    Ok(paths.iter().map(|p| Some(p.clone())).collect())
}
