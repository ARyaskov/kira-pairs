//! `kira-pairs sort`.

use std::path::PathBuf;
use std::sync::Arc;
use std::time::Instant;

use clap::Args;

use crate::cli::ResourceOpts;
use crate::cli::common::{Context, Resources, append_pg, new_dict, stdio_path};
use crate::error::Result;
use crate::metrics::{Metrics, Progress};
use crate::pairs::reader::PairsReader;
use crate::sort::external::{ExternalSorter, SortConfig};

/// Arguments for `sort`.
#[derive(Debug, Args)]
pub struct SortArgs {
    /// Input .pairs file (plain, .gz/BGZF or .lz4); `-` or omitted = stdin.
    pub input: Option<PathBuf>,
    /// Output file (`.gz` -> BGZF, `.lz4` -> LZ4); `-` or omitted = stdout.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// Extra column (name or 0-based index) used as an additional sort key
    /// after pair_type. May be repeated.
    #[arg(long = "extra-col", value_name = "COLUMN")]
    pub extra_col: Vec<String>,
    /// Maximum number of runs merged at once.
    #[arg(long, default_value_t = 64, value_name = "N")]
    pub max_fan_in: usize,
    /// Do not LZ4-compress temporary runs.
    #[arg(long)]
    pub no_compress_runs: bool,
    /// Resource options.
    #[command(flatten)]
    pub res: ResourceOpts,
}

/// Run `sort`.
pub fn run(a: SortArgs, ctx: &Context) -> Result<()> {
    let res = Resources::from_opts(&a.res)?;
    let mut metrics = Metrics::start();
    let input = stdio_path(&a.input);
    let reader = PairsReader::open(input, res.io_threads)?;
    let (mut header, cols, body) = reader.into_parts();
    let name = body.name().to_string();
    append_pg(&mut header, "sort", ctx)?;
    if header.is_valid_pairs_header() {
        header.mark_sorted()?;
    } else if !header.is_empty() {
        header.set_field("sorted", crate::pairs::header::SORTED_VALUE);
    }
    let extra_cols: Vec<usize> = a
        .extra_col
        .iter()
        .map(|c| cols.require(c))
        .collect::<Result<_>>()?;
    if cols.pair_type.is_none() {
        log::warn!("no pair_type column: sorting by chromosomes and positions only");
    }
    let mut cfg = SortConfig::new(res.threads, res.budget);
    cfg.tmpdir = res.tmpdir.clone();
    cfg.compress_runs = !a.no_compress_runs;
    cfg.max_fan_in = a.max_fan_in.max(2);
    cfg.extra_cols = extra_cols;
    cfg.input_name = name;
    let dict = new_dict();
    let mut sorter = ExternalSorter::new(cfg, cols, Arc::clone(&dict))?;
    let mut writer = res.open_writer(stdio_path(&a.output))?;
    writer.write_header(&header)?;

    let mut progress = Progress::new(ctx.progress, std::time::Duration::from_secs(5));
    let t0 = Instant::now();
    let mut body = body;
    let mut lines_in = 0u64;
    while let Some(block) = body.next_block()? {
        lines_in += block.n_lines;
        progress.tick(lines_in, body.bytes_read());
        sorter.push_block(block)?;
    }
    let bytes_read = body.bytes_read();
    let (mut stream, sm) = sorter.finish()?;
    let t1 = Instant::now();
    let mut written = 0u64;
    stream.for_each(|_k, line| {
        written += 1;
        writer.write_line(line)
    })?;
    let bytes_written = writer.finish()?;
    progress.finish(written, bytes_read);
    if ctx.metrics {
        metrics.set("records_read", sm.records);
        metrics.set("records_written", written);
        metrics.set("bytes_read", bytes_read);
        metrics.set("bytes_written", bytes_written);
        metrics.set_f("parse_and_run_seconds", t1.duration_since(t0).as_secs_f64());
        metrics.set_f("sort_run_seconds", sm.run_seconds);
        metrics.set_f("merge_seconds", t1.elapsed().as_secs_f64());
        metrics.set_f("merge_pass_seconds", sm.merge_pass_seconds);
        metrics.set("merge_passes", sm.merge_passes);
        metrics.set("number_of_runs", sm.runs);
        metrics.set("temporary_bytes_written", sm.temp_bytes);
        metrics.set("peak_records_buffered", sm.peak_records_buffered);
        metrics.finalize();
        metrics.print();
    }
    Ok(())
}
