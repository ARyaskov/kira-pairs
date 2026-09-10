//! `kira-pairs stats`.

use std::io::{Read, Write};
use std::path::PathBuf;
use std::sync::Arc;

use clap::Args;

use crate::cli::ResourceOpts;
use crate::cli::common::{Context, Resources, new_dict, stdio_path};
use crate::error::{KiraError, Result};
use crate::metrics::{Metrics, Progress};
use crate::pairs::parallel::OrderedParser;
use crate::pairs::reader::PairsReader;
use crate::stats::format::{parse_tsv, parse_yaml, render};
use crate::stats::{DistBins, StatsAccumulator, StatsFormat, StatsSnapshot};

/// Arguments for `stats`.
#[derive(Debug, Args)]
pub struct StatsArgs {
    /// Input .pairs file (or, with --merge, stats files); `-`/omitted = stdin.
    pub input: Vec<PathBuf>,
    /// Output file; `-` or omitted = stdout.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// Merge stats files given as inputs instead of computing statistics.
    #[arg(long)]
    pub merge: bool,
    /// YAML output (and YAML inputs for --merge).
    #[arg(long)]
    pub yaml: bool,
    /// JSON output (kira-pairs extension).
    #[arg(long, conflicts_with = "yaml")]
    pub json: bool,
    /// Distance bins per decade.
    #[arg(long, default_value_t = 8, value_name = "N")]
    pub n_dist_bins_decade: usize,
    /// Do not add `chromsizes/*` from the header.
    #[arg(long)]
    pub no_chromsizes: bool,
    /// Resource options.
    #[command(flatten)]
    pub res: ResourceOpts,
}

fn format_of(a: &StatsArgs) -> StatsFormat {
    if a.yaml {
        StatsFormat::Yaml
    } else if a.json {
        StatsFormat::Json
    } else {
        StatsFormat::Tsv
    }
}

fn write_output(a: &StatsArgs, res: &Resources, text: &str) -> Result<()> {
    let mut w = crate::io::compression::open_output(stdio_path(&a.output), res.output_options())?;
    w.write_all(text.as_bytes())?;
    w.finish()?;
    Ok(())
}

/// Run `stats`.
pub fn run(a: StatsArgs, ctx: &Context) -> Result<()> {
    let res = Resources::from_opts(&a.res)?;
    if a.merge {
        if a.input.is_empty() {
            return Err(KiraError::arg("--merge requires at least one stats file"));
        }
        let mut merged: Option<StatsSnapshot> = None;
        for p in &a.input {
            let mut src =
                crate::io::compression::open_input(stdio_path(&Some(p.clone())), res.io_threads)?;
            let mut text = String::new();
            src.reader
                .read_to_string(&mut text)
                .map_err(|e| KiraError::io(p, e))?;
            let name = p.display().to_string();
            let snap = if a.yaml {
                parse_yaml(&text, &name)?
            } else {
                parse_tsv(&text, &name)?
            };
            match merged.as_mut() {
                None => merged = Some(snap),
                Some(m) => m
                    .merge(&snap)
                    .map_err(|msg| KiraError::format(msg, crate::error::Location::file(&name)))?,
            }
        }
        let merged =
            merged.unwrap_or_else(|| StatsSnapshot::empty(DistBins::new(a.n_dist_bins_decade)));
        return write_output(&a, &res, &render(&merged, format_of(&a)));
    }
    if a.input.len() > 1 {
        return Err(KiraError::arg(
            "stats takes a single input file (use --merge for several stats files)",
        ));
    }
    let mut metrics = Metrics::start();
    let input = a.input.first().cloned();
    let reader = PairsReader::open(stdio_path(&input), res.io_threads)?;
    let (header, cols, body) = reader.into_parts();
    let dict = new_dict();
    let mut acc = StatsAccumulator::new(
        DistBins::new(a.n_dist_bins_decade),
        dict.get(crate::chroms::UNMAPPED_CHROM),
    );
    if !a.no_chromsizes {
        acc.set_chromsizes(header.chromsizes()?);
    }
    let mut parser = OrderedParser::with_depth(
        body,
        Arc::new(cols),
        Arc::clone(&dict),
        res.threads,
        res.budget.channel_depth(res.threads),
    );
    let mut progress = Progress::new(ctx.progress, std::time::Duration::from_secs(5));
    let mut records = 0u64;
    while let Some(chunk) = parser.next_chunk()? {
        for e in &chunk.entries {
            acc.observe_plain(&e.key);
        }
        records += chunk.entries.len() as u64;
        progress.tick(records, parser.bytes_read());
    }
    let bytes_read = parser.bytes_read();
    drop(parser);
    let text = render(&acc.snapshot(&dict), format_of(&a));
    write_output(&a, &res, &text)?;
    progress.finish(records, bytes_read);
    if ctx.metrics {
        metrics.set("records_read", records);
        metrics.set("bytes_read", bytes_read);
        metrics.finalize();
        metrics.print();
    }
    Ok(())
}
