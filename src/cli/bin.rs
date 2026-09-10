//! `kira-pairs bin`.

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use clap::Args;

use crate::binning::aggregate::BinTable;
use crate::binning::{BinConfig, BinLayout, Binner};
use crate::chroms::ChromSizes;
use crate::cli::common::{Context, Resources, new_dict, parse_pair_types, stdio_path};
use crate::cli::{BinFormatArg, ResourceOpts};
use crate::error::{KiraError, Result};
use crate::metrics::{Metrics, Progress};
use crate::pairs::parallel::OrderedParser;
use crate::pairs::reader::PairsReader;
use crate::pairs::record::split_fields;
use crate::util::int::{parse_u64, write_u64};

/// Arguments for `bin`.
#[derive(Debug, Args)]
pub struct BinArgs {
    /// Input .pairs file; `-` or omitted = stdin.
    pub input: Option<PathBuf>,
    /// Output file(s). One per --resolution (in order), or a single path
    /// containing `{res}`; `-` or omitted = stdout (single resolution only).
    #[arg(short, long, value_name = "FILE")]
    pub output: Vec<PathBuf>,
    /// Chromosome sizes (defines the bin table).
    #[arg(short = 'c', long, value_name = "FILE")]
    pub chroms_path: PathBuf,
    /// Bin size in bp. May be repeated to produce several tables in one pass.
    #[arg(short = 'r', long, value_name = "N", required = true)]
    pub resolution: Vec<u64>,
    /// Minimal MAPQ on both sides (requires mapq1/mapq2 columns).
    #[arg(long, value_name = "N")]
    pub min_mapq: Option<u64>,
    /// Comma-separated pair types to accept (default: all mapped pairs).
    #[arg(long, value_name = "LIST")]
    pub pair_types: Option<String>,
    /// Output format.
    #[arg(long, value_enum, default_value_t = BinFormatArg::Coo)]
    pub format: BinFormatArg,
    /// Treat positions as 0-based instead of 1-based.
    #[arg(long)]
    pub zero_based: bool,
    /// Also write the bin table(s) (`chrom start end`), one per resolution
    /// or a path with `{res}`.
    #[arg(long, value_name = "FILE")]
    pub bins_out: Vec<PathBuf>,
    /// Resource options.
    #[command(flatten)]
    pub res: ResourceOpts,
}

fn resolve_paths(
    paths: &[PathBuf],
    resolutions: &[u64],
    what: &str,
) -> Result<Vec<Option<PathBuf>>> {
    if paths.is_empty() {
        if resolutions.len() == 1 {
            return Ok(vec![None]);
        }
        return Err(KiraError::arg(format!(
            "{what}: several resolutions need one path per resolution or a path with {{res}}"
        )));
    }
    if paths.len() == 1 && paths[0].to_string_lossy().contains("{res}") {
        let t = paths[0].to_string_lossy().to_string();
        return Ok(resolutions
            .iter()
            .map(|r| Some(PathBuf::from(t.replace("{res}", &r.to_string()))))
            .collect());
    }
    if paths.len() != resolutions.len() {
        return Err(KiraError::arg(format!(
            "{what}: got {} path(s) for {} resolution(s)",
            paths.len(),
            resolutions.len()
        )));
    }
    Ok(paths
        .iter()
        .map(|p| {
            if p.as_os_str() == "-" {
                None
            } else {
                Some(p.clone())
            }
        })
        .collect())
}

/// Write one table.
pub fn write_table(
    table: &mut BinTable,
    layout: &BinLayout,
    cs: &ChromSizes,
    format: BinFormatArg,
    path: Option<&Path>,
    res: &Resources,
) -> Result<u64> {
    let mut w = crate::io::compression::open_output(path, res.output_options())?;
    let mut buf = Vec::with_capacity(1 << 20);
    let mut rows = 0u64;
    while let Some((b1, b2, n)) = table.next_row()? {
        match format {
            BinFormatArg::Coo => {
                write_u64(&mut buf, b1);
                buf.push(b'\t');
                write_u64(&mut buf, b2);
                buf.push(b'\t');
                write_u64(&mut buf, n);
                buf.push(b'\n');
            }
            BinFormatArg::Bg2 => {
                let (c1, s1, e1) = layout
                    .locate(b1, cs)
                    .ok_or_else(|| KiraError::format("bin id out of range", Default::default()))?;
                let (c2, s2, e2) = layout
                    .locate(b2, cs)
                    .ok_or_else(|| KiraError::format("bin id out of range", Default::default()))?;
                buf.extend_from_slice(cs.names()[c1].as_bytes());
                buf.push(b'\t');
                write_u64(&mut buf, s1);
                buf.push(b'\t');
                write_u64(&mut buf, e1);
                buf.push(b'\t');
                buf.extend_from_slice(cs.names()[c2].as_bytes());
                buf.push(b'\t');
                write_u64(&mut buf, s2);
                buf.push(b'\t');
                write_u64(&mut buf, e2);
                buf.push(b'\t');
                write_u64(&mut buf, n);
                buf.push(b'\n');
            }
        }
        rows += 1;
        if buf.len() >= 1 << 20 {
            w.write_all(&buf)?;
            buf.clear();
        }
    }
    w.write_all(&buf)?;
    w.finish()?;
    Ok(rows)
}

/// Write a bin table (`chrom start end`).
pub fn write_bins(
    layout: &BinLayout,
    cs: &ChromSizes,
    path: Option<&Path>,
    res: &Resources,
) -> Result<()> {
    let mut w = crate::io::compression::open_output(path, res.output_options())?;
    let mut buf = Vec::with_capacity(1 << 20);
    for (ci, (name, size)) in cs.iter().enumerate() {
        let n = layout.chrom_bins(ci);
        for b in 0..n {
            let start = b * layout.resolution;
            let end = (start + layout.resolution).min(size);
            buf.extend_from_slice(name.as_bytes());
            buf.push(b'\t');
            write_u64(&mut buf, start);
            buf.push(b'\t');
            write_u64(&mut buf, end);
            buf.push(b'\n');
            if buf.len() >= 1 << 20 {
                w.write_all(&buf)?;
                buf.clear();
            }
        }
    }
    w.write_all(&buf)?;
    w.finish()?;
    Ok(())
}

/// Run `bin`.
pub fn run(a: BinArgs, ctx: &Context) -> Result<()> {
    let res = Resources::from_opts(&a.res)?;
    let mut metrics = Metrics::start();
    let cs = ChromSizes::from_path(&a.chroms_path)?;
    let outputs = resolve_paths(&a.output, &a.resolution, "--output")?;
    let bins_out = if a.bins_out.is_empty() {
        Vec::new()
    } else {
        resolve_paths(&a.bins_out, &a.resolution, "--bins-out")?
    };
    let reader = PairsReader::open(stdio_path(&a.input), res.io_threads)?;
    let (_header, cols, body) = reader.into_parts();
    let mapq_cols = match (cols.index_of("mapq1"), cols.index_of("mapq2")) {
        (Some(x), Some(y)) => Some((x, y)),
        _ => None,
    };
    if a.min_mapq.is_some() && mapq_cols.is_none() {
        return Err(KiraError::arg(
            "--min-mapq requires mapq1 and mapq2 columns in the input",
        ));
    }
    let dict = new_dict();
    let cfg = BinConfig {
        resolutions: a.resolution.clone(),
        chromsizes: cs.clone(),
        min_mapq: a.min_mapq,
        pair_types: a.pair_types.as_deref().map(parse_pair_types),
        zero_based: a.zero_based,
        memory_bytes: res.budget.records,
        tmpdir: res.tmpdir.clone(),
    };
    let mut binner = Binner::new(cfg, Arc::clone(&dict))?;
    let mut parser = OrderedParser::new(body, Arc::new(cols), Arc::clone(&dict), res.threads);
    let mut progress = Progress::new(ctx.progress, std::time::Duration::from_secs(5));
    let mut ends = Vec::with_capacity(32);
    let mut records = 0u64;
    let need_mapq = a.min_mapq.is_some();
    while let Some(chunk) = parser.next_chunk()? {
        for e in &chunk.entries {
            let mapq = if need_mapq {
                let line = chunk.line(e);
                split_fields(line, &mut ends);
                let (m1, m2) = mapq_cols.unwrap_or((0, 0));
                let rec = crate::pairs::record::PairRecordRef::new(line, &ends);
                let a1 = rec.field(m1).and_then(parse_u64);
                let a2 = rec.field(m2).and_then(parse_u64);
                match (a1, a2) {
                    (Some(x), Some(y)) => Some((x, y)),
                    _ => None,
                }
            } else {
                None
            };
            binner.observe(&e.key, mapq)?;
        }
        records += chunk.entries.len() as u64;
        progress.tick(records, parser.bytes_read());
    }
    let bytes_read = parser.bytes_read();
    drop(parser);
    let (tables, bm) = binner.finish()?;
    let mut rows_total = 0u64;
    for (i, (layout, mut table)) in tables.into_iter().enumerate() {
        rows_total += write_table(
            &mut table,
            &layout,
            &cs,
            a.format,
            outputs[i].as_deref(),
            &res,
        )?;
        if let Some(b) = bins_out.get(i) {
            write_bins(&layout, &cs, b.as_deref(), &res)?;
        }
    }
    progress.finish(records, bytes_read);
    if ctx.metrics {
        metrics.set("records_read", records);
        metrics.set("pairs_binned", bm.accepted);
        metrics.set("pairs_rejected", bm.rejected);
        metrics.set("pairs_unknown_chrom", bm.unknown_chrom);
        metrics.set("bin_rows_written", rows_total);
        metrics.set("bytes_read", bytes_read);
        metrics.finalize();
        metrics.print();
    }
    Ok(())
}
