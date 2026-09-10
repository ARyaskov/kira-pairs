//! `kira-pairs select`.

use std::path::PathBuf;

use clap::Args;

use crate::cli::ResourceOpts;
use crate::cli::common::{Context, Resources, append_pg, stdio_path};
use crate::error::{KiraError, Result};
use crate::io::buffered::lines;
use crate::metrics::{Metrics, Progress};
use crate::pairs::columns::ColumnMap;
use crate::pairs::reader::PairsReader;
use crate::pairs::record::{PairRecordRef, split_fields};
use crate::select::eval::parse_column_type;
use crate::select::{ColumnType, Filter};

/// Arguments for `select`.
#[derive(Debug, Args)]
pub struct SelectArgs {
    /// Filter expression, e.g. '(pair_type == "UU") and (abs(pos1-pos2) < 1000)'.
    pub condition: String,
    /// Input .pairs file; `-` or omitted = stdin.
    pub input: Option<PathBuf>,
    /// Output file for selected pairs; `-` or omitted = stdout.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// Output file for pairs that do not match.
    #[arg(long, value_name = "FILE")]
    pub output_rest: Option<PathBuf>,
    /// Restrict to chromosomes listed in this file (first column); also
    /// subsets #chromsize/#chromosomes header fields.
    #[arg(long, value_name = "FILE")]
    pub chrom_subset: Option<PathBuf>,
    /// Cast a column to a type (int, float, str) for the expression. May be repeated.
    #[arg(short = 't', long = "type-cast", num_args = 2, value_names = ["COLUMN", "TYPE"])]
    pub type_cast: Vec<String>,
    /// Comma-separated columns to drop from the output.
    #[arg(long, value_name = "COLUMNS")]
    pub remove_columns: Option<String>,
    /// Resource options.
    #[command(flatten)]
    pub res: ResourceOpts,
}

/// Run `select`.
pub fn run(a: SelectArgs, ctx: &Context) -> Result<()> {
    let res = Resources::from_opts(&a.res)?;
    let mut metrics = Metrics::start();
    let reader = PairsReader::open(stdio_path(&a.input), res.io_threads)?;
    let (mut header, cols, body) = reader.into_parts();
    append_pg(&mut header, "select", ctx)?;
    // Column removal.
    let mut scheme: Option<Vec<usize>> = None;
    if let Some(rc) = &a.remove_columns {
        let remove: Vec<&str> = rc
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .collect();
        for c in &remove {
            if crate::pairs::header::STANDARD_COLUMNS.contains(c) {
                log::warn!(
                    "removing required {c} column for .pairs format; output is not .pairs anymore"
                );
            }
        }
        let keep: Vec<(usize, String)> = cols
            .names()
            .iter()
            .enumerate()
            .filter(|(_, n)| !remove.contains(&n.as_str()))
            .map(|(i, n)| (i, n.clone()))
            .collect();
        if keep.len() == cols.len() {
            log::warn!("column(s) {rc} not in the file, --remove-columns has no effect");
        } else {
            header.set_columns(&keep.iter().map(|(_, n)| n.clone()).collect::<Vec<_>>());
            scheme = Some(keep.into_iter().map(|(i, _)| i).collect());
        }
    }
    let mut condition = a.condition.trim().to_string();
    if let Some(p) = &a.chrom_subset {
        let cs = crate::chroms::ChromSizes::from_path(p).or_else(|_| {
            // Accept a plain list of names too.
            let text = std::fs::read_to_string(p).map_err(|e| KiraError::io(p, e))?;
            let mut cs = crate::chroms::ChromSizes::default();
            for l in text.lines() {
                let n = l.split('\t').next().unwrap_or("").trim();
                if !n.is_empty() {
                    cs.push(n, 0);
                }
            }
            Ok::<_, KiraError>(cs)
        })?;
        let names: Vec<String> = cs.names().to_vec();
        header.subset_chromosomes(&names);
        let list = names
            .iter()
            .map(|n| format!("\"{}\"", n.replace('"', "\\\"")))
            .collect::<Vec<_>>()
            .join(", ");
        condition = format!("({condition}) and (chrom1 in [{list}]) and (chrom2 in [{list}])");
    }
    let casts: Vec<(String, ColumnType)> = a
        .type_cast
        .chunks(2)
        .map(|c| {
            if c.len() != 2 {
                return Err(KiraError::arg("--type-cast needs COLUMN TYPE"));
            }
            Ok((c[0].clone(), parse_column_type(&c[1])?))
        })
        .collect::<Result<_>>()?;
    let filter = Filter::compile(&condition, &cols, &casts)?;
    let mut writer = res.open_writer(stdio_path(&a.output))?;
    writer.write_header(&header)?;
    let mut rest = match &a.output_rest {
        Some(p) => {
            let mut w = res.open_writer(stdio_path(&Some(p.clone())))?;
            w.write_header(&header)?;
            Some(w)
        }
        None => None,
    };
    let mut progress = Progress::new(ctx.progress, std::time::Duration::from_secs(5));
    let mut body = body;
    let mut ends = Vec::with_capacity(32);
    let mut out = Vec::with_capacity(256);
    let mut records = 0u64;
    let mut selected = 0u64;
    let name = body.name().to_string();
    while let Some(block) = body.next_block()? {
        for (i, line) in lines(&block.data).enumerate() {
            if line.is_empty() {
                continue;
            }
            records += 1;
            split_fields(line, &mut ends);
            let rec = PairRecordRef::new(line, &ends);
            let passed = filter.eval(&rec).map_err(|e| {
                KiraError::Expression(format!(
                    "{e} ({} line {})",
                    name,
                    block.first_line + i as u64
                ))
            })?;
            let line_out: &[u8] = match &scheme {
                Some(s) => {
                    out.clear();
                    for (j, &c) in s.iter().enumerate() {
                        if j > 0 {
                            out.push(b'\t');
                        }
                        out.extend_from_slice(rec.field(c).unwrap_or(b""));
                    }
                    &out
                }
                None => line,
            };
            if passed {
                selected += 1;
                writer.write_line(line_out)?;
            } else if let Some(r) = rest.as_mut() {
                r.write_line(line_out)?;
            }
        }
        progress.tick(records, body.bytes_read());
    }
    let _ = ColumnMap::standard();
    let bytes_read = body.bytes_read();
    let bytes_written = writer.finish()?;
    if let Some(r) = rest {
        r.finish()?;
    }
    progress.finish(records, bytes_read);
    if ctx.metrics {
        metrics.set("records_read", records);
        metrics.set("records_written", selected);
        metrics.set("bytes_read", bytes_read);
        metrics.set("bytes_written", bytes_written);
        metrics.finalize();
        metrics.print();
    }
    Ok(())
}
