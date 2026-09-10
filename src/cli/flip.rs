//! `kira-pairs flip`.

use std::path::PathBuf;
use std::sync::Arc;

use clap::Args;

use crate::chroms::{ChromOrder, ChromSizes};
use crate::cli::ResourceOpts;
use crate::cli::common::{Context, Resources, append_pg, new_dict, stdio_path};
use crate::error::Result;
use crate::flip::Flipper;
use crate::metrics::{Metrics, Progress};
use crate::pairs::parallel::OrderedParser;
use crate::pairs::reader::PairsReader;
use crate::pairs::record::{PairRecordRef, split_fields};

/// Arguments for `flip`.
#[derive(Debug, Args)]
pub struct FlipArgs {
    /// Input .pairs file; `-` or omitted = stdin.
    pub input: Option<PathBuf>,
    /// Output file; `-` or omitted = stdout.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// Chromosome order: a chrom.sizes-like file whose first column lists
    /// chromosome names.
    #[arg(short = 'c', long, value_name = "FILE")]
    pub chroms_path: PathBuf,
    /// Resource options.
    #[command(flatten)]
    pub res: ResourceOpts,
}

/// Run `flip`.
pub fn run(a: FlipArgs, ctx: &Context) -> Result<()> {
    let res = Resources::from_opts(&a.res)?;
    let mut metrics = Metrics::start();
    let cs = ChromSizes::from_path(&a.chroms_path)?;
    let order = ChromOrder::from_chromsizes(&cs);
    let reader = PairsReader::open(stdio_path(&a.input), res.io_threads)?;
    let (mut header, cols, body) = reader.into_parts();
    append_pg(&mut header, "flip", ctx)?;
    let dict = new_dict();
    let mut flipper = Flipper::new(order, &cols, Arc::clone(&dict));
    let mut writer = res.open_writer(stdio_path(&a.output))?;
    writer.write_header(&header)?;
    let mut parser = OrderedParser::with_depth(
        body,
        Arc::new(cols),
        Arc::clone(&dict),
        res.threads,
        res.budget.channel_depth(res.threads),
    );
    let mut progress = Progress::new(ctx.progress, std::time::Duration::from_secs(5));
    let mut ends = Vec::with_capacity(32);
    let mut out = Vec::with_capacity(256);
    let mut records = 0u64;
    let mut flipped = 0u64;
    while let Some(chunk) = parser.next_chunk()? {
        for e in &chunk.entries {
            let line = chunk.line(e);
            if flipper.needs_flip(&e.key) {
                split_fields(line, &mut ends);
                flipper.flip_line(&PairRecordRef::new(line, &ends), &mut out);
                writer.write_line(&out)?;
                flipped += 1;
            } else {
                writer.write_line(line)?;
            }
        }
        records += chunk.entries.len() as u64;
        progress.tick(records, parser.bytes_read());
    }
    let bytes_read = parser.bytes_read();
    drop(parser);
    let bytes_written = writer.finish()?;
    progress.finish(records, bytes_read);
    if ctx.metrics {
        metrics.set("records_read", records);
        metrics.set("records_flipped", flipped);
        metrics.set("bytes_read", bytes_read);
        metrics.set("bytes_written", bytes_written);
        metrics.finalize();
        metrics.print();
    }
    Ok(())
}
