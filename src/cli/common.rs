//! Helpers shared by the command implementations.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::cli::ResourceOpts;
use crate::error::{KiraError, Result};
use crate::io::compression::OutputOptions;
use crate::memory::{MemoryBudget, parse_size};
use crate::pairs::header::Header;
use crate::pairs::writer::PairsWriter;

/// Per-invocation context.
#[derive(Debug, Clone)]
pub struct Context {
    /// Print metrics at exit.
    pub metrics: bool,
    /// Print progress periodically.
    pub progress: bool,
    /// Command line for `@PG CL:`.
    pub argv: String,
}

/// Resolved resource settings.
#[derive(Debug, Clone)]
pub struct Resources {
    /// Worker threads.
    pub threads: usize,
    /// Input decompression threads.
    pub io_threads: usize,
    /// Output compression threads.
    pub compression_threads: usize,
    /// Compression level.
    pub compression_level: u32,
    /// Memory budget.
    pub budget: MemoryBudget,
    /// Temporary directory.
    pub tmpdir: Option<PathBuf>,
}

/// Default thread count: available parallelism.
pub fn default_threads() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

impl Resources {
    /// Resolve options.
    pub fn from_opts(o: &ResourceOpts) -> Result<Self> {
        let threads = o.threads.max(1);
        let io_threads = o.io_threads.unwrap_or(threads.min(4)).max(1);
        let compression_threads = o.compression_threads.unwrap_or(threads).max(1);
        if !(1..=9).contains(&o.compression_level) {
            return Err(KiraError::arg(
                "--compression-level must be between 1 and 9",
            ));
        }
        let total = parse_size(&o.memory)?;
        let budget = MemoryBudget::new(total, threads, io_threads)?;
        Ok(Self {
            threads,
            io_threads,
            compression_threads,
            compression_level: o.compression_level,
            budget,
            tmpdir: o.tmpdir.clone(),
        })
    }

    /// Output options for pairs writers.
    pub fn output_options(&self) -> OutputOptions {
        OutputOptions {
            compression: None,
            level: self.compression_level,
            threads: self.compression_threads,
        }
    }

    /// Open a pairs writer.
    pub fn open_writer(&self, path: Option<&Path>) -> Result<PairsWriter> {
        PairsWriter::create(path, self.output_options())
    }
}

/// Normalise an optional path argument: `None` and `-` mean stdio.
pub fn stdio_path(p: &Option<PathBuf>) -> Option<&Path> {
    match p {
        None => None,
        Some(x) if x.as_os_str() == "-" || x.as_os_str().is_empty() => None,
        Some(x) => Some(x.as_path()),
    }
}

/// Append the kira-pairs `@PG` record for a command.
pub fn append_pg(header: &mut Header, command: &str, ctx: &Context) -> Result<()> {
    if header.is_valid_pairs_header() {
        header.append_pg(
            &format!("kira-pairs_{command}"),
            "kira-pairs",
            &ctx.argv,
            crate::VERSION,
        )?;
    } else if !header.is_empty() {
        log::warn!("input header lacks '## pairs format' line; provenance (@PG) not recorded");
    }
    Ok(())
}

/// Parse a comma-separated list of pair types.
pub fn parse_pair_types(s: &str) -> std::collections::HashSet<Vec<u8>> {
    s.split(',')
        .map(str::trim)
        .filter(|x| !x.is_empty())
        .map(|x| x.as_bytes().to_vec())
        .collect()
}

/// Shared chromosome dictionary factory.
pub fn new_dict() -> Arc<crate::chroms::ChromDict> {
    let d = crate::chroms::ChromDict::new();
    d.intern(crate::chroms::UNMAPPED_CHROM);
    Arc::new(d)
}
