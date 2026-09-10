//! `kira-pairs generate`.

use std::path::PathBuf;

use clap::Args;

use crate::cli::ResourceOpts;
use crate::cli::common::{Context, Resources, stdio_path};
use crate::error::Result;
use crate::generate::{GenerateConfig, generate};

/// Arguments for `generate`.
#[derive(Debug, Args)]
pub struct GenerateArgs {
    /// Output file; `-` or omitted = stdout.
    #[arg(short, long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// Number of records.
    #[arg(long, default_value_t = 1_000_000)]
    pub records: u64,
    /// RNG seed.
    #[arg(long, default_value_t = 42)]
    pub seed: u64,
    /// Number of chromosomes.
    #[arg(long, default_value_t = 24)]
    pub chromosomes: usize,
    /// Chromosome length (bp).
    #[arg(long, default_value_t = 100_000_000)]
    pub chrom_length: u64,
    /// Fraction of cis pairs.
    #[arg(long, default_value_t = 0.75)]
    pub cis_fraction: f64,
    /// Fraction of near-duplicate records.
    #[arg(long, default_value_t = 0.1)]
    pub duplicate_rate: f64,
    /// Maximum per-side offset of duplicates (bp).
    #[arg(long, default_value_t = 2)]
    pub duplicate_radius: u64,
    /// Fraction of records with an unmapped side.
    #[arg(long, default_value_t = 0.02)]
    pub unmapped_fraction: f64,
    /// Number of extra integer columns.
    #[arg(long, default_value_t = 0)]
    pub extra_columns: usize,
    /// Read ID length.
    #[arg(long, default_value_t = 24)]
    pub readid_length: usize,
    /// Emit block-sorted, upper-triangular records (kept in memory).
    #[arg(long)]
    pub sorted: bool,
    /// Resource options (compression threads apply to the output).
    #[command(flatten)]
    pub res: ResourceOpts,
}

/// Run `generate`.
pub fn run(a: GenerateArgs, _ctx: &Context) -> Result<()> {
    let res = Resources::from_opts(&a.res)?;
    let cfg = GenerateConfig {
        records: a.records,
        seed: a.seed,
        chromosomes: a.chromosomes,
        chrom_length: a.chrom_length,
        cis_fraction: a.cis_fraction,
        duplicate_rate: a.duplicate_rate,
        duplicate_radius: a.duplicate_radius,
        unmapped_fraction: a.unmapped_fraction,
        extra_columns: a.extra_columns,
        readid_length: a.readid_length,
        sorted: a.sorted,
    };
    let mut out = crate::io::compression::open_output(stdio_path(&a.output), res.output_options())?;
    let n = generate(&cfg, &mut out)?;
    out.finish()?;
    log::info!("generated {n} records");
    Ok(())
}
