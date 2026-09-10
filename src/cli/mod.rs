//! Command-line interface (clap derive) and command dispatch.

pub mod bin;
pub mod common;
pub mod dedup;
pub mod flip;
pub mod generate;
pub mod parse;
pub mod process;
pub mod select;
pub mod sort;
pub mod stats;

use std::path::PathBuf;

use clap::{Args, Parser, Subcommand, ValueEnum};

use crate::error::Result;

/// kira-pairs: fast, memory-bounded Hi-C `.pairs` processing compatible
/// with pairtools.
#[derive(Debug, Parser)]
#[command(name = "kira-pairs", version, about, long_about = None, propagate_version = true)]
pub struct Cli {
    /// Increase diagnostic verbosity (-v, -vv).
    #[arg(short, long, action = clap::ArgAction::Count, global = true)]
    pub verbose: u8,
    /// Only print errors.
    #[arg(short, long, global = true)]
    pub quiet: bool,
    /// Print runtime metrics to stderr when done.
    #[arg(long, global = true)]
    pub metrics: bool,
    /// Print periodic progress to stderr.
    #[arg(long, global = true)]
    pub progress: bool,
    /// Subcommand.
    #[command(subcommand)]
    pub command: Command,
}

/// Resource options shared by all commands.
#[derive(Debug, Clone, Args)]
pub struct ResourceOpts {
    /// Worker threads for parsing, sorting and compression.
    #[arg(long, default_value_t = common::default_threads(), value_name = "N")]
    pub threads: usize,
    /// Memory budget for buffered data (e.g. 512M, 8G).
    #[arg(long, default_value = "2G", value_name = "SIZE")]
    pub memory: String,
    /// Directory for temporary files.
    #[arg(long, value_name = "PATH")]
    pub tmpdir: Option<PathBuf>,
    /// Threads used for reading/decompressing input [default: min(threads, 4)].
    #[arg(long, value_name = "N")]
    pub io_threads: Option<usize>,
    /// Threads used for compressing output [default: threads].
    #[arg(long, value_name = "N")]
    pub compression_threads: Option<usize>,
    /// gzip/BGZF compression level (1-9).
    #[arg(long, default_value_t = 6, value_name = "LEVEL")]
    pub compression_level: u32,
}

/// Duplicate distance metric.
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum MethodArg {
    /// max(|dpos1|, |dpos2|) <= max-mismatch
    Max,
    /// |dpos1| + |dpos2| <= max-mismatch
    Sum,
}

/// Duplicate clustering mode.
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum ClusteringArg {
    /// Connected components (pairtools scipy/sklearn backends)
    Transitive,
    /// Compare only with earlier non-duplicates (pairtools cython backend)
    Greedy,
}

/// pairtools backend names accepted for compatibility.
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum BackendArg {
    /// Transitive clustering.
    Scipy,
    /// Transitive clustering.
    Sklearn,
    /// Greedy clustering.
    Cython,
}

/// Walks policy for `parse`.
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum WalksPolicyArg {
    /// Mask unrescuable walks (chrom `!`, pos 0, strand `-`, type `WW`).
    Mask,
    /// Report the 5'-most alignment on each side.
    #[value(name = "5any")]
    FiveAny,
    /// Report the 5'-most unique alignment on each side.
    #[value(name = "5unique")]
    FiveUnique,
    /// Report the 3'-most alignment on each side.
    #[value(name = "3any")]
    ThreeAny,
    /// Report the 3'-most unique alignment on each side.
    #[value(name = "3unique")]
    ThreeUnique,
}

/// Bin output format.
#[derive(Debug, Clone, Copy, ValueEnum, PartialEq, Eq)]
pub enum BinFormatArg {
    /// `bin1_id<TAB>bin2_id<TAB>count` (cooler `load -f coo`)
    Coo,
    /// `chrom1 start1 end1 chrom2 start2 end2 count` (cooler `load -f bg2`)
    Bg2,
}

/// Subcommands.
#[derive(Debug, Subcommand)]
pub enum Command {
    /// Sort pairs (chrom1, chrom2 lexicographic; pos1, pos2 numeric; pair_type).
    Sort(sort::SortArgs),
    /// Flip pairs into upper-triangular order.
    Flip(flip::FlipArgs),
    /// Find and remove PCR/optical duplicates in sorted pairs.
    Dedup(dedup::DedupArgs),
    /// Calculate pairtools-compatible statistics (or merge stats files).
    Stats(stats::StatsArgs),
    /// Select pairs matching a filter expression.
    Select(select::SelectArgs),
    /// Aggregate pairs into fixed-resolution contact bins.
    Bin(bin::BinArgs),
    /// Parse paired-end Hi-C alignments (SAM/BAM) into pairs.
    Parse(parse::ParseArgs),
    /// Fused parse -> flip -> sort -> dedup -> stats/bin pipeline.
    Process(process::ProcessArgs),
    /// Generate a deterministic synthetic .pairs dataset.
    Generate(generate::GenerateArgs),
}

/// Run the CLI.
pub fn run(cli: Cli) -> Result<()> {
    let ctx = common::Context {
        metrics: cli.metrics,
        progress: cli.progress,
        argv: std::env::args().collect::<Vec<_>>().join(" "),
    };
    match cli.command {
        Command::Sort(a) => sort::run(a, &ctx),
        Command::Flip(a) => flip::run(a, &ctx),
        Command::Dedup(a) => dedup::run(a, &ctx),
        Command::Stats(a) => stats::run(a, &ctx),
        Command::Select(a) => select::run(a, &ctx),
        Command::Bin(a) => bin::run(a, &ctx),
        Command::Parse(a) => parse::run(a, &ctx),
        Command::Process(a) => process::run(a, &ctx),
        Command::Generate(a) => generate::run(a, &ctx),
    }
}
