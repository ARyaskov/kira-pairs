//! # kira-pairs
//!
//! A high-performance, memory-bounded implementation of the core
//! [pairtools](https://github.com/open2c/pairtools) Hi-C `.pairs` workflow:
//! sorting, duplicate removal, statistics, flipping, filtering, binning and
//! standard paired-end BAM parsing, plus a fused `process` pipeline.
//!
//! The crate is a library first; the `kira-pairs` binary is a thin CLI
//! wrapper. Public building blocks:
//!
//! * [`pairs::PairsReader`] / [`pairs::PairsWriter`] / [`pairs::Header`]
//! * [`sort::ExternalSorter`]
//! * [`dedup::Deduper`]
//! * [`stats::StatsAccumulator`]
//! * [`binning::Binner`]
//! * [`select::Filter`]
//! * [`parse::HicParser`]
//!
//! The private run format used for temporary files (`sort::run`) is
//! **not** a stable public format and may change between versions.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod binning;
pub mod chroms;
pub mod cli;
pub mod dedup;
pub mod error;
pub mod flip;
pub mod generate;
pub mod io;
pub mod logging;
pub mod memory;
pub mod metrics;
pub mod pairs;
pub mod parse;
pub mod pipeline;
pub mod select;
pub mod sort;
pub mod stats;
pub mod util;

pub use error::{KiraError, Result};

/// Crate version string.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
/// Program name used in `@PG` records.
pub const PROGRAM_NAME: &str = "kira-pairs";
