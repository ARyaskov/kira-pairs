//! Streaming pairtools-compatible statistics.
//!
//! [`StatsAccumulator::observe`] is a cheap per-record update (no
//! allocation, no locking) so that it can run inside other pipelines;
//! accumulators from worker threads are combined with
//! [`StatsAccumulator::merge`].

pub mod accumulator;
pub mod format;
pub mod lambertw;

pub use accumulator::{DistBins, StatsAccumulator, StatsSnapshot};
pub use format::{StatsFormat, Value};
