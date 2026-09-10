//! Parallel, memory-bounded external merge sort with pairtools block-sort
//! semantics (`chrom1`, `chrom2` lexicographic; `pos1`, `pos2` numeric;
//! `pair_type` lexicographic; input order for ties).

pub mod external;
pub mod key;
pub mod merge;
pub mod run;

pub use external::{ExternalSorter, SortConfig, SortMetrics, SortedStream};
pub use key::{ParsedChunk, SortEntry, SortKeyContext};
