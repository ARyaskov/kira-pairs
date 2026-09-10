//! The 4DN/pairtools `.pairs` format: header handling, column maps, record
//! views and streaming readers/writers.

pub mod columns;
pub mod header;
pub mod parallel;
pub mod reader;
pub mod record;
pub mod writer;

pub use columns::ColumnMap;
pub use header::Header;
pub use parallel::OrderedParser;
pub use reader::{BodyBlocks, PairsReader};
pub use record::{PairKey, PairRecordRef, parse_line, split_fields};
pub use writer::PairsWriter;
