//! A safe filter-expression engine for `kira-pairs select` (a subset of the
//! Python expressions accepted by `pairtools select`).
//!
//! Expressions are parsed once into an AST with column references resolved
//! to indices, then evaluated per record without allocation for the common
//! comparison/arithmetic cases.

pub mod eval;
pub mod parser;

pub use eval::{ColumnType, Filter, Value};
