//! Column name resolution for `.pairs` bodies.

use crate::error::{KiraError, Location, Result};
use crate::pairs::header::STANDARD_COLUMNS;

/// Resolved indices of the standard columns plus a name → index table.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ColumnMap {
    names: Vec<String>,
    /// Index of `readID`.
    pub readid: usize,
    /// Index of `chrom1`.
    pub chrom1: usize,
    /// Index of `pos1`.
    pub pos1: usize,
    /// Index of `chrom2`.
    pub chrom2: usize,
    /// Index of `pos2`.
    pub pos2: usize,
    /// Index of `strand1`.
    pub strand1: usize,
    /// Index of `strand2`.
    pub strand2: usize,
    /// `pair_type` may be absent in minimal 7-column files.
    pub pair_type: Option<usize>,
    /// Highest index among the standard columns; rows must have at least
    /// `min_fields` fields.
    pub min_fields: usize,
    /// Whether the names came from a `#columns:` header line.
    pub from_header: bool,
}

impl ColumnMap {
    /// Column map from header column names; falls back to the standard eight
    /// when the header has no `#columns:` line.
    pub fn from_names(names: Vec<String>) -> Result<Self> {
        let (names, from_header) = if names.is_empty() {
            (
                STANDARD_COLUMNS.iter().map(|s| s.to_string()).collect(),
                false,
            )
        } else {
            (names, true)
        };
        let find = |n: &str| names.iter().position(|c| c == n);
        let req = |n: &str| {
            find(n).ok_or_else(|| KiraError::MissingColumn {
                column: n.to_string(),
                location: Location::default(),
            })
        };
        let map = Self {
            readid: req("readID")?,
            chrom1: req("chrom1")?,
            pos1: req("pos1")?,
            chrom2: req("chrom2")?,
            pos2: req("pos2")?,
            strand1: req("strand1")?,
            strand2: req("strand2")?,
            pair_type: find("pair_type"),
            min_fields: 0,
            from_header,
            names,
        };
        let mut min = [
            map.readid,
            map.chrom1,
            map.pos1,
            map.chrom2,
            map.pos2,
            map.strand1,
            map.strand2,
        ]
        .into_iter()
        .max()
        .unwrap_or(0)
            + 1;
        if let Some(pt) = map.pair_type
            && from_header
        {
            min = min.max(pt + 1);
        }
        Ok(Self {
            min_fields: min,
            ..map
        })
    }

    /// Standard eight-column map.
    pub fn standard() -> Self {
        // Standard names always resolve.
        Self::from_names(Vec::new()).unwrap_or_else(|_| unreachable!())
    }

    /// All column names.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Number of named columns.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// True when no columns are named.
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// Index of a column by name or 0-based numeric string.
    pub fn index_of(&self, name: &str) -> Option<usize> {
        if let Ok(i) = name.parse::<usize>() {
            return Some(i);
        }
        self.names.iter().position(|c| c == name)
    }

    /// Index of a column or an error naming it.
    pub fn require(&self, name: &str) -> Result<usize> {
        self.index_of(name).ok_or_else(|| KiraError::MissingColumn {
            column: name.to_string(),
            location: Location::default(),
        })
    }

    /// Pairs of side-specific columns `(xxx1, xxx2)` present in the map, in
    /// column order (pairtools `flip` swaps exactly these).
    pub fn side_pairs(&self) -> Vec<(usize, usize)> {
        self.names
            .iter()
            .enumerate()
            .filter_map(|(i, n)| {
                let stem = n.strip_suffix('1')?;
                let partner = format!("{stem}2");
                self.names
                    .iter()
                    .position(|c| *c == partner)
                    .map(|j| (i, j))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_map() {
        let m = ColumnMap::standard();
        assert_eq!(m.chrom1, 1);
        assert_eq!(m.pair_type, Some(7));
        assert_eq!(m.min_fields, 7);
        assert!(!m.from_header);
    }

    #[test]
    fn custom_order_and_side_pairs() {
        let names: Vec<String> = [
            "readID",
            "chrom1",
            "pos1",
            "chrom2",
            "pos2",
            "strand1",
            "strand2",
            "pair_type",
            "mapq1",
            "mapq2",
            "phase1",
            "other",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let m = ColumnMap::from_names(names).unwrap();
        assert_eq!(m.min_fields, 8);
        assert_eq!(m.side_pairs(), vec![(1, 3), (2, 4), (5, 6), (8, 9)]);
        assert_eq!(m.index_of("9"), Some(9));
        assert_eq!(m.index_of("mapq2"), Some(9));
        assert!(m.require("nope").is_err());
    }

    #[test]
    fn missing_required() {
        let names: Vec<String> = ["readID", "chrom1"].iter().map(|s| s.to_string()).collect();
        assert!(ColumnMap::from_names(names).is_err());
    }
}
