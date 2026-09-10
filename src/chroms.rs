//! Chromosome sizes, chromosome ordering and the process-wide chromosome
//! name dictionary.

use std::collections::HashMap;
use std::path::Path;
use std::sync::RwLock;

use crate::error::{KiraError, Location, Result};

/// The chromosome name used by pairtools for an unmapped side.
pub const UNMAPPED_CHROM: &[u8] = b"!";
/// Position reported for an unmapped side.
pub const UNMAPPED_POS: u64 = 0;
/// Strand reported for an unmapped side.
pub const UNMAPPED_STRAND: u8 = b'-';

/// An ordered list of chromosomes with lengths (a `chrom.sizes` file).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ChromSizes {
    names: Vec<String>,
    sizes: Vec<u64>,
    index: HashMap<String, usize>,
}

impl ChromSizes {
    /// Read a UCSC style `chrom.sizes` file: `name<TAB>length` per line.
    ///
    /// Blank lines are skipped. A missing or non-numeric length is an error,
    /// as is a duplicated chromosome name.
    pub fn from_path(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path).map_err(|e| KiraError::io(path, e))?;
        Self::parse(&text, &path.display().to_string())
    }

    /// Parse from text; `name` is used in error messages.
    pub fn parse(text: &str, name: &str) -> Result<Self> {
        let mut cs = Self::default();
        for (i, raw) in text.lines().enumerate() {
            let line = raw.trim_end_matches('\r');
            if line.trim().is_empty() {
                continue;
            }
            let loc = || Location::file(name).at_line(i as u64 + 1);
            let mut it = line.split('\t');
            let chrom = it.next().unwrap_or("").trim();
            if chrom.is_empty() {
                return Err(KiraError::ChromSizes {
                    message: "empty chromosome name".into(),
                    location: loc(),
                });
            }
            let size_str = it
                .next()
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .ok_or_else(|| KiraError::ChromSizes {
                    message: format!("missing length for {chrom}"),
                    location: loc().at_column(2),
                })?;
            let size: u64 = size_str.parse().map_err(|_| KiraError::ChromSizes {
                message: format!("invalid length for {chrom}"),
                location: loc().at_column(2).with_value(size_str.as_bytes()),
            })?;
            if !cs.push(chrom, size) {
                return Err(KiraError::ChromSizes {
                    message: format!("duplicate chromosome {chrom}"),
                    location: loc(),
                });
            }
        }
        Ok(cs)
    }

    /// Append a chromosome. Returns false (and does nothing) on duplicates.
    pub fn push(&mut self, name: &str, size: u64) -> bool {
        if self.index.contains_key(name) {
            return false;
        }
        self.index.insert(name.to_string(), self.names.len());
        self.names.push(name.to_string());
        self.sizes.push(size);
        true
    }

    /// Number of chromosomes.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// True when there are no chromosomes.
    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    /// Chromosome names in file order.
    pub fn names(&self) -> &[String] {
        &self.names
    }

    /// Chromosome lengths in file order.
    pub fn sizes(&self) -> &[u64] {
        &self.sizes
    }

    /// Iterate `(name, size)` pairs in file order.
    pub fn iter(&self) -> impl Iterator<Item = (&str, u64)> + '_ {
        self.names
            .iter()
            .map(String::as_str)
            .zip(self.sizes.iter().copied())
    }

    /// Position of a chromosome in the file order.
    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.index.get(name).copied()
    }

    /// Length of a chromosome.
    pub fn size_of(&self, name: &str) -> Option<u64> {
        self.index_of(name).map(|i| self.sizes[i])
    }
}

/// Process-wide dictionary mapping chromosome names to dense ids.
///
/// Ids are assigned in first-seen order and never change. Callers that need
/// lexicographic ordering ask for a [`ChromRanks`] snapshot.
#[derive(Debug, Default)]
pub struct ChromDict {
    inner: RwLock<DictInner>,
}

#[derive(Debug, Default)]
struct DictInner {
    names: Vec<Box<[u8]>>,
    index: HashMap<Box<[u8]>, u32>,
}

impl ChromDict {
    /// Empty dictionary.
    pub fn new() -> Self {
        Self::default()
    }

    /// Dictionary pre-populated with the given names (ids follow the order).
    pub fn with_names<I, S>(names: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<[u8]>,
    {
        let d = Self::new();
        for n in names {
            d.intern(n.as_ref());
        }
        d
    }

    /// Look up the id of a name without inserting.
    pub fn get(&self, name: &[u8]) -> Option<u32> {
        let inner = self.inner.read().unwrap_or_else(|e| e.into_inner());
        inner.index.get(name).copied()
    }

    /// Return the id of a name, inserting it if unseen.
    pub fn intern(&self, name: &[u8]) -> u32 {
        if let Some(id) = self.get(name) {
            return id;
        }
        let mut inner = self.inner.write().unwrap_or_else(|e| e.into_inner());
        if let Some(id) = inner.index.get(name) {
            return *id;
        }
        let id = inner.names.len() as u32;
        let boxed: Box<[u8]> = name.into();
        inner.names.push(boxed.clone());
        inner.index.insert(boxed, id);
        id
    }

    /// Number of names in the dictionary.
    pub fn len(&self) -> usize {
        self.inner
            .read()
            .unwrap_or_else(|e| e.into_inner())
            .names
            .len()
    }

    /// True when the dictionary is empty.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Owned copy of a chromosome name.
    pub fn name(&self, id: u32) -> Vec<u8> {
        let inner = self.inner.read().unwrap_or_else(|e| e.into_inner());
        inner.names[id as usize].to_vec()
    }

    /// Snapshot of all names in id order.
    pub fn names(&self) -> Vec<Box<[u8]>> {
        self.inner
            .read()
            .unwrap_or_else(|e| e.into_inner())
            .names
            .clone()
    }

    /// Lexicographic (byte-wise, `LC_ALL=C`) ranks of all ids.
    pub fn ranks(&self) -> ChromRanks {
        let names = self.names();
        ChromRanks::from_names(&names)
    }
}

/// Lexicographic rank per chromosome id, used for pairtools-compatible sorting.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ChromRanks {
    ranks: Vec<u32>,
}

impl ChromRanks {
    /// Build ranks for names indexed by id.
    pub fn from_names(names: &[Box<[u8]>]) -> Self {
        let mut order: Vec<u32> = (0..names.len() as u32).collect();
        order.sort_by(|a, b| names[*a as usize].cmp(&names[*b as usize]));
        let mut ranks = vec![0u32; names.len()];
        for (rank, id) in order.into_iter().enumerate() {
            ranks[id as usize] = rank as u32;
        }
        Self { ranks }
    }

    /// Rank of an id. Ids beyond the snapshot sort last (they cannot occur in
    /// data that was produced before the snapshot).
    #[inline]
    pub fn rank(&self, id: u32) -> u32 {
        self.ranks.get(id as usize).copied().unwrap_or(u32::MAX)
    }

    /// Number of known ids.
    pub fn len(&self) -> usize {
        self.ranks.len()
    }

    /// True when empty.
    pub fn is_empty(&self) -> bool {
        self.ranks.is_empty()
    }

    /// Raw rank table.
    pub fn as_slice(&self) -> &[u32] {
        &self.ranks
    }
}

/// Chromosome enumeration used by `flip` and `parse` (pairtools
/// `headerops.get_chrom_order` + `chrom_enum`): `!` is 0, then the listed
/// chromosomes in order.
#[derive(Debug, Clone, Default)]
pub struct ChromOrder {
    order: HashMap<Vec<u8>, u32>,
    names: Vec<String>,
}

impl ChromOrder {
    /// Order from a chromosome sizes file (first column only matters).
    pub fn from_chromsizes(cs: &ChromSizes) -> Self {
        Self::from_names(cs.names().iter().map(String::as_str))
    }

    /// Order from names; `!` gets enumeration 0, names get 1..n.
    pub fn from_names<'a, I: IntoIterator<Item = &'a str>>(names: I) -> Self {
        let mut order = HashMap::new();
        order.insert(UNMAPPED_CHROM.to_vec(), 0u32);
        let mut list = Vec::new();
        let mut i = 1u32;
        for n in names {
            if n.is_empty() || order.contains_key(n.as_bytes()) {
                continue;
            }
            order.insert(n.as_bytes().to_vec(), i);
            list.push(n.to_string());
            i += 1;
        }
        Self { order, names: list }
    }

    /// pairtools `get_chrom_order(chroms_file, sam_chroms)`: listed
    /// chromosomes present in `present` keep file order, remaining present
    /// chromosomes follow in lexicographic order.
    pub fn from_names_restricted<'a, I, J>(listed: I, present: J) -> Self
    where
        I: IntoIterator<Item = &'a str>,
        J: IntoIterator<Item = &'a str>,
    {
        let present: Vec<&str> = present.into_iter().collect();
        let present_set: std::collections::HashSet<&str> = present.iter().copied().collect();
        let mut names: Vec<&str> = listed
            .into_iter()
            .filter(|n| !n.is_empty() && present_set.contains(n))
            .collect();
        let listed_set: std::collections::HashSet<&str> = names.iter().copied().collect();
        let mut remaining: Vec<&str> = present
            .iter()
            .copied()
            .filter(|n| !listed_set.contains(n))
            .collect();
        remaining.sort_unstable();
        remaining.dedup();
        names.extend(remaining);
        Self::from_names(names)
    }

    /// Enumeration value of a chromosome name, if annotated.
    #[inline]
    pub fn get(&self, name: &[u8]) -> Option<u32> {
        self.order.get(name).copied()
    }

    /// Annotated chromosome names, excluding `!`, in enumeration order.
    pub fn names(&self) -> &[String] {
        &self.names
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_chromsizes() {
        let cs = ChromSizes::parse("chr1\t100\nchr2\t50\n\nchrM\t7\n", "t").unwrap();
        assert_eq!(cs.len(), 3);
        assert_eq!(cs.size_of("chr2"), Some(50));
        assert_eq!(cs.index_of("chrM"), Some(2));
        assert!(ChromSizes::parse("chr1\tabc\n", "t").is_err());
        assert!(ChromSizes::parse("chr1\n", "t").is_err());
        assert!(ChromSizes::parse("chr1\t1\nchr1\t2\n", "t").is_err());
    }

    #[test]
    fn dict_and_ranks_are_lexicographic() {
        let d = ChromDict::new();
        let c2 = d.intern(b"chr2");
        let c10 = d.intern(b"chr10");
        let c1 = d.intern(b"chr1");
        assert_eq!(d.intern(b"chr2"), c2);
        let r = d.ranks();
        assert!(r.rank(c1) < r.rank(c10));
        assert!(r.rank(c10) < r.rank(c2));
        assert_eq!(d.name(c10), b"chr10");
    }

    #[test]
    fn chrom_order_matches_pairtools() {
        let o = ChromOrder::from_names(["chr1", "chr2"]);
        assert_eq!(o.get(b"!"), Some(0));
        assert_eq!(o.get(b"chr1"), Some(1));
        assert_eq!(o.get(b"chr2"), Some(2));
        assert_eq!(o.get(b"chrX"), None);
        let o = ChromOrder::from_names_restricted(
            ["chr2", "chrY", "chr1"],
            ["chr1", "chr2", "chrM", "chrA"],
        );
        assert_eq!(o.names(), &["chr2", "chr1", "chrA", "chrM"]);
    }
}
