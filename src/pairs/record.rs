//! Record representations: a zero-copy view over a line and a compact key.

use memchr::memchr_iter;

use crate::chroms::{ChromDict, UNMAPPED_CHROM};
use crate::error::{KiraError, Location, Result};
use crate::pairs::columns::ColumnMap;
use crate::util::int::parse_u64;

/// Number of pair-type bytes stored inline in a [`PairKey`].
pub const PAIR_TYPE_INLINE: usize = 8;

/// Compact, fixed-size key of a pair record.
///
/// Holds every field the hot paths (sort, dedup, stats, bin, flip) need so
/// that the text line is never re-parsed by those stages.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(C)]
pub struct PairKey {
    /// Monotonic input sequence number (0-based), used for stable ties.
    pub seq: u64,
    /// Position on side 1 (1-based, 0 when unmapped).
    pub pos1: u64,
    /// Position on side 2.
    pub pos2: u64,
    /// Dictionary id of chrom1.
    pub chrom1: u32,
    /// Dictionary id of chrom2.
    pub chrom2: u32,
    /// Pair type, zero-padded. Longer pair types keep their first 8 bytes
    /// and set [`PairKey::PT_TRUNCATED`].
    pub pair_type: [u8; PAIR_TYPE_INLINE],
    /// Strand byte of side 1 (`+`, `-`).
    pub strand1: u8,
    /// Strand byte of side 2.
    pub strand2: u8,
    /// Flag bits.
    pub flags: u8,
    /// Length of the pair type (0 when the column is absent).
    pub pair_type_len: u8,
}

impl PairKey {
    /// Set when the pair type did not fit inline.
    pub const PT_TRUNCATED: u8 = 1;

    /// Inline pair type bytes (possibly truncated).
    #[inline]
    pub fn pair_type_bytes(&self) -> &[u8] {
        let n = (self.pair_type_len as usize).min(PAIR_TYPE_INLINE);
        &self.pair_type[..n]
    }

    /// True when the pair type had to be truncated.
    #[inline]
    pub fn pair_type_truncated(&self) -> bool {
        self.flags & Self::PT_TRUNCATED != 0
    }

    /// True when the pair type is exactly `DD`.
    #[inline]
    pub fn is_dd(&self) -> bool {
        self.pair_type_len == 2 && self.pair_type[0] == b'D' && self.pair_type[1] == b'D'
    }

    /// Strand pair as a 2-byte array, e.g. `b"+-"`.
    #[inline]
    pub fn strands(&self) -> [u8; 2] {
        [self.strand1, self.strand2]
    }
}

/// Fill `ends` with the exclusive end offset of every TAB-separated field.
///
/// The number of fields is `ends.len()`; field `i` spans
/// `start(i)..ends[i]` where `start(0) = 0` and `start(i) = ends[i-1] + 1`.
#[inline]
pub fn split_fields(line: &[u8], ends: &mut Vec<u32>) {
    ends.clear();
    for p in memchr_iter(b'\t', line) {
        ends.push(p as u32);
    }
    ends.push(line.len() as u32);
}

/// Zero-copy view over one body line and its field boundaries.
#[derive(Debug, Clone, Copy)]
pub struct PairRecordRef<'a> {
    line: &'a [u8],
    ends: &'a [u32],
}

impl<'a> PairRecordRef<'a> {
    /// Build a view from a line (without newline) and its field ends.
    #[inline]
    pub fn new(line: &'a [u8], ends: &'a [u32]) -> Self {
        Self { line, ends }
    }

    /// Full line bytes without trailing newline.
    #[inline]
    pub fn line(&self) -> &'a [u8] {
        self.line
    }

    /// Number of fields.
    #[inline]
    pub fn n_fields(&self) -> usize {
        self.ends.len()
    }

    /// Byte range of field `i`, if present.
    #[inline]
    pub fn field_range(&self, i: usize) -> Option<(usize, usize)> {
        let end = *self.ends.get(i)? as usize;
        let start = if i == 0 {
            0
        } else {
            self.ends[i - 1] as usize + 1
        };
        Some((start, end))
    }

    /// Bytes of field `i`, if present.
    #[inline]
    pub fn field(&self, i: usize) -> Option<&'a [u8]> {
        self.field_range(i).map(|(s, e)| &self.line[s..e])
    }

    /// Field ends slice.
    #[inline]
    pub fn ends(&self) -> &'a [u32] {
        self.ends
    }
}

/// Parse a body line into a [`PairKey`] using the column map and dictionary.
///
/// `ends` is a scratch buffer that receives field boundaries and can be
/// reused across calls. `loc` builds the error location lazily.
pub fn parse_line(
    line: &[u8],
    cols: &ColumnMap,
    dict: &ChromDict,
    seq: u64,
    ends: &mut Vec<u32>,
    loc: impl Fn() -> Location,
) -> Result<PairKey> {
    split_fields(line, ends);
    let rec = PairRecordRef::new(line, ends);
    parse_record(&rec, cols, dict, seq, loc)
}

/// Parse a [`PairKey`] from an already-split record.
pub fn parse_record(
    rec: &PairRecordRef<'_>,
    cols: &ColumnMap,
    dict: &ChromDict,
    seq: u64,
    loc: impl Fn() -> Location,
) -> Result<PairKey> {
    if rec.n_fields() < cols.min_fields {
        return Err(KiraError::format(
            format!(
                "expected at least {} fields, found {}",
                cols.min_fields,
                rec.n_fields()
            ),
            loc(),
        ));
    }
    let field = |i: usize| rec.field(i).unwrap_or(b"");
    let pos = |i: usize| -> Result<u64> {
        let f = field(i);
        parse_u64(f).ok_or_else(|| {
            KiraError::format("invalid position", loc().at_column(i + 1).with_value(f))
        })
    };
    let strand = |i: usize| -> Result<u8> {
        let f = field(i);
        match f {
            [b'+'] | [b'-'] => Ok(f[0]),
            // pairtools itself only ever writes '+' or '-', but it never
            // validates strands on input; '.' is common in third-party files.
            [b'.'] => Ok(b'.'),
            _ => Err(KiraError::format(
                "invalid strand (expected '+' or '-')",
                loc().at_column(i + 1).with_value(f),
            )),
        }
    };
    let chrom = |i: usize| -> Result<u32> {
        let f = field(i);
        if f.is_empty() {
            return Err(KiraError::format(
                "empty chromosome name",
                loc().at_column(i + 1),
            ));
        }
        Ok(dict.intern(f))
    };
    let chrom1 = chrom(cols.chrom1)?;
    let chrom2 = chrom(cols.chrom2)?;
    let pos1 = pos(cols.pos1)?;
    let pos2 = pos(cols.pos2)?;
    let strand1 = strand(cols.strand1)?;
    let strand2 = strand(cols.strand2)?;
    let mut key = PairKey {
        seq,
        pos1,
        pos2,
        chrom1,
        chrom2,
        pair_type: [0; PAIR_TYPE_INLINE],
        strand1,
        strand2,
        flags: 0,
        pair_type_len: 0,
    };
    if let Some(pt) = cols.pair_type
        && let Some(f) = rec.field(pt)
    {
        let n = f.len().min(PAIR_TYPE_INLINE);
        key.pair_type[..n].copy_from_slice(&f[..n]);
        key.pair_type_len = f.len().min(u8::MAX as usize) as u8;
        if f.len() > PAIR_TYPE_INLINE {
            key.flags |= PairKey::PT_TRUNCATED;
        }
    }
    Ok(key)
}

/// True when the chromosome name denotes an unmapped side.
#[inline]
pub fn is_unmapped_name(name: &[u8]) -> bool {
    name == UNMAPPED_CHROM
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splits_fields() {
        let mut ends = Vec::new();
        split_fields(b"a\tbb\t\tc", &mut ends);
        let r = PairRecordRef::new(b"a\tbb\t\tc", &ends);
        assert_eq!(r.n_fields(), 4);
        assert_eq!(r.field(0), Some(&b"a"[..]));
        assert_eq!(r.field(1), Some(&b"bb"[..]));
        assert_eq!(r.field(2), Some(&b""[..]));
        assert_eq!(r.field(3), Some(&b"c"[..]));
        assert_eq!(r.field(4), None);
    }

    #[test]
    fn parses_standard_line() {
        let dict = ChromDict::new();
        let cols = ColumnMap::standard();
        let mut ends = Vec::new();
        let k = parse_line(
            b"r1\tchr1\t100\tchr2\t200\t+\t-\tUU\tx",
            &cols,
            &dict,
            7,
            &mut ends,
            Location::default,
        )
        .unwrap();
        assert_eq!(k.pos1, 100);
        assert_eq!(k.pos2, 200);
        assert_eq!(k.strand1, b'+');
        assert_eq!(k.strand2, b'-');
        assert_eq!(k.pair_type_bytes(), b"UU");
        assert_eq!(k.seq, 7);
        assert_eq!(dict.name(k.chrom1), b"chr1");
        assert_eq!(dict.name(k.chrom2), b"chr2");
    }

    #[test]
    fn rejects_bad_lines() {
        let dict = ChromDict::new();
        let cols = ColumnMap::standard();
        let mut ends = Vec::new();
        let e = parse_line(
            b"r1\tchr1\t1x\tchr2\t200\t+\t-\tUU",
            &cols,
            &dict,
            0,
            &mut ends,
            || Location::file("f").at_line(3),
        )
        .unwrap_err();
        let msg = e.to_string();
        assert!(
            msg.contains("line 3") && msg.contains("column 3") && msg.contains("1x"),
            "{msg}"
        );
        assert!(
            parse_line(
                b"r1\tchr1\t1\tchr2\t200\t*\t-\tUU",
                &cols,
                &dict,
                0,
                &mut ends,
                Location::default
            )
            .is_err()
        );
        assert!(
            parse_line(
                b"r1\tchr1\t1",
                &cols,
                &dict,
                0,
                &mut ends,
                Location::default
            )
            .is_err()
        );
        // seven columns are acceptable (pair_type optional without header).
        let k = parse_line(
            b"r1\tchr1\t1\tchr2\t2\t+\t+",
            &cols,
            &dict,
            0,
            &mut ends,
            Location::default,
        )
        .unwrap();
        assert_eq!(k.pair_type_len, 0);
    }

    #[test]
    fn long_pair_type_is_truncated() {
        let dict = ChromDict::new();
        let cols = ColumnMap::standard();
        let mut ends = Vec::new();
        let k = parse_line(
            b"r\t!\t0\t!\t0\t-\t-\tabcdefghij",
            &cols,
            &dict,
            0,
            &mut ends,
            Location::default,
        )
        .unwrap();
        assert!(k.pair_type_truncated());
        assert_eq!(k.pair_type_bytes(), b"abcdefgh");
        assert_eq!(k.pair_type_len, 10);
    }
}
