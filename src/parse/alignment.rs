//! Conversion of SAM/BAM records into pairtools-style alignments.

use noodles_sam as sam;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::{Tag, Value};

use crate::error::{KiraError, Location, Result};

/// A parsed alignment (pairtools' alignment dict).
#[derive(Debug, Clone, PartialEq)]
pub struct Alignment {
    /// Reference sequence index (`None` = unmapped `!`).
    pub chrom: Option<usize>,
    /// 5' end (1-based).
    pub pos5: u64,
    /// 3' end (1-based).
    pub pos3: u64,
    /// Reported position.
    pub pos: u64,
    /// `+`/`-`.
    pub strand: u8,
    /// MAPQ.
    pub mapq: u8,
    /// Mapped.
    pub is_mapped: bool,
    /// MAPQ >= min_mapq.
    pub is_unique: bool,
    /// No `SA` tag.
    pub is_linear: bool,
    /// Clipped bases at the read's 5' end.
    pub dist_to_5: u64,
    /// Clipped bases at the read's 3' end.
    pub dist_to_3: u64,
    /// CIGAR text (`*` for empty alignments, `None` when a record had no
    /// CIGAR, mirroring pysam's `cigarstring`).
    pub cigar_text: String,
    /// Reference span.
    pub algn_ref_span: u64,
    /// Read span.
    pub algn_read_span: u64,
    /// Matched bases.
    pub matched_bp: u64,
    /// Leading clip.
    pub clip5_ref: u64,
    /// Trailing clip.
    pub clip3_ref: u64,
    /// Read length.
    pub read_len: u64,
    /// Type letter (`N`, `M`, `U`, `X`, `R`, `W`).
    pub kind: u8,
    /// Requested SAM tag values (one per requested tag).
    pub tags: Vec<String>,
    /// Read sequence when requested.
    pub seq: Option<Vec<u8>>,
    /// Side (1/2) after normalisation.
    pub read_side: Option<u64>,
    /// Index on its side after normalisation.
    pub algn_idx: Option<u64>,
    /// Number of alignments on its side.
    pub same_side_count: Option<u64>,
}

impl Alignment {
    /// pairtools `empty_alignment()`.
    pub fn empty() -> Self {
        Self {
            chrom: None,
            pos5: 0,
            pos3: 0,
            pos: 0,
            strand: b'-',
            mapq: 0,
            is_mapped: false,
            is_unique: false,
            is_linear: true,
            dist_to_5: 0,
            dist_to_3: 0,
            cigar_text: "*".into(),
            algn_ref_span: 0,
            algn_read_span: 0,
            matched_bp: 0,
            clip5_ref: 0,
            clip3_ref: 0,
            read_len: 0,
            kind: b'N',
            tags: Vec::new(),
            seq: None,
            read_side: None,
            algn_idx: None,
            same_side_count: None,
        }
    }

    /// Reset coordinates to the unmapped state (pairtools `mask_alignment`).
    pub fn mask(&mut self) {
        self.chrom = None;
        self.pos5 = 0;
        self.pos3 = 0;
        self.pos = 0;
        self.strand = b'-';
    }
}

/// Value of the `i`-th requested tag as pysam's `str(value)`.
pub fn tag_value_string(a: &Alignment, i: usize) -> String {
    a.tags.get(i).cloned().unwrap_or_default()
}

/// Leading soft clip length (pysam `query_alignment_start`).
pub fn leading_soft_clip<R: sam::alignment::Record>(r: &R) -> usize {
    let cigar = r.cigar();
    let mut n = 0usize;
    for op in cigar.iter() {
        let Ok(op) = op else { break };
        match op.kind() {
            Kind::SoftClip => n += op.len(),
            Kind::HardClip => {}
            _ => break,
        }
    }
    n
}

fn format_tag_value(v: &Value<'_>) -> String {
    match v {
        Value::Character(c) => (*c as char).to_string(),
        Value::Int8(i) => i.to_string(),
        Value::UInt8(i) => i.to_string(),
        Value::Int16(i) => i.to_string(),
        Value::UInt16(i) => i.to_string(),
        Value::Int32(i) => i.to_string(),
        Value::UInt32(i) => i.to_string(),
        Value::Float(f) => crate::util::pyfloat::format_py_float(f64::from(*f)),
        Value::String(s) | Value::Hex(s) => String::from_utf8_lossy(s).into_owned(),
        Value::Array(arr) => {
            use noodles_sam::alignment::record::data::field::value::Array;
            let (code, items): (char, Vec<String>) = match arr {
                Array::Int8(v) => ('b', v.iter().flatten().map(|x| x.to_string()).collect()),
                Array::UInt8(v) => ('B', v.iter().flatten().map(|x| x.to_string()).collect()),
                Array::Int16(v) => ('h', v.iter().flatten().map(|x| x.to_string()).collect()),
                Array::UInt16(v) => ('H', v.iter().flatten().map(|x| x.to_string()).collect()),
                Array::Int32(v) => ('i', v.iter().flatten().map(|x| x.to_string()).collect()),
                Array::UInt32(v) => ('I', v.iter().flatten().map(|x| x.to_string()).collect()),
                Array::Float(v) => (
                    'f',
                    v.iter()
                        .flatten()
                        .map(|x| crate::util::pyfloat::format_py_float(f64::from(x)))
                        .collect(),
                ),
            };
            format!("array('{code}', [{}])", items.join(", "))
        }
    }
}

/// pairtools `parse_pysam_entry`.
pub fn parse_alignment<R: sam::alignment::Record>(
    rec: &R,
    header: &sam::Header,
    min_mapq: u8,
    sam_tags: &[[u8; 2]],
    store_seq: bool,
) -> Result<Alignment> {
    let err = |msg: String| KiraError::Alignment {
        message: msg,
        location: Location::default(),
    };
    let flags = rec
        .flags()
        .map_err(|e| err(format!("bad flags: {e}")))?
        .bits();
    let is_mapped = flags & 0x04 == 0;
    let mapq = rec
        .mapping_quality()
        .transpose()
        .map_err(|e| err(format!("bad MAPQ: {e}")))?
        .map(u8::from)
        .unwrap_or(255);
    let is_unique = mapq >= min_mapq;
    let data = rec.data();
    let is_linear = data.get(&Tag::new(b'S', b'A')).is_none();
    // CIGAR dict.
    let mut matched_bp = 0u64;
    let mut algn_ref_span = 0u64;
    let mut algn_read_span = 0u64;
    let mut read_len = 0u64;
    let mut clip5_ref = 0u64;
    let mut clip3_ref = 0u64;
    let mut cigar_text = String::new();
    let cigar = rec.cigar();
    for op in cigar.iter() {
        let op = op.map_err(|e| err(format!("bad CIGAR: {e}")))?;
        let len = op.len() as u64;
        let code = match op.kind() {
            Kind::Match => b'M',
            Kind::Insertion => b'I',
            Kind::Deletion => b'D',
            Kind::Skip => b'N',
            Kind::SoftClip => b'S',
            Kind::HardClip => b'H',
            Kind::Pad => b'P',
            Kind::SequenceMatch => b'=',
            Kind::SequenceMismatch => b'X',
        };
        cigar_text.push_str(&len.to_string());
        cigar_text.push(code as char);
        match op.kind() {
            Kind::Match => {
                matched_bp += len;
                algn_ref_span += len;
                algn_read_span += len;
                read_len += len;
            }
            Kind::Insertion => {
                algn_read_span += len;
                read_len += len;
            }
            Kind::Deletion => {
                algn_ref_span += len;
            }
            Kind::SoftClip | Kind::HardClip => {
                read_len += len;
                if matched_bp == 0 {
                    clip5_ref = len;
                } else {
                    clip3_ref = len;
                }
            }
            _ => {}
        }
    }
    if cigar_text.is_empty() {
        cigar_text.push_str("None");
    }
    let (chrom, strand, pos5, pos3, dist_to_5, dist_to_3);
    if is_mapped {
        let mapped_strand = if flags & 0x10 == 0 { b'+' } else { b'-' };
        if mapped_strand == b'+' {
            dist_to_5 = clip5_ref;
            dist_to_3 = clip3_ref;
        } else {
            dist_to_5 = clip3_ref;
            dist_to_3 = clip5_ref;
        }
        if is_unique {
            let ref_id = rec
                .reference_sequence_id(header)
                .transpose()
                .map_err(|e| err(format!("bad reference id: {e}")))?;
            let start = rec
                .alignment_start()
                .transpose()
                .map_err(|e| err(format!("bad position: {e}")))?
                .map(|p| usize::from(p) as u64)
                .unwrap_or(0);
            // pysam reference_start is 0-based: start here is 1-based.
            let reference_start = start.saturating_sub(1);
            chrom = ref_id;
            strand = mapped_strand;
            if mapped_strand == b'+' {
                pos5 = reference_start + 1;
                pos3 = reference_start + algn_ref_span;
            } else {
                pos5 = reference_start + algn_ref_span;
                pos3 = reference_start + 1;
            }
        } else {
            chrom = None;
            strand = b'-';
            pos5 = 0;
            pos3 = 0;
        }
    } else {
        chrom = None;
        strand = b'-';
        pos5 = 0;
        pos3 = 0;
        dist_to_5 = 0;
        dist_to_3 = 0;
    }
    let kind = if !is_mapped {
        b'N'
    } else if !is_unique {
        b'M'
    } else {
        b'U'
    };
    let mut tags = Vec::with_capacity(sam_tags.len());
    for t in sam_tags {
        let v = data
            .get(&Tag::from(*t))
            .transpose()
            .map_err(|e| err(format!("bad tag {}{}: {e}", t[0] as char, t[1] as char)))?;
        tags.push(v.map(|v| format_tag_value(&v)).unwrap_or_default());
    }
    let seq = if store_seq {
        let s = rec.sequence();
        Some(s.iter().collect::<Vec<u8>>())
    } else {
        None
    };
    Ok(Alignment {
        chrom,
        pos5,
        pos3,
        pos: pos5,
        strand,
        mapq,
        is_mapped,
        is_unique,
        is_linear,
        dist_to_5,
        dist_to_3,
        cigar_text,
        algn_ref_span,
        algn_read_span,
        matched_bp,
        clip5_ref,
        clip3_ref,
        read_len,
        kind,
        tags,
        seq,
        read_side: None,
        algn_idx: None,
        same_side_count: None,
    })
}
