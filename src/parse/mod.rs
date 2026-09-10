//! Standard paired-end Hi-C parsing of SAM/BAM alignments into `.pairs`
//! (pairtools `parse` semantics for the common BWA-MEM case).

pub mod alignment;
pub mod hic;

use std::io::{self, BufReader, Read, Write};
use std::num::NonZero;
use std::path::Path;
use std::sync::Arc;

use noodles_bam as bam;
use noodles_bgzf as bgzf;
use noodles_sam as sam;

use crate::chroms::{ChromDict, ChromOrder, ChromSizes, UNMAPPED_CHROM};
use crate::error::{KiraError, Location, Result};
use crate::io::compression::{Compression, open_input};
use crate::pairs::header::Header;
use crate::pairs::record::{PAIR_TYPE_INLINE, PairKey};
use crate::stats::StatsAccumulator;
use alignment::{Alignment, parse_alignment, tag_value_string};
use hic::{PairIndex, WalksPolicy, check_pair_order, parse_read};

/// pairtools' separator inside SAM columns.
pub const SAM_SEP: u8 = 0x19;
/// pairtools' separator between SAM entries of one side.
pub const INTER_SAM_SEP: &[u8] = b"\x19NEXT_SAM\x19";

/// Known `--add-columns` names besides two-letter SAM tags.
pub const EXTRA_COLUMNS: [&str; 15] = [
    "mapq",
    "pos5",
    "pos3",
    "cigar",
    "read_len",
    "matched_bp",
    "algn_ref_span",
    "algn_read_span",
    "dist_to_5",
    "dist_to_3",
    "seq",
    "mismatches",
    "read_side",
    "algn_idx",
    "same_side_algn_count",
];

/// Parser configuration (mirrors `pairtools parse` options).
#[derive(Debug, Clone)]
pub struct ParseConfig {
    /// Minimal MAPQ for a uniquely mapped alignment.
    pub min_mapq: u8,
    /// Maximal Hi-C molecule size for rescuing single ligations.
    pub max_molecule_size: u64,
    /// Gaps longer than this become null alignments.
    pub max_inter_align_gap: Option<u64>,
    /// Walks policy.
    pub walks_policy: WalksPolicy,
    /// Report the 3' alignment end instead of the 5' end.
    pub report_3_end: bool,
    /// Flip pairs into upper-triangular order.
    pub flip: bool,
    /// Replace read IDs with `.`.
    pub drop_readid: bool,
    /// Blank sequences/qualities in SAM columns.
    pub drop_seq: bool,
    /// Omit `sam1`/`sam2` columns.
    pub drop_sam: bool,
    /// Add `walk_pair_index`/`walk_pair_type` columns.
    pub add_pair_index: bool,
    /// Extra columns.
    pub add_columns: Vec<String>,
    /// Genome assembly name for the header.
    pub assembly: Option<String>,
}

impl Default for ParseConfig {
    fn default() -> Self {
        Self {
            min_mapq: 1,
            max_molecule_size: 750,
            max_inter_align_gap: Some(20),
            walks_policy: WalksPolicy::FiveUnique,
            report_3_end: false,
            flip: true,
            drop_readid: false,
            drop_seq: false,
            drop_sam: false,
            add_pair_index: false,
            add_columns: Vec::new(),
            assembly: None,
        }
    }
}

impl ParseConfig {
    /// Validate `add_columns` names.
    pub fn validate(&self) -> Result<()> {
        for c in &self.add_columns {
            let is_tag = c.len() == 2 && c.bytes().all(|b| b.is_ascii_uppercase());
            if !(EXTRA_COLUMNS.contains(&c.as_str()) || is_tag) {
                return Err(KiraError::arg(format!("{c} is not a valid extra column")));
            }
            if c == "mismatches" {
                return Err(KiraError::Unsupported(
                    "--add-columns mismatches is not implemented in this version".into(),
                ));
            }
        }
        if self.walks_policy == WalksPolicy::All {
            return Err(KiraError::Unsupported(
                "--walks-policy all (complex walk parsing) is not implemented in this version"
                    .into(),
            ));
        }
        Ok(())
    }

    /// Output column names.
    pub fn columns(&self) -> Vec<String> {
        let mut cols: Vec<String> = crate::pairs::header::STANDARD_COLUMNS
            .iter()
            .map(|s| s.to_string())
            .collect();
        if !self.drop_sam {
            cols.push("sam1".into());
            cols.push("sam2".into());
        }
        if self.add_pair_index {
            cols.push("walk_pair_index".into());
            cols.push("walk_pair_type".into());
        }
        for c in &self.add_columns {
            cols.push(format!("{c}1"));
            cols.push(format!("{c}2"));
        }
        cols
    }
}

/// Alignment input (BAM or SAM text).
pub enum AlignmentSource {
    /// BAM via a multithreaded BGZF reader.
    Bam(bam::io::Reader<bgzf::io::MultithreadedReader<Box<dyn Read + Send>>>),
    /// SAM text.
    Sam(sam::io::Reader<BufReader<Box<dyn Read + Send>>>),
}

/// Open a SAM/BAM input and read its header.
pub fn open_alignments(
    path: Option<&Path>,
    io_threads: usize,
) -> Result<(AlignmentSource, sam::Header, String)> {
    let src = open_input_raw(path)?;
    let name = src.0;
    let compression = src.1;
    let raw = src.2;
    match compression {
        Compression::Bgzf | Compression::Gzip => {
            let workers = NonZero::new(io_threads.max(1)).unwrap_or(NonZero::<usize>::MIN);
            let inner = bgzf::io::MultithreadedReader::with_worker_count(workers, raw);
            let mut reader = bam::io::Reader::from(inner);
            let header = reader.read_header().map_err(|e| KiraError::Alignment {
                message: format!("cannot read BAM header: {e}"),
                location: Location::file(&name),
            })?;
            Ok((AlignmentSource::Bam(reader), header, name))
        }
        _ => {
            let mut reader = sam::io::Reader::new(BufReader::with_capacity(1 << 20, raw));
            let header = reader.read_header().map_err(|e| KiraError::Alignment {
                message: format!("cannot read SAM header: {e}"),
                location: Location::file(&name),
            })?;
            Ok((AlignmentSource::Sam(reader), header, name))
        }
    }
}

fn open_input_raw(path: Option<&Path>) -> Result<(String, Compression, Box<dyn Read + Send>)> {
    // We need the raw (still compressed) stream for BAM, so peek ourselves.
    let (name, mut raw): (String, Box<dyn Read + Send>) = match path {
        None => ("-".into(), Box::new(io::stdin())),
        Some(p) if p.as_os_str() == "-" => ("-".into(), Box::new(io::stdin())),
        Some(p) => {
            let f = std::fs::File::open(p).map_err(|e| KiraError::io(p, e))?;
            (p.display().to_string(), Box::new(f))
        }
    };
    let mut head = vec![0u8; 18];
    let mut n = 0;
    while n < head.len() {
        match raw.read(&mut head[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(KiraError::io(Path::new(&name), e)),
        }
    }
    head.truncate(n);
    let compression = Compression::detect(&head);
    let chained: Box<dyn Read + Send> = Box::new(io::Cursor::new(head).chain(raw));
    let _ = open_input; // plain-text SAM could also be lz4; not supported for alignments
    Ok((name, compression, chained))
}

/// The Hi-C parser state.
pub struct HicParser {
    cfg: ParseConfig,
    sam_header: sam::Header,
    ref_names: Vec<Vec<u8>>,
    ref_enum: Vec<u32>,
    ref_dict_ids: Vec<u32>,
    unmapped_id: u32,
    dict: Arc<ChromDict>,
    header: Header,
    sam_tags: Vec<[u8; 2]>,
    store_seq: bool,
    seq: u64,
    sam_writer: sam::io::Writer<Vec<u8>>,
    line: Vec<u8>,
    scratch1: Vec<Alignment>,
    scratch2: Vec<Alignment>,
    /// Records processed.
    pub records_in: u64,
    /// Pairs written.
    pub pairs_out: u64,
}

impl HicParser {
    /// Build a parser from the SAM header and chromosome order file.
    pub fn new(
        cfg: ParseConfig,
        sam_header: sam::Header,
        chroms: &ChromSizes,
        dict: Arc<ChromDict>,
    ) -> Result<Self> {
        cfg.validate()?;
        let ref_names: Vec<Vec<u8>> = sam_header
            .reference_sequences()
            .keys()
            .map(|k| k.to_vec())
            .collect();
        let ref_lengths: Vec<u64> = sam_header
            .reference_sequences()
            .values()
            .map(|m| m.length().get() as u64)
            .collect();
        if ref_names.is_empty() {
            return Err(KiraError::Alignment {
                message: "the input SAM/BAM header has no @SQ reference sequences".into(),
                location: Location::default(),
            });
        }
        let ref_name_strs: Vec<String> = ref_names
            .iter()
            .map(|n| String::from_utf8_lossy(n).into_owned())
            .collect();
        let order = ChromOrder::from_names_restricted(
            chroms.names().iter().map(String::as_str),
            ref_name_strs.iter().map(String::as_str),
        );
        let ref_enum: Vec<u32> = ref_names
            .iter()
            .map(|n| order.get(n).unwrap_or(u32::MAX))
            .collect();
        let ref_dict_ids: Vec<u32> = ref_names.iter().map(|n| dict.intern(n)).collect();
        let unmapped_id = dict.intern(UNMAPPED_CHROM);
        // Header: chromsizes in enumeration order.
        let mut cs = ChromSizes::default();
        for name in order.names() {
            if let Some(i) = ref_name_strs.iter().position(|r| r == name) {
                cs.push(name, ref_lengths[i]);
            }
        }
        let shape = if cfg.flip {
            "upper triangle"
        } else {
            "whole matrix"
        };
        let mut header =
            Header::standard(cfg.assembly.as_deref(), Some(&cs), &cfg.columns(), shape);
        let mut buf = Vec::new();
        {
            let mut w = sam::io::Writer::new(&mut buf);
            w.write_header(&sam_header)
                .map_err(|e| KiraError::Alignment {
                    message: format!("cannot serialise SAM header: {e}"),
                    location: Location::default(),
                })?;
        }
        let sam_lines: Vec<String> = String::from_utf8_lossy(&buf)
            .lines()
            .filter(|l| !l.is_empty())
            .map(String::from)
            .collect();
        header.insert_samheader(&sam_lines);
        let sam_tags: Vec<[u8; 2]> = cfg
            .add_columns
            .iter()
            .filter(|c| c.len() == 2 && c.bytes().all(|b| b.is_ascii_uppercase()))
            .map(|c| [c.as_bytes()[0], c.as_bytes()[1]])
            .collect();
        let store_seq = cfg.add_columns.iter().any(|c| c == "seq");
        Ok(Self {
            cfg,
            sam_header,
            ref_names,
            ref_enum,
            ref_dict_ids,
            unmapped_id,
            dict,
            header,
            sam_tags,
            store_seq,
            seq: 0,
            sam_writer: sam::io::Writer::new(Vec::new()),
            line: Vec::with_capacity(512),
            scratch1: Vec::new(),
            scratch2: Vec::new(),
            records_in: 0,
            pairs_out: 0,
        })
    }

    /// The `.pairs` header for the output (without the kira-pairs `@PG`).
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Mutable access to the header (to append `@PG`).
    pub fn header_mut(&mut self) -> &mut Header {
        &mut self.header
    }

    /// Chromosome dictionary.
    pub fn dict(&self) -> &Arc<ChromDict> {
        &self.dict
    }

    /// Parse one read (all records sharing a query name) and emit pairs.
    pub fn parse_group<R, F>(
        &mut self,
        name: &[u8],
        records: &mut [R],
        mut sink: F,
        stats: Option<&mut StatsAccumulator>,
    ) -> Result<()>
    where
        R: sam::alignment::Record,
        F: FnMut(&PairKey, &[u8]) -> Result<()>,
    {
        self.records_in += records.len() as u64;
        // pairtools: sort by (is_read2, query_alignment_start), stable.
        let mut order: Vec<usize> = (0..records.len()).collect();
        let keys: Vec<(bool, usize)> = records
            .iter()
            .map(|r| {
                let flags = r.flags().map(|f| f.bits()).unwrap_or(0);
                (flags & 0x80 != 0, alignment::leading_soft_clip(r))
            })
            .collect();
        order.sort_by_key(|&i| keys[i]);
        let mut sams1: Vec<usize> = Vec::new();
        let mut sams2: Vec<usize> = Vec::new();
        for &i in &order {
            let flags = records[i].flags().map(|f| f.bits()).unwrap_or(0);
            if flags & 0x40 != 0 {
                sams1.push(i);
            } else {
                sams2.push(i);
            }
        }
        let mut algns1 = std::mem::take(&mut self.scratch1);
        let mut algns2 = std::mem::take(&mut self.scratch2);
        algns1.clear();
        algns2.clear();
        let is_empty = sams1.is_empty() || sams2.is_empty();
        let (mut hic1, mut hic2, pair_index): (Alignment, Alignment, PairIndex);
        if is_empty {
            let mut a = Alignment::empty();
            a.kind = b'X';
            let mut b = Alignment::empty();
            b.kind = b'X';
            hic1 = a;
            hic2 = b;
            pair_index = PairIndex::R12;
        } else {
            for &i in &sams1 {
                algns1.push(parse_alignment(
                    &records[i],
                    &self.sam_header,
                    self.cfg.min_mapq,
                    &self.sam_tags,
                    self.store_seq,
                )?);
            }
            for &i in &sams2 {
                algns2.push(parse_alignment(
                    &records[i],
                    &self.sam_header,
                    self.cfg.min_mapq,
                    &self.sam_tags,
                    self.store_seq,
                )?);
            }
            let (h1, h2, pi) = parse_read(
                &mut algns1,
                &mut algns2,
                self.cfg.max_molecule_size,
                self.cfg.max_inter_align_gap,
                self.cfg.walks_policy,
            );
            hic1 = h1;
            hic2 = h2;
            pair_index = pi;
        }
        hic1.pos = if self.cfg.report_3_end {
            hic1.pos3
        } else {
            hic1.pos5
        };
        hic2.pos = if self.cfg.report_3_end {
            hic2.pos3
        } else {
            hic2.pos5
        };
        let (mut side1, mut side2) = (&sams1, &sams2);
        if self.cfg.flip && !check_pair_order(&hic1, &hic2, &self.ref_enum) {
            std::mem::swap(&mut hic1, &mut hic2);
            std::mem::swap(&mut side1, &mut side2);
        }
        let pair_type = [hic1.kind, hic2.kind];
        // Build the line.
        self.line.clear();
        if self.cfg.drop_readid {
            self.line.push(b'.');
        } else {
            self.line.extend_from_slice(name);
        }
        self.line.push(b'\t');
        self.push_chrom(hic1.chrom);
        self.line.push(b'\t');
        crate::util::int::write_u64(&mut self.line, hic1.pos);
        self.line.push(b'\t');
        self.push_chrom(hic2.chrom);
        self.line.push(b'\t');
        crate::util::int::write_u64(&mut self.line, hic2.pos);
        self.line.push(b'\t');
        self.line.push(hic1.strand);
        self.line.push(b'\t');
        self.line.push(hic2.strand);
        self.line.push(b'\t');
        self.line.extend_from_slice(&pair_type);
        if !self.cfg.drop_sam {
            for side in [side1, side2] {
                self.line.push(b'\t');
                for (k, &i) in side.iter().enumerate() {
                    if k > 0 {
                        self.line.extend_from_slice(INTER_SAM_SEP);
                    }
                    self.push_sam(&records[i], &pair_type)?;
                }
            }
        }
        if self.cfg.add_pair_index {
            self.line.push(b'\t');
            self.line.push(b'1');
            self.line.push(b'\t');
            self.line.extend_from_slice(pair_index.label().as_bytes());
        }
        for (ci, col) in self.cfg.add_columns.iter().enumerate() {
            for a in [&hic1, &hic2] {
                self.line.push(b'\t');
                push_extra_column(&mut self.line, a, col, &self.sam_tags, ci);
            }
        }
        let key = PairKey {
            seq: self.seq,
            pos1: hic1.pos,
            pos2: hic2.pos,
            chrom1: self.dict_id(hic1.chrom),
            chrom2: self.dict_id(hic2.chrom),
            pair_type: {
                let mut p = [0u8; PAIR_TYPE_INLINE];
                p[0] = pair_type[0];
                p[1] = pair_type[1];
                p
            },
            strand1: hic1.strand,
            strand2: hic2.strand,
            flags: 0,
            pair_type_len: 2,
        };
        self.seq += 1;
        self.pairs_out += 1;
        if let Some(s) = stats {
            s.observe(&key, None, false);
        }
        sink(&key, &self.line)?;
        self.scratch1 = algns1;
        self.scratch2 = algns2;
        Ok(())
    }

    #[inline]
    fn dict_id(&self, chrom: Option<usize>) -> u32 {
        match chrom {
            Some(i) => self.ref_dict_ids[i],
            None => self.unmapped_id,
        }
    }

    #[inline]
    fn push_chrom(&mut self, chrom: Option<usize>) {
        match chrom {
            Some(i) => self.line.extend_from_slice(&self.ref_names[i]),
            None => self.line.extend_from_slice(UNMAPPED_CHROM),
        }
    }

    fn push_sam<R: sam::alignment::Record>(&mut self, rec: &R, pair_type: &[u8; 2]) -> Result<()> {
        self.sam_writer.get_mut().clear();
        use sam::alignment::io::Write as _;
        self.sam_writer
            .write_alignment_record(&self.sam_header, rec)
            .map_err(|e| KiraError::Alignment {
                message: format!("cannot serialise SAM record: {e}"),
                location: Location::default(),
            })?;
        let buf = self.sam_writer.get_ref();
        let text = buf.strip_suffix(b"\n").unwrap_or(buf);
        if self.cfg.drop_seq {
            let mut field = 0;
            for b in text {
                if *b == b'\t' {
                    field += 1;
                    self.line.push(SAM_SEP);
                    if field == 9 || field == 10 {
                        self.line.push(b'*');
                    }
                } else if field == 9 || field == 10 {
                    continue;
                } else {
                    self.line.push(*b);
                }
            }
        } else {
            for b in text {
                self.line.push(if *b == b'\t' { SAM_SEP } else { *b });
            }
        }
        self.line.push(SAM_SEP);
        self.line.extend_from_slice(b"Yt:Z:");
        self.line.extend_from_slice(pair_type);
        Ok(())
    }
}

fn push_extra_column(
    line: &mut Vec<u8>,
    a: &Alignment,
    col: &str,
    sam_tags: &[[u8; 2]],
    _ci: usize,
) {
    use crate::util::int::write_u64;
    match col {
        "mapq" => write_u64(line, u64::from(a.mapq)),
        "pos5" => write_u64(line, a.pos5),
        "pos3" => write_u64(line, a.pos3),
        "cigar" => line.extend_from_slice(a.cigar_text.as_bytes()),
        "read_len" => write_u64(line, a.read_len),
        "matched_bp" => write_u64(line, a.matched_bp),
        "algn_ref_span" => write_u64(line, a.algn_ref_span),
        "algn_read_span" => write_u64(line, a.algn_read_span),
        "dist_to_5" => write_u64(line, a.dist_to_5),
        "dist_to_3" => write_u64(line, a.dist_to_3),
        "seq" => {
            if let Some(s) = &a.seq {
                line.extend_from_slice(s);
            }
        }
        "mismatches" => {}
        "read_side" => {
            if let Some(s) = a.read_side {
                write_u64(line, s);
            }
        }
        "algn_idx" => {
            if let Some(s) = a.algn_idx {
                write_u64(line, s);
            }
        }
        "same_side_algn_count" => {
            if let Some(s) = a.same_side_count {
                write_u64(line, s);
            }
        }
        tag => {
            if let Some(i) = sam_tags.iter().position(|t| t == tag.as_bytes()) {
                line.extend_from_slice(tag_value_string(a, i).as_bytes());
            }
        }
    }
}

/// Drive a parser over a whole alignment source, grouping by read name.
pub fn drive<F>(
    source: AlignmentSource,
    parser: &mut HicParser,
    mut stats: Option<&mut StatsAccumulator>,
    name: &str,
    mut sink: F,
) -> Result<()>
where
    F: FnMut(&PairKey, &[u8]) -> Result<()>,
{
    match source {
        AlignmentSource::Bam(mut reader) => {
            let mut group: Vec<bam::Record> = Vec::new();
            let mut spare = bam::Record::default();
            let mut current: Vec<u8> = Vec::new();
            let mut n = 0u64;
            loop {
                let got = reader
                    .read_record(&mut spare)
                    .map_err(|e| KiraError::Alignment {
                        message: format!("corrupt BAM record: {e}"),
                        location: Location::file(name).at_line(n + 1),
                    })?;
                if got == 0 {
                    break;
                }
                n += 1;
                let qname = spare.name().map(|s| s.to_vec()).unwrap_or_default();
                if !group.is_empty() && qname != current {
                    parser.parse_group(&current, &mut group, &mut sink, stats.as_deref_mut())?;
                    group.clear();
                }
                current = qname;
                group.push(std::mem::take(&mut spare));
            }
            if !group.is_empty() {
                parser.parse_group(&current, &mut group, &mut sink, stats)?;
            }
        }
        AlignmentSource::Sam(mut reader) => {
            let mut group: Vec<sam::Record> = Vec::new();
            let mut spare = sam::Record::default();
            let mut current: Vec<u8> = Vec::new();
            let mut n = 0u64;
            loop {
                let got = reader
                    .read_record(&mut spare)
                    .map_err(|e| KiraError::Alignment {
                        message: format!("malformed SAM record: {e}"),
                        location: Location::file(name).at_line(n + 1),
                    })?;
                if got == 0 {
                    break;
                }
                n += 1;
                let qname = spare.name().map(|s| s.to_vec()).unwrap_or_default();
                if !group.is_empty() && qname != current {
                    parser.parse_group(&current, &mut group, &mut sink, stats.as_deref_mut())?;
                    group.clear();
                }
                current = qname;
                group.push(std::mem::take(&mut spare));
            }
            if !group.is_empty() {
                parser.parse_group(&current, &mut group, &mut sink, stats)?;
            }
        }
    }
    Ok(())
}

/// Write pairs lines produced by [`drive`] to a writer (helper for tests).
pub fn write_line<W: Write>(w: &mut W, line: &[u8]) -> io::Result<()> {
    w.write_all(line)?;
    w.write_all(b"\n")
}
