//! Private on-disk run format for external sorting.
//!
//! **This format is not stable.** It is an implementation detail of one
//! kira-pairs version and may change without notice.
//!
//! Layout:
//!
//! ```text
//! magic "KPRUN1\0\0" | version u32 | flags u32 | n_records u64
//! block*: comp_len u32 | raw_len u32 | n_records u32 | payload[comp_len]
//! terminator: comp_len = 0, raw_len = 0, n_records = 0
//! ```
//!
//! A raw block is a sequence of `[key: 48 bytes][line_len: u32][line]`.
//! Payloads are LZ4 block-compressed when flag bit 0 is set.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use crate::error::{KiraError, Result};
use crate::pairs::record::{PAIR_TYPE_INLINE, PairKey};

const MAGIC: &[u8; 8] = b"KPRUN1\0\0";
const VERSION: u32 = 1;
const FLAG_LZ4: u32 = 1;
/// Serialized size of a [`PairKey`].
pub const KEY_BYTES: usize = 48;
/// Default raw block size for runs.
pub const DEFAULT_RUN_BLOCK: usize = 1 << 20;

/// Serialise a key into 48 little-endian bytes.
#[inline]
pub fn encode_key(k: &PairKey, out: &mut Vec<u8>) {
    out.extend_from_slice(&k.seq.to_le_bytes());
    out.extend_from_slice(&k.pos1.to_le_bytes());
    out.extend_from_slice(&k.pos2.to_le_bytes());
    out.extend_from_slice(&k.chrom1.to_le_bytes());
    out.extend_from_slice(&k.chrom2.to_le_bytes());
    out.extend_from_slice(&k.pair_type);
    out.push(k.strand1);
    out.push(k.strand2);
    out.push(k.flags);
    out.push(k.pair_type_len);
    out.extend_from_slice(&[0u8; 4]);
}

/// Decode a key from 48 bytes.
#[inline]
pub fn decode_key(b: &[u8]) -> PairKey {
    let u64_at = |i: usize| {
        u64::from_le_bytes([
            b[i],
            b[i + 1],
            b[i + 2],
            b[i + 3],
            b[i + 4],
            b[i + 5],
            b[i + 6],
            b[i + 7],
        ])
    };
    let u32_at = |i: usize| u32::from_le_bytes([b[i], b[i + 1], b[i + 2], b[i + 3]]);
    let mut pt = [0u8; PAIR_TYPE_INLINE];
    pt.copy_from_slice(&b[32..40]);
    PairKey {
        seq: u64_at(0),
        pos1: u64_at(8),
        pos2: u64_at(16),
        chrom1: u32_at(24),
        chrom2: u32_at(28),
        pair_type: pt,
        strand1: b[40],
        strand2: b[41],
        flags: b[42],
        pair_type_len: b[43],
    }
}

/// Encode a raw block of records into `raw`.
pub struct BlockEncoder {
    raw: Vec<u8>,
    n: u32,
    limit: usize,
}

impl BlockEncoder {
    /// Encoder with a target raw block size.
    pub fn new(limit: usize) -> Self {
        Self {
            raw: Vec::with_capacity(limit + 4096),
            n: 0,
            limit,
        }
    }

    /// Append a record. Returns true when the block is full and should be
    /// taken with [`BlockEncoder::take`].
    #[inline]
    pub fn push(&mut self, key: &PairKey, line: &[u8]) -> bool {
        encode_key(key, &mut self.raw);
        self.raw
            .extend_from_slice(&(line.len() as u32).to_le_bytes());
        self.raw.extend_from_slice(line);
        self.n += 1;
        self.raw.len() >= self.limit
    }

    /// Take the accumulated raw block (bytes, record count).
    pub fn take(&mut self) -> (Vec<u8>, u32) {
        let raw = std::mem::replace(&mut self.raw, Vec::with_capacity(self.limit + 4096));
        let n = self.n;
        self.n = 0;
        (raw, n)
    }

    /// Whether anything is buffered.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }
}

/// Compress a raw block payload (LZ4 block format, size-prepended).
pub fn compress_payload(raw: &[u8]) -> Vec<u8> {
    lz4_flex::block::compress_prepend_size(raw)
}

/// Decompress a payload produced by [`compress_payload`].
pub fn decompress_payload(comp: &[u8], expected: usize) -> io::Result<Vec<u8>> {
    let out = lz4_flex::block::decompress_size_prepended(comp)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e.to_string()))?;
    if out.len() != expected {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "run block size mismatch",
        ));
    }
    Ok(out)
}

/// Sequential writer of a run file.
pub struct RunWriter {
    path: PathBuf,
    out: BufWriter<File>,
    compress: bool,
    n_records: u64,
    bytes: u64,
}

impl RunWriter {
    /// Create a run file.
    pub fn create(path: &Path, compress: bool) -> Result<Self> {
        let f = File::create(path).map_err(|e| KiraError::io(path, e))?;
        let mut out = BufWriter::with_capacity(1 << 20, f);
        let mut header = Vec::with_capacity(24);
        header.extend_from_slice(MAGIC);
        header.extend_from_slice(&VERSION.to_le_bytes());
        header.extend_from_slice(&(if compress { FLAG_LZ4 } else { 0 }).to_le_bytes());
        header.extend_from_slice(&0u64.to_le_bytes());
        out.write_all(&header).map_err(|e| KiraError::io(path, e))?;
        Ok(Self {
            path: path.to_path_buf(),
            out,
            compress,
            n_records: 0,
            bytes: 24,
        })
    }

    /// Whether payloads are compressed.
    pub fn compress(&self) -> bool {
        self.compress
    }

    /// Write a (possibly already compressed) block.
    pub fn write_block(&mut self, payload: &[u8], raw_len: u32, n: u32) -> Result<()> {
        let mut hdr = [0u8; 12];
        hdr[..4].copy_from_slice(&(payload.len() as u32).to_le_bytes());
        hdr[4..8].copy_from_slice(&raw_len.to_le_bytes());
        hdr[8..].copy_from_slice(&n.to_le_bytes());
        self.out
            .write_all(&hdr)
            .map_err(|e| KiraError::io(&self.path, e))?;
        self.out
            .write_all(payload)
            .map_err(|e| KiraError::io(&self.path, e))?;
        self.n_records += u64::from(n);
        self.bytes += 12 + payload.len() as u64;
        Ok(())
    }

    /// Encode and write a raw block (compressing if enabled).
    pub fn write_raw_block(&mut self, raw: &[u8], n: u32) -> Result<()> {
        if self.compress {
            let c = compress_payload(raw);
            self.write_block(&c, raw.len() as u32, n)
        } else {
            self.write_block(raw, raw.len() as u32, n)
        }
    }

    /// Write the terminator, patch the record count and close.
    pub fn finish(mut self) -> Result<(PathBuf, u64, u64)> {
        self.out
            .write_all(&[0u8; 12])
            .map_err(|e| KiraError::io(&self.path, e))?;
        self.bytes += 12;
        self.out.flush().map_err(|e| KiraError::io(&self.path, e))?;
        let mut f = self
            .out
            .into_inner()
            .map_err(|e| KiraError::io(&self.path, e.into_error()))?;
        use std::io::Seek;
        f.seek(io::SeekFrom::Start(16))
            .map_err(|e| KiraError::io(&self.path, e))?;
        f.write_all(&self.n_records.to_le_bytes())
            .map_err(|e| KiraError::io(&self.path, e))?;
        f.sync_data().ok();
        Ok((self.path, self.n_records, self.bytes))
    }
}

/// Sequential reader of a run file that hands out one record at a time.
pub struct RunReader {
    path: PathBuf,
    inp: BufReader<File>,
    compress: bool,
    n_records: u64,
    seen: u64,
    block: Vec<u8>,
    pos: usize,
    block_remaining: u32,
    current_key: PairKey,
    current_line: (usize, usize),
    done: bool,
}

impl RunReader {
    /// Open a run file and position on the first record.
    pub fn open(path: &Path) -> Result<Self> {
        let f = File::open(path).map_err(|e| KiraError::io(path, e))?;
        let mut inp = BufReader::with_capacity(1 << 20, f);
        let mut hdr = [0u8; 24];
        inp.read_exact(&mut hdr).map_err(|e| KiraError::RunFile {
            path: path.to_path_buf(),
            message: format!("cannot read header: {e}"),
        })?;
        if &hdr[..8] != MAGIC {
            return Err(KiraError::RunFile {
                path: path.to_path_buf(),
                message: "bad magic".into(),
            });
        }
        let version = u32::from_le_bytes([hdr[8], hdr[9], hdr[10], hdr[11]]);
        if version != VERSION {
            return Err(KiraError::RunFile {
                path: path.to_path_buf(),
                message: format!("unsupported run version {version}"),
            });
        }
        let flags = u32::from_le_bytes([hdr[12], hdr[13], hdr[14], hdr[15]]);
        let n_records = u64::from_le_bytes([
            hdr[16], hdr[17], hdr[18], hdr[19], hdr[20], hdr[21], hdr[22], hdr[23],
        ]);
        let mut r = Self {
            path: path.to_path_buf(),
            inp,
            compress: flags & FLAG_LZ4 != 0,
            n_records,
            seen: 0,
            block: Vec::new(),
            pos: 0,
            block_remaining: 0,
            current_key: decode_key(&[0u8; KEY_BYTES]),
            current_line: (0, 0),
            done: false,
        };
        r.advance()?;
        Ok(r)
    }

    /// Records declared in the header.
    pub fn n_records(&self) -> u64 {
        self.n_records
    }

    /// True when all records have been consumed.
    #[inline]
    pub fn is_done(&self) -> bool {
        self.done
    }

    /// Current key (valid when not done).
    #[inline]
    pub fn key(&self) -> &PairKey {
        &self.current_key
    }

    /// Current line (valid when not done).
    #[inline]
    pub fn line(&self) -> &[u8] {
        &self.block[self.current_line.0..self.current_line.1]
    }

    fn corrupt(&self, msg: &str) -> KiraError {
        KiraError::RunFile {
            path: self.path.clone(),
            message: msg.to_string(),
        }
    }

    fn load_block(&mut self) -> Result<bool> {
        let mut hdr = [0u8; 12];
        self.inp
            .read_exact(&mut hdr)
            .map_err(|e| self.corrupt(&format!("truncated block header: {e}")))?;
        let comp_len = u32::from_le_bytes([hdr[0], hdr[1], hdr[2], hdr[3]]) as usize;
        let raw_len = u32::from_le_bytes([hdr[4], hdr[5], hdr[6], hdr[7]]) as usize;
        let n = u32::from_le_bytes([hdr[8], hdr[9], hdr[10], hdr[11]]);
        if comp_len == 0 && raw_len == 0 && n == 0 {
            return Ok(false);
        }
        let mut payload = vec![0u8; comp_len];
        self.inp
            .read_exact(&mut payload)
            .map_err(|e| self.corrupt(&format!("truncated block payload: {e}")))?;
        self.block = if self.compress {
            decompress_payload(&payload, raw_len).map_err(|e| self.corrupt(&e.to_string()))?
        } else {
            if payload.len() != raw_len {
                return Err(self.corrupt("block length mismatch"));
            }
            payload
        };
        self.pos = 0;
        self.block_remaining = n;
        Ok(true)
    }

    /// Move to the next record. Returns false at end of run.
    pub fn advance(&mut self) -> Result<bool> {
        if self.done {
            return Ok(false);
        }
        while self.block_remaining == 0 {
            if !self.load_block()? {
                self.done = true;
                if self.seen != self.n_records {
                    return Err(self.corrupt(&format!(
                        "record count mismatch: header says {}, found {}",
                        self.n_records, self.seen
                    )));
                }
                return Ok(false);
            }
        }
        let p = self.pos;
        if p + KEY_BYTES + 4 > self.block.len() {
            return Err(self.corrupt("record header beyond block"));
        }
        self.current_key = decode_key(&self.block[p..p + KEY_BYTES]);
        let len = u32::from_le_bytes([
            self.block[p + KEY_BYTES],
            self.block[p + KEY_BYTES + 1],
            self.block[p + KEY_BYTES + 2],
            self.block[p + KEY_BYTES + 3],
        ]) as usize;
        let start = p + KEY_BYTES + 4;
        let end = start + len;
        if end > self.block.len() {
            return Err(self.corrupt("record line beyond block"));
        }
        self.current_line = (start, end);
        self.pos = end;
        self.block_remaining -= 1;
        self.seen += 1;
        Ok(true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn key(i: u64) -> PairKey {
        PairKey {
            seq: i,
            pos1: i * 3,
            pos2: i * 7,
            chrom1: (i % 5) as u32,
            chrom2: (i % 3) as u32,
            pair_type: *b"UU\0\0\0\0\0\0",
            strand1: b'+',
            strand2: b'-',
            flags: 0,
            pair_type_len: 2,
        }
    }

    #[test]
    fn key_roundtrip() {
        let k = key(42);
        let mut v = Vec::new();
        encode_key(&k, &mut v);
        assert_eq!(v.len(), KEY_BYTES);
        assert_eq!(decode_key(&v), k);
    }

    #[test]
    fn run_roundtrip() {
        for compress in [false, true] {
            let dir = tempfile::tempdir().unwrap();
            let p = dir.path().join("r.kprun");
            let mut w = RunWriter::create(&p, compress).unwrap();
            let mut enc = BlockEncoder::new(4096);
            let n = 5000u64;
            for i in 0..n {
                let line = format!("r{i}\tchr1\t{}\tchr1\t{}\t+\t-\tUU", i * 3, i * 7);
                if enc.push(&key(i), line.as_bytes()) {
                    let (raw, c) = enc.take();
                    w.write_raw_block(&raw, c).unwrap();
                }
            }
            if !enc.is_empty() {
                let (raw, c) = enc.take();
                w.write_raw_block(&raw, c).unwrap();
            }
            let (path, count, _bytes) = w.finish().unwrap();
            assert_eq!(count, n);
            let mut r = RunReader::open(&path).unwrap();
            assert_eq!(r.n_records(), n);
            let mut i = 0u64;
            while !r.is_done() {
                assert_eq!(r.key(), &key(i));
                assert!(r.line().starts_with(format!("r{i}\t").as_bytes()));
                i += 1;
                r.advance().unwrap();
            }
            assert_eq!(i, n);
        }
    }

    #[test]
    fn detects_corruption() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("r.kprun");
        let mut w = RunWriter::create(&p, true).unwrap();
        let mut enc = BlockEncoder::new(1 << 20);
        for i in 0..100u64 {
            enc.push(&key(i), b"line");
        }
        let (raw, c) = enc.take();
        w.write_raw_block(&raw, c).unwrap();
        let (path, _, _) = w.finish().unwrap();
        let mut data = std::fs::read(&path).unwrap();
        let n = data.len();
        data.truncate(n - 20);
        std::fs::write(&path, &data).unwrap();
        let r = RunReader::open(&path);
        assert!(
            r.is_err() || {
                let mut r = r.unwrap();
                let mut failed = false;
                while !r.is_done() {
                    if r.advance().is_err() {
                        failed = true;
                        break;
                    }
                }
                failed
            }
        );
        std::fs::write(&path, b"not a run file at all, definitely not").unwrap();
        assert!(RunReader::open(&path).is_err());
    }
}
