//! Streaming `.pairs` reader: header first, then large newline-aligned
//! body blocks suitable for parallel parsing.

use std::io::{self, Read};
use std::path::Path;

use crate::error::{KiraError, Location, Result};
use crate::io::buffered::{BlockReader, DEFAULT_BLOCK_SIZE, LineBlock};
use crate::io::compression::{Compression, open_input};
use crate::pairs::columns::ColumnMap;
use crate::pairs::header::Header;

/// Iterator over body blocks of a `.pairs` stream.
pub struct BodyBlocks {
    inner: BlockReader<Box<dyn Read + Send>>,
    name: String,
}

impl BodyBlocks {
    /// Next block of complete lines.
    pub fn next_block(&mut self) -> Result<Option<LineBlock>> {
        self.inner
            .next_block()
            .map_err(|e| KiraError::io(Path::new(&self.name), e))
    }

    /// Bytes consumed from the (decompressed) input so far.
    pub fn bytes_read(&self) -> u64 {
        self.inner.bytes_read()
    }

    /// Input name for diagnostics.
    pub fn name(&self) -> &str {
        &self.name
    }
}

impl Iterator for BodyBlocks {
    type Item = Result<LineBlock>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_block().transpose()
    }
}

/// A `.pairs` input with its parsed header.
pub struct PairsReader {
    header: Header,
    columns: ColumnMap,
    name: String,
    compression: Compression,
    body: Option<BodyBlocks>,
}

impl PairsReader {
    /// Open a path (`None`/`-` = stdin) with automatic decompression.
    pub fn open(path: Option<&Path>, threads: usize) -> Result<Self> {
        Self::open_with_block_size(path, threads, DEFAULT_BLOCK_SIZE)
    }

    /// Open with an explicit body block size (see
    /// [`crate::memory::MemoryBudget::block_size`]).
    pub fn open_with_block_size(
        path: Option<&Path>,
        threads: usize,
        block_size: usize,
    ) -> Result<Self> {
        let src = open_input(path, threads)?;
        Self::from_reader(src.reader, &src.name, src.compression, block_size)
    }

    /// Wrap an arbitrary reader.
    pub fn from_reader(
        mut reader: Box<dyn Read + Send>,
        name: &str,
        compression: Compression,
        block_size: usize,
    ) -> Result<Self> {
        let (header_lines, rest, n_header_lines) = read_header(&mut reader, name)?;
        let header = Header::from_lines(header_lines)?;
        if header.is_empty() {
            log::warn!(
                "{name}: headerless input, assuming standard columns (readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type)"
            );
        }
        let columns = ColumnMap::from_names(header.columns()).map_err(|e| match e {
            KiraError::MissingColumn { column, .. } => KiraError::MissingColumn {
                column,
                location: Location::file(name),
            },
            other => other,
        })?;
        let inner = BlockReader::with_prefix(reader, block_size, rest, n_header_lines + 1);
        Ok(Self {
            header,
            columns,
            name: name.to_string(),
            compression,
            body: Some(BodyBlocks {
                inner,
                name: name.to_string(),
            }),
        })
    }

    /// The parsed header.
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// Column map derived from the header.
    pub fn columns(&self) -> &ColumnMap {
        &self.columns
    }

    /// Input name for diagnostics.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Detected input compression.
    pub fn compression(&self) -> Compression {
        self.compression
    }

    /// Take ownership of the body block iterator.
    pub fn into_body(mut self) -> BodyBlocks {
        // `body` is only ever taken once; construct guarantees presence.
        self.body
            .take()
            .unwrap_or_else(|| unreachable!("body already taken"))
    }

    /// Split into header, columns and body.
    pub fn into_parts(mut self) -> (Header, ColumnMap, BodyBlocks) {
        let body = self
            .body
            .take()
            .unwrap_or_else(|| unreachable!("body already taken"));
        (self.header, self.columns, body)
    }

    /// Read the next body block without consuming the reader.
    pub fn next_block(&mut self) -> Result<Option<LineBlock>> {
        match self.body.as_mut() {
            Some(b) => b.next_block(),
            None => Ok(None),
        }
    }
}

/// Read header lines from a stream. Returns the header lines, the leftover
/// bytes that belong to the body and the number of header lines.
fn read_header(reader: &mut dyn Read, name: &str) -> Result<(Vec<String>, Vec<u8>, u64)> {
    let mut buf: Vec<u8> = Vec::with_capacity(64 * 1024);
    let mut lines = Vec::new();
    let mut scan_from = 0usize;
    let mut eof = false;
    loop {
        // Find complete lines in buf[scan_from..].
        while let Some(nl) = memchr::memchr(b'\n', &buf[scan_from..]) {
            let end = scan_from + nl;
            let line = &buf[scan_from..end];
            if line.first() == Some(&b'#') {
                lines.push(decode_line(line, name, lines.len() as u64 + 1)?);
                scan_from = end + 1;
            } else {
                let rest = buf.split_off(scan_from);
                let n = lines.len() as u64;
                return Ok((lines, rest, n));
            }
        }
        // Partial line at buf[scan_from..]: if it does not start with '#',
        // it is body (or empty), stop here.
        if buf.len() > scan_from && buf[scan_from] != b'#' {
            let rest = buf.split_off(scan_from);
            let n = lines.len() as u64;
            return Ok((lines, rest, n));
        }
        if eof {
            if buf.len() > scan_from {
                // Unterminated final header line.
                let line = &buf[scan_from..];
                lines.push(decode_line(line, name, lines.len() as u64 + 1)?);
            }
            let n = lines.len() as u64;
            return Ok((lines, Vec::new(), n));
        }
        // Compact consumed header bytes and read more.
        if scan_from > 0 {
            buf.drain(..scan_from);
            scan_from = 0;
        }
        let old = buf.len();
        buf.resize(old + 64 * 1024, 0);
        match reader.read(&mut buf[old..]) {
            Ok(0) => {
                buf.truncate(old);
                eof = true;
            }
            Ok(n) => buf.truncate(old + n),
            Err(e) if e.kind() == io::ErrorKind::Interrupted => buf.truncate(old),
            Err(e) => return Err(KiraError::io(Path::new(name), e)),
        }
    }
}

fn decode_line(line: &[u8], name: &str, lineno: u64) -> Result<String> {
    let line = line.strip_suffix(b"\r").unwrap_or(line);
    String::from_utf8(line.to_vec()).map_err(|_| KiraError::Header {
        message: "header line is not valid UTF-8".into(),
        location: Location::file(name).at_line(lineno),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn open(data: &[u8]) -> PairsReader {
        let r: Box<dyn Read + Send> = Box::new(io::Cursor::new(data.to_vec()));
        PairsReader::from_reader(r, "test", Compression::None, 64 * 1024).unwrap()
    }

    #[test]
    fn splits_header_and_body() {
        let data = b"## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\nr1\tchr1\t1\tchr1\t2\t+\t-\tUU\nr2\tchr1\t3\tchr1\t4\t+\t-\tUU";
        let r = open(data);
        assert_eq!(r.header().lines().len(), 2);
        assert_eq!(r.columns().pair_type, Some(7));
        let mut body = r.into_body();
        let b = body.next_block().unwrap().unwrap();
        assert_eq!(b.first_line, 3);
        assert_eq!(b.n_lines, 2);
        assert!(b.data.ends_with(b"UU\n"));
        assert!(body.next_block().unwrap().is_none());
    }

    #[test]
    fn header_only_and_empty() {
        let r = open(
            b"## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2\n",
        );
        assert_eq!(r.header().lines().len(), 2);
        let mut body = r.into_body();
        assert!(body.next_block().unwrap().is_none());
        let r = open(b"");
        assert!(r.header().is_empty());
        let r = open(
            b"## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2",
        );
        assert_eq!(r.header().lines().len(), 2);
    }

    #[test]
    fn large_header_spanning_reads() {
        let mut data = b"## pairs format v1.0.0\n".to_vec();
        for i in 0..5000 {
            data.extend_from_slice(format!("#chromsize: scaffold_{i} {}\n", 1000 + i).as_bytes());
        }
        data.extend_from_slice(
            b"#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n",
        );
        data.extend_from_slice(b"r1\tscaffold_1\t1\tscaffold_2\t2\t+\t-\tUU\n");
        let r = open(&data);
        assert_eq!(r.header().chromsizes().unwrap().len(), 5000);
        let mut body = r.into_body();
        let b = body.next_block().unwrap().unwrap();
        assert_eq!(b.n_lines, 1);
        assert_eq!(b.first_line, 5003);
    }
}
