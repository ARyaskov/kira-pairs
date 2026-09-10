//! Large block reads that are split on newline boundaries so that blocks can
//! be parsed independently and in parallel.

use std::io::{self, Read};

use memchr::memrchr;

/// Default block size used by streaming readers and channels.
pub const DEFAULT_BLOCK_SIZE: usize = 4 * 1024 * 1024;

/// A block of complete lines from the input stream.
#[derive(Debug, Default, Clone)]
pub struct LineBlock {
    /// Raw bytes. Every line ends with `\n` (the reader appends one to an
    /// unterminated final line).
    pub data: Vec<u8>,
    /// 1-based line number of the first line in this block.
    pub first_line: u64,
    /// Number of lines in this block.
    pub n_lines: u64,
    /// 0-based ordinal of this block in the stream.
    pub index: u64,
}

/// Reads an input stream in large blocks and yields [`LineBlock`]s whose
/// boundaries always fall on newlines.
pub struct BlockReader<R: Read> {
    inner: R,
    block_size: usize,
    carry: Vec<u8>,
    next_line: u64,
    next_index: u64,
    eof: bool,
    bytes_read: u64,
}

impl<R: Read> BlockReader<R> {
    /// Wrap a reader with the given target block size.
    pub fn new(inner: R, block_size: usize) -> Self {
        Self {
            inner,
            block_size: block_size.max(64 * 1024),
            carry: Vec::new(),
            next_line: 1,
            next_index: 0,
            eof: false,
            bytes_read: 0,
        }
    }

    /// Wrap a reader and seed it with bytes that were already consumed from
    /// the underlying stream (e.g. after header parsing). `first_line` is the
    /// line number of the first line in `prefix`.
    pub fn with_prefix(inner: R, block_size: usize, prefix: Vec<u8>, first_line: u64) -> Self {
        let mut r = Self::new(inner, block_size);
        r.carry = prefix;
        r.next_line = first_line;
        r
    }

    /// Bytes read from the underlying reader so far.
    pub fn bytes_read(&self) -> u64 {
        self.bytes_read
    }

    /// Read the next block, or `None` at end of stream.
    pub fn next_block(&mut self) -> io::Result<Option<LineBlock>> {
        if self.eof && self.carry.is_empty() {
            return Ok(None);
        }
        let mut buf = std::mem::take(&mut self.carry);
        let target = self.block_size;
        if buf.capacity() < target + 1 {
            buf.reserve(target + 1 - buf.len());
        }
        // Fill until we have at least `target` bytes or hit EOF.
        while !self.eof && buf.len() < target {
            let old = buf.len();
            buf.resize(target.max(old + 1), 0);
            match self.inner.read(&mut buf[old..]) {
                Ok(0) => {
                    buf.truncate(old);
                    self.eof = true;
                }
                Ok(n) => {
                    buf.truncate(old + n);
                    self.bytes_read += n as u64;
                }
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {
                    buf.truncate(old);
                }
                Err(e) => return Err(e),
            }
        }
        if buf.is_empty() {
            return Ok(None);
        }
        // Split at the last newline; the remainder is carried over. At end
        // of input everything is emitted (an unterminated final line gets a
        // newline appended).
        if self.eof {
            if !buf.ends_with(b"\n") {
                buf.push(b'\n');
            }
        } else {
            match memrchr(b'\n', &buf) {
                Some(p) => {
                    let cut = p + 1;
                    if cut < buf.len() {
                        self.carry.extend_from_slice(&buf[cut..]);
                        buf.truncate(cut);
                    }
                }
                None => {
                    // A single line longer than the block: keep growing.
                    self.carry = buf;
                    self.block_size *= 2;
                    return self.next_block();
                }
            }
        }
        let n_lines = memchr::memchr_iter(b'\n', &buf).count() as u64;
        let block = LineBlock {
            data: buf,
            first_line: self.next_line,
            n_lines,
            index: self.next_index,
        };
        self.next_line += n_lines;
        self.next_index += 1;
        Ok(Some(block))
    }
}

impl<R: Read> Iterator for BlockReader<R> {
    type Item = io::Result<LineBlock>;

    fn next(&mut self) -> Option<Self::Item> {
        self.next_block().transpose()
    }
}

/// Iterate over the lines of a block (without the trailing `\n`).
pub fn lines(block: &[u8]) -> impl Iterator<Item = &[u8]> + '_ {
    let mut start = 0usize;
    memchr::memchr_iter(b'\n', block).map(move |p| {
        let line = &block[start..p];
        start = p + 1;
        line.strip_suffix(b"\r").unwrap_or(line)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn splits_on_newlines() {
        let data = b"line1\nline2\nline3\nlast".to_vec();
        let mut r = BlockReader::new(&data[..], 64 * 1024);
        let mut blocks = Vec::new();
        while let Some(b) = r.next_block().unwrap() {
            blocks.push(b);
        }
        let joined: Vec<u8> = blocks.iter().flat_map(|b| b.data.clone()).collect();
        assert_eq!(joined, b"line1\nline2\nline3\nlast\n");
        let total: u64 = blocks.iter().map(|b| b.n_lines).sum();
        assert_eq!(total, 4);
        assert_eq!(blocks[0].first_line, 1);
    }

    #[test]
    fn long_lines_grow_block() {
        let long = vec![b'x'; 200_000];
        let mut data = long.clone();
        data.push(b'\n');
        data.extend_from_slice(b"short\n");
        let mut r = BlockReader::new(&data[..], 64 * 1024);
        let b = r.next_block().unwrap().unwrap();
        assert!(b.data.starts_with(&long));
        assert_eq!(b.n_lines, 2);
        assert!(r.next_block().unwrap().is_none());
    }

    #[test]
    fn line_iterator_strips_cr() {
        let v: Vec<&[u8]> = lines(b"a\r\nb\n").collect();
        assert_eq!(v, vec![&b"a"[..], &b"b"[..]]);
    }
}
