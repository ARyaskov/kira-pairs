//! Parallel BGZF (blocked gzip) writer and reader.
//!
//! BGZF is a series of independent gzip members of at most 64 KiB of
//! uncompressed data, which makes both compression and decompression
//! embarrassingly parallel. `.pairs.gz` files written by pairtools are BGZF
//! (bgzip), and BGZF output remains readable by any gzip tool.

use std::collections::BTreeMap;
use std::io::{self, Read, Write};
use std::thread::JoinHandle;

use crossbeam_channel::{Receiver, Sender, bounded};
use flate2::{Compress, Compression, Crc, Decompress, FlushCompress, FlushDecompress};

/// Maximum uncompressed payload per BGZF block (as used by htslib).
pub const BGZF_BLOCK_DATA: usize = 0xff00;
/// Maximum total size of a BGZF block.
pub const BGZF_MAX_BLOCK: usize = 0x10000;
/// The canonical 28-byte BGZF end-of-file marker block.
pub const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];
const HEADER_LEN: usize = 18;
const FOOTER_LEN: usize = 8;

/// True when `head` (at least 18 bytes) looks like a BGZF block header.
pub fn is_bgzf_header(head: &[u8]) -> bool {
    head.len() >= HEADER_LEN
        && head[0] == 0x1f
        && head[1] == 0x8b
        && head[2] == 8
        && head[3] & 4 != 0
        && head[12] == b'B'
        && head[13] == b'C'
        && head[14] == 2
        && head[15] == 0
}

/// Compress one block of data into a complete BGZF member.
pub fn compress_block(data: &[u8], level: u32, scratch: &mut Compress) -> io::Result<Vec<u8>> {
    debug_assert!(data.len() <= BGZF_BLOCK_DATA);
    let mut out = Vec::with_capacity(HEADER_LEN + data.len() + 64);
    out.extend_from_slice(&BGZF_EOF[..HEADER_LEN]);
    scratch.reset();
    let _ = level;
    loop {
        let consumed = scratch.total_in() as usize;
        let before = out.len();
        let status = scratch.compress_vec(&data[consumed..], &mut out, FlushCompress::Finish)?;
        match status {
            flate2::Status::StreamEnd => break,
            flate2::Status::Ok | flate2::Status::BufError => {
                if out.len() == before
                    && scratch.total_in() as usize == consumed
                    && out.capacity() > out.len()
                {
                    return Err(io::Error::other("deflate made no progress on BGZF block"));
                }
                out.reserve(16 * 1024);
            }
        }
    }
    let mut crc = Crc::new();
    crc.update(data);
    out.extend_from_slice(&crc.sum().to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    if out.len() > BGZF_MAX_BLOCK {
        return Err(io::Error::other("BGZF block exceeds 64 KiB"));
    }
    let bsize = (out.len() - 1) as u16;
    out[16..18].copy_from_slice(&bsize.to_le_bytes());
    Ok(out)
}

/// Decompress one complete BGZF member into `out` (cleared first).
pub fn decompress_block(
    block: &[u8],
    out: &mut Vec<u8>,
    scratch: &mut Decompress,
) -> io::Result<()> {
    if block.len() < HEADER_LEN + FOOTER_LEN || !is_bgzf_header(block) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "malformed BGZF block header",
        ));
    }
    let xlen = u16::from_le_bytes([block[10], block[11]]) as usize;
    let data_start = 12 + xlen;
    let n = block.len();
    let isize =
        u32::from_le_bytes([block[n - 4], block[n - 3], block[n - 2], block[n - 1]]) as usize;
    let crc_expected = u32::from_le_bytes([block[n - 8], block[n - 7], block[n - 6], block[n - 5]]);
    if data_start > n - FOOTER_LEN {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "malformed BGZF block",
        ));
    }
    let payload = &block[data_start..n - FOOTER_LEN];
    out.clear();
    out.reserve(isize);
    scratch.reset(false);
    loop {
        let consumed = scratch.total_in() as usize;
        let before = out.len();
        let status = scratch.decompress_vec(&payload[consumed..], out, FlushDecompress::Finish)?;
        match status {
            flate2::Status::StreamEnd => break,
            flate2::Status::Ok | flate2::Status::BufError => {
                if out.len() > isize
                    || (out.len() == before && scratch.total_in() as usize == consumed)
                {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        "truncated or corrupt BGZF block",
                    ));
                }
                out.reserve(isize.max(1024));
            }
        }
    }
    if out.len() != isize {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "BGZF block size mismatch",
        ));
    }
    let mut crc = Crc::new();
    crc.update(out);
    if crc.sum() != crc_expected {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "BGZF block CRC mismatch",
        ));
    }
    Ok(())
}

enum Job {
    Block(u64, Vec<u8>),
    Done,
}

/// Streaming BGZF writer that compresses blocks on a worker pool.
///
/// With `threads <= 1` compression happens inline on the calling thread.
pub struct BgzfWriter<W: Write + Send + 'static> {
    buf: Vec<u8>,
    next_index: u64,
    level: u32,
    // Inline mode.
    inline: Option<(W, Compress)>,
    // Threaded mode.
    job_tx: Option<Sender<Job>>,
    workers: Vec<JoinHandle<()>>,
    writer_thread: Option<JoinHandle<io::Result<W>>>,
    error_rx: Option<Receiver<io::Error>>,
    finished: bool,
}

impl<W: Write + Send + 'static> BgzfWriter<W> {
    /// Create a writer. `threads` is the number of compression workers.
    pub fn new(inner: W, level: u32, threads: usize) -> Self {
        let level = level.min(9);
        let compression = Compression::new(level);
        if threads <= 1 {
            return Self {
                buf: Vec::with_capacity(BGZF_BLOCK_DATA),
                next_index: 0,
                level,
                inline: Some((inner, Compress::new(compression, false))),
                job_tx: None,
                workers: Vec::new(),
                writer_thread: None,
                error_rx: None,
                finished: false,
            };
        }
        let (job_tx, job_rx) = bounded::<Job>(threads * 2);
        let (res_tx, res_rx) = bounded::<(u64, io::Result<Vec<u8>>)>(threads * 2);
        let (err_tx, err_rx) = bounded::<io::Error>(1);
        let mut workers = Vec::with_capacity(threads);
        for _ in 0..threads {
            let rx = job_rx.clone();
            let tx = res_tx.clone();
            workers.push(std::thread::spawn(move || {
                let mut scratch = Compress::new(compression, false);
                while let Ok(job) = rx.recv() {
                    match job {
                        Job::Block(i, data) => {
                            let r = compress_block(&data, level, &mut scratch);
                            if tx.send((i, r)).is_err() {
                                break;
                            }
                        }
                        Job::Done => break,
                    }
                }
            }));
        }
        drop(res_tx);
        let writer_thread = std::thread::spawn(move || {
            let mut inner = inner;
            let mut pending: BTreeMap<u64, Vec<u8>> = BTreeMap::new();
            let mut next = 0u64;
            let mut failed: Option<io::Error> = None;
            while let Ok((i, r)) = res_rx.recv() {
                if failed.is_some() {
                    continue;
                }
                match r {
                    Ok(block) => {
                        pending.insert(i, block);
                        while let Some(b) = pending.remove(&next) {
                            if let Err(e) = inner.write_all(&b) {
                                let _ = err_tx.try_send(io::Error::new(e.kind(), e.to_string()));
                                failed = Some(e);
                                break;
                            }
                            next += 1;
                        }
                    }
                    Err(e) => {
                        let _ = err_tx.try_send(io::Error::new(e.kind(), e.to_string()));
                        failed = Some(e);
                    }
                }
            }
            match failed {
                Some(e) => Err(e),
                None => {
                    inner.write_all(&BGZF_EOF)?;
                    inner.flush()?;
                    Ok(inner)
                }
            }
        });
        Self {
            buf: Vec::with_capacity(BGZF_BLOCK_DATA),
            next_index: 0,
            level,
            inline: None,
            job_tx: Some(job_tx),
            workers,
            writer_thread: Some(writer_thread),
            error_rx: Some(err_rx),
            finished: false,
        }
    }

    fn check_error(&self) -> io::Result<()> {
        if let Some(rx) = &self.error_rx
            && let Ok(e) = rx.try_recv()
        {
            return Err(e);
        }
        Ok(())
    }

    fn emit_block(&mut self) -> io::Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        let data = std::mem::replace(&mut self.buf, Vec::with_capacity(BGZF_BLOCK_DATA));
        let idx = self.next_index;
        self.next_index += 1;
        if let Some((inner, scratch)) = self.inline.as_mut() {
            let block = compress_block(&data, self.level, scratch)?;
            inner.write_all(&block)?;
            return Ok(());
        }
        self.check_error()?;
        if let Some(tx) = &self.job_tx
            && tx.send(Job::Block(idx, data)).is_err()
        {
            self.check_error()?;
            return Err(io::Error::other("BGZF compression workers stopped"));
        }
        Ok(())
    }

    /// Flush remaining data, write the EOF block and return the inner writer.
    pub fn finish(mut self) -> io::Result<W> {
        self.finish_inner()
    }

    fn finish_inner(&mut self) -> io::Result<W> {
        if self.finished {
            return Err(io::Error::other("BGZF writer already finished"));
        }
        self.finished = true;
        self.emit_block()?;
        if let Some((mut inner, _)) = self.inline.take() {
            inner.write_all(&BGZF_EOF)?;
            inner.flush()?;
            return Ok(inner);
        }
        if let Some(tx) = self.job_tx.take() {
            for _ in 0..self.workers.len() {
                let _ = tx.send(Job::Done);
            }
            drop(tx);
        }
        for w in self.workers.drain(..) {
            let _ = w.join();
        }
        match self.writer_thread.take() {
            Some(h) => h
                .join()
                .map_err(|_| io::Error::other("BGZF writer thread panicked"))?,
            None => Err(io::Error::other("BGZF writer thread missing")),
        }
    }
}

impl<W: Write + Send + 'static> Write for BgzfWriter<W> {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        let mut rest = data;
        while !rest.is_empty() {
            let room = BGZF_BLOCK_DATA - self.buf.len();
            let take = room.min(rest.len());
            self.buf.extend_from_slice(&rest[..take]);
            rest = &rest[take..];
            if self.buf.len() >= BGZF_BLOCK_DATA {
                self.emit_block()?;
            }
        }
        Ok(data.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        // Data is flushed in whole blocks at finish; a mid-stream flush only
        // pushes the current partial block out.
        self.emit_block()
    }
}

impl<W: Write + Send + 'static> Drop for BgzfWriter<W> {
    fn drop(&mut self) {
        if !self.finished {
            // Best effort: callers should use `finish` to observe errors.
            let _ = self.finish_inner();
        }
    }
}

/// Parallel BGZF reader: a reader thread splits the stream into blocks,
/// workers inflate them, and blocks are handed out in order.
pub struct BgzfReader {
    rx: Receiver<(u64, io::Result<Vec<u8>>)>,
    pending: BTreeMap<u64, Vec<u8>>,
    next: u64,
    current: Vec<u8>,
    pos: usize,
    done: bool,
    threads: Vec<JoinHandle<()>>,
}

impl BgzfReader {
    /// Start reading BGZF data from `inner` using `threads` inflate workers.
    pub fn new<R: Read + Send + 'static>(inner: R, threads: usize) -> Self {
        let threads_n = threads.max(1);
        let (raw_tx, raw_rx) = bounded::<(u64, io::Result<Vec<u8>>)>(threads_n * 4);
        let (out_tx, out_rx) = bounded::<(u64, io::Result<Vec<u8>>)>(threads_n * 4);
        let mut handles = Vec::new();
        handles.push(std::thread::spawn(move || {
            let mut r = inner;
            let mut idx = 0u64;
            loop {
                match read_raw_block(&mut r) {
                    Ok(Some(block)) => {
                        if raw_tx.send((idx, Ok(block))).is_err() {
                            break;
                        }
                        idx += 1;
                    }
                    Ok(None) => break,
                    Err(e) => {
                        let _ = raw_tx.send((idx, Err(e)));
                        break;
                    }
                }
            }
        }));
        for _ in 0..threads_n {
            let rx = raw_rx.clone();
            let tx = out_tx.clone();
            handles.push(std::thread::spawn(move || {
                let mut scratch = Decompress::new(false);
                let mut out = Vec::new();
                while let Ok((i, r)) = rx.recv() {
                    let res = match r {
                        Ok(block) => decompress_block(&block, &mut out, &mut scratch)
                            .map(|()| std::mem::take(&mut out)),
                        Err(e) => Err(e),
                    };
                    if tx.send((i, res)).is_err() {
                        break;
                    }
                }
            }));
        }
        Self {
            rx: out_rx,
            pending: BTreeMap::new(),
            next: 0,
            current: Vec::new(),
            pos: 0,
            done: false,
            threads: handles,
        }
    }

    fn fill(&mut self) -> io::Result<bool> {
        loop {
            if let Some(b) = self.pending.remove(&self.next) {
                self.next += 1;
                self.current = b;
                self.pos = 0;
                if self.current.is_empty() {
                    continue;
                }
                return Ok(true);
            }
            if self.done {
                return Ok(false);
            }
            match self.rx.recv() {
                Ok((i, Ok(b))) => {
                    self.pending.insert(i, b);
                }
                Ok((_, Err(e))) => {
                    self.done = true;
                    return Err(e);
                }
                Err(_) => {
                    self.done = true;
                }
            }
        }
    }
}

impl Read for BgzfReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.pos >= self.current.len() && !self.fill()? {
            return Ok(0);
        }
        let n = (self.current.len() - self.pos).min(buf.len());
        buf[..n].copy_from_slice(&self.current[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

impl Drop for BgzfReader {
    fn drop(&mut self) {
        // Disconnect so that producer threads exit, then join them.
        let (_tx, rx) = bounded(0);
        self.rx = rx;
        for h in self.threads.drain(..) {
            let _ = h.join();
        }
    }
}

fn read_exact_or_eof<R: Read>(r: &mut R, buf: &mut [u8]) -> io::Result<usize> {
    let mut n = 0;
    while n < buf.len() {
        match r.read(&mut buf[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(n)
}

/// Read one raw BGZF block (header + payload + footer). Returns `None` at a
/// clean end of stream; a partial block is a truncation error.
fn read_raw_block<R: Read>(r: &mut R) -> io::Result<Option<Vec<u8>>> {
    let mut head = [0u8; HEADER_LEN];
    let n = read_exact_or_eof(r, &mut head)?;
    if n == 0 {
        return Ok(None);
    }
    if n < HEADER_LEN || !is_bgzf_header(&head) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "truncated or non-BGZF gzip block (input is not valid BGZF)",
        ));
    }
    let bsize = u16::from_le_bytes([head[16], head[17]]) as usize + 1;
    if bsize < HEADER_LEN + FOOTER_LEN {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "invalid BGZF block size",
        ));
    }
    let mut block = vec![0u8; bsize];
    block[..HEADER_LEN].copy_from_slice(&head);
    let got = read_exact_or_eof(r, &mut block[HEADER_LEN..])?;
    if got != bsize - HEADER_LEN {
        return Err(io::Error::new(
            io::ErrorKind::UnexpectedEof,
            "truncated BGZF block (file is incomplete)",
        ));
    }
    Ok(Some(block))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(threads_w: usize, threads_r: usize, size: usize) {
        let data: Vec<u8> = (0..size).map(|i| (i % 251) as u8).collect();
        let mut w = BgzfWriter::new(Vec::new(), 6, threads_w);
        w.write_all(&data).unwrap();
        let compressed = w.finish().unwrap();
        assert!(compressed.ends_with(&BGZF_EOF));
        assert!(is_bgzf_header(&compressed));
        // Any gzip reader must accept it.
        let mut plain = Vec::new();
        flate2::read::MultiGzDecoder::new(&compressed[..])
            .read_to_end(&mut plain)
            .unwrap();
        assert_eq!(plain, data);
        let mut r = BgzfReader::new(std::io::Cursor::new(compressed), threads_r);
        let mut back = Vec::new();
        r.read_to_end(&mut back).unwrap();
        assert_eq!(back, data);
    }

    #[test]
    fn roundtrips_inline_and_threaded() {
        roundtrip(1, 1, 0);
        roundtrip(1, 1, 10);
        roundtrip(1, 2, 300_000);
        roundtrip(4, 4, 1_000_003);
    }

    #[test]
    fn detects_truncation() {
        let mut w = BgzfWriter::new(Vec::new(), 6, 1);
        w.write_all(&vec![7u8; 100_000]).unwrap();
        let compressed = w.finish().unwrap();
        let cut = &compressed[..compressed.len() - 40];
        let mut r = BgzfReader::new(std::io::Cursor::new(cut.to_vec()), 2);
        let mut back = Vec::new();
        assert!(r.read_to_end(&mut back).is_err());
    }
}
