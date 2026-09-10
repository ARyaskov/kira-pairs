//! Compression auto-detection for inputs and extension-driven compression
//! for outputs.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::thread::JoinHandle;

use crossbeam_channel::{Receiver, bounded};

use crate::error::{KiraError, Result};
use crate::io::bgzf::{BgzfReader, BgzfWriter, is_bgzf_header};
use crate::io::buffered::DEFAULT_BLOCK_SIZE;

/// Compression formats understood by kira-pairs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Compression {
    /// Plain text.
    None,
    /// Single- or multi-member gzip that is not BGZF.
    Gzip,
    /// Blocked gzip (bgzip output).
    Bgzf,
    /// LZ4 frame format (lz4c output).
    Lz4,
}

impl Compression {
    /// Detect a format from the first bytes of a stream.
    pub fn detect(head: &[u8]) -> Self {
        if head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b {
            if is_bgzf_header(head) {
                Self::Bgzf
            } else {
                Self::Gzip
            }
        } else if head.len() >= 4 && head[..4] == [0x04, 0x22, 0x4d, 0x18] {
            Self::Lz4
        } else {
            Self::None
        }
    }
}

/// Output compression chosen from a file extension.
pub fn output_compression_for_path(path: &Path) -> Compression {
    let name = path.to_string_lossy().to_ascii_lowercase();
    if name.ends_with(".gz") || name.ends_with(".bgz") || name.ends_with(".bgzf") {
        Compression::Bgzf
    } else if name.ends_with(".lz4") {
        Compression::Lz4
    } else {
        Compression::None
    }
}

/// A source of decompressed bytes with its detected format.
pub struct InputSource {
    /// Human-readable name for diagnostics (`-` for stdin).
    pub name: String,
    /// Detected compression.
    pub compression: Compression,
    /// Decompressed byte stream.
    pub reader: Box<dyn Read + Send>,
}

/// A reader that pumps another reader on a background thread so that
/// decompression overlaps with parsing.
pub struct ThreadedReader {
    rx: Receiver<io::Result<Vec<u8>>>,
    current: Vec<u8>,
    pos: usize,
    handle: Option<JoinHandle<()>>,
}

impl ThreadedReader {
    /// Spawn the pump thread.
    pub fn new<R: Read + Send + 'static>(mut inner: R, block: usize, depth: usize) -> Self {
        let (tx, rx) = bounded(depth.max(1));
        let handle = std::thread::spawn(move || {
            loop {
                let mut buf = vec![0u8; block];
                let mut n = 0;
                let mut err = None;
                while n < block {
                    match inner.read(&mut buf[n..]) {
                        Ok(0) => break,
                        Ok(k) => n += k,
                        Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                        Err(e) => {
                            err = Some(e);
                            break;
                        }
                    }
                }
                buf.truncate(n);
                if n > 0 && tx.send(Ok(buf)).is_err() {
                    return;
                }
                if let Some(e) = err {
                    let _ = tx.send(Err(e));
                    return;
                }
                if n == 0 {
                    return;
                }
            }
        });
        Self {
            rx,
            current: Vec::new(),
            pos: 0,
            handle: Some(handle),
        }
    }
}

impl Read for ThreadedReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.pos >= self.current.len() {
            match self.rx.recv() {
                Ok(Ok(b)) => {
                    self.current = b;
                    self.pos = 0;
                }
                Ok(Err(e)) => return Err(e),
                Err(_) => return Ok(0),
            }
        }
        let n = (self.current.len() - self.pos).min(buf.len());
        buf[..n].copy_from_slice(&self.current[self.pos..self.pos + n]);
        self.pos += n;
        Ok(n)
    }
}

impl Drop for ThreadedReader {
    fn drop(&mut self) {
        let (_tx, rx) = bounded(0);
        self.rx = rx;
        if let Some(h) = self.handle.take() {
            let _ = h.join();
        }
    }
}

fn open_raw(path: Option<&Path>) -> Result<(String, Box<dyn Read + Send>)> {
    match path {
        None => Ok(("-".to_string(), Box::new(io::stdin()))),
        Some(p) if p.as_os_str() == "-" => Ok(("-".to_string(), Box::new(io::stdin()))),
        Some(p) => {
            let f = File::open(p).map_err(|e| KiraError::io(p, e))?;
            Ok((p.display().to_string(), Box::new(f)))
        }
    }
}

/// Open an input path (`None` or `-` = stdin), detect compression by magic
/// bytes and return a decompressed stream. `threads` bounds the number of
/// parallel BGZF inflate workers.
pub fn open_input(path: Option<&Path>, threads: usize) -> Result<InputSource> {
    let (name, mut raw) = open_raw(path)?;
    let mut head = vec![0u8; 18];
    let mut n = 0;
    while n < head.len() {
        match raw.read(&mut head[n..]) {
            Ok(0) => break,
            Ok(k) => n += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(KiraError::io(PathBuf::from(&name), e)),
        }
    }
    head.truncate(n);
    let compression = Compression::detect(&head);
    let chained: Box<dyn Read + Send> = Box::new(io::Cursor::new(head).chain(raw));
    let reader: Box<dyn Read + Send> = match compression {
        Compression::None => Box::new(BufReader::with_capacity(DEFAULT_BLOCK_SIZE, chained)),
        Compression::Bgzf => Box::new(BgzfReader::new(chained, threads.max(1))),
        Compression::Gzip => {
            let dec = flate2::read::MultiGzDecoder::new(BufReader::with_capacity(1 << 20, chained));
            Box::new(ThreadedReader::new(dec, DEFAULT_BLOCK_SIZE, 4))
        }
        Compression::Lz4 => {
            let dec =
                lz4_flex::frame::FrameDecoder::new(BufReader::with_capacity(1 << 20, chained));
            Box::new(ThreadedReader::new(dec, DEFAULT_BLOCK_SIZE, 4))
        }
    };
    log::debug!("input {name}: {compression:?}");
    Ok(InputSource {
        name,
        compression,
        reader,
    })
}

/// A writer that must be explicitly finished to flush trailers.
pub trait FinishWrite: Write + Send {
    /// Flush everything, write any trailer, and close.
    fn finish(self: Box<Self>) -> io::Result<()>;
}

struct PlainOut<W: Write + Send>(BufWriter<W>);

impl<W: Write + Send> Write for PlainOut<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.write(buf)
    }
    fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        self.0.write_all(buf)
    }
    fn flush(&mut self) -> io::Result<()> {
        self.0.flush()
    }
}

impl<W: Write + Send> FinishWrite for PlainOut<W> {
    fn finish(mut self: Box<Self>) -> io::Result<()> {
        self.0.flush()
    }
}

struct BgzfOut<W: Write + Send + 'static>(Option<BgzfWriter<BufWriter<W>>>);

impl<W: Write + Send + 'static> Write for BgzfOut<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self.0.as_mut() {
            Some(w) => w.write(buf),
            None => Err(io::Error::other("writer finished")),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self.0.as_mut() {
            Some(w) => w.flush(),
            None => Ok(()),
        }
    }
}

impl<W: Write + Send + 'static> FinishWrite for BgzfOut<W> {
    fn finish(mut self: Box<Self>) -> io::Result<()> {
        if let Some(w) = self.0.take() {
            let mut inner = w.finish()?;
            inner.flush()?;
        }
        Ok(())
    }
}

struct Lz4Out<W: Write + Send>(Option<lz4_flex::frame::FrameEncoder<BufWriter<W>>>);

impl<W: Write + Send> Write for Lz4Out<W> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self.0.as_mut() {
            Some(w) => w.write(buf),
            None => Err(io::Error::other("writer finished")),
        }
    }
    fn flush(&mut self) -> io::Result<()> {
        match self.0.as_mut() {
            Some(w) => w.flush(),
            None => Ok(()),
        }
    }
}

impl<W: Write + Send> FinishWrite for Lz4Out<W> {
    fn finish(mut self: Box<Self>) -> io::Result<()> {
        if let Some(w) = self.0.take() {
            let mut inner = w.finish().map_err(io::Error::other)?;
            inner.flush()?;
        }
        Ok(())
    }
}

/// Options controlling output compression.
#[derive(Debug, Clone, Copy)]
pub struct OutputOptions {
    /// Format; `None` means "choose from the extension".
    pub compression: Option<Compression>,
    /// Compression level (gzip/BGZF 1-9).
    pub level: u32,
    /// Compression worker threads.
    pub threads: usize,
}

impl Default for OutputOptions {
    fn default() -> Self {
        Self {
            compression: None,
            level: 6,
            threads: 1,
        }
    }
}

/// Open an output path (`None` or `-` = stdout) with the compression
/// implied by its extension.
pub fn open_output(path: Option<&Path>, opts: OutputOptions) -> Result<Box<dyn FinishWrite>> {
    let is_stdout = path.is_none_or(|p| p.as_os_str() == "-");
    let compression = opts.compression.unwrap_or_else(|| {
        path.map(output_compression_for_path)
            .unwrap_or(Compression::None)
    });
    let raw: Box<dyn Write + Send> = if is_stdout {
        Box::new(io::stdout())
    } else {
        let p = path.unwrap_or_else(|| Path::new("-"));
        Box::new(File::create(p).map_err(|e| KiraError::io(p, e))?)
    };
    let buffered = BufWriter::with_capacity(1 << 20, raw);
    let w: Box<dyn FinishWrite> = match compression {
        Compression::None => Box::new(PlainOut(buffered)),
        Compression::Bgzf | Compression::Gzip => Box::new(BgzfOut(Some(BgzfWriter::new(
            buffered,
            opts.level,
            opts.threads,
        )))),
        Compression::Lz4 => Box::new(Lz4Out(Some(lz4_flex::frame::FrameEncoder::new(buffered)))),
    };
    Ok(w)
}

/// Open several output paths, sharing a single writer for paths that refer
/// to the same file (pairtools semantics for `--output-dups -`, etc.).
pub fn same_output(a: Option<&Path>, b: Option<&Path>) -> bool {
    match (a, b) {
        (Some(a), Some(b)) => {
            let sa = a.as_os_str() == "-";
            let sb = b.as_os_str() == "-";
            if sa || sb {
                return sa && sb;
            }
            let ca = std::fs::canonicalize(a).unwrap_or_else(|_| a.to_path_buf());
            let cb = std::fs::canonicalize(b).unwrap_or_else(|_| b.to_path_buf());
            ca == cb || (a == b)
        }
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_formats() {
        assert_eq!(Compression::detect(b"## pairs"), Compression::None);
        assert_eq!(
            Compression::detect(&[0x1f, 0x8b, 8, 0, 0, 0]),
            Compression::Gzip
        );
        assert_eq!(
            Compression::detect(&crate::io::bgzf::BGZF_EOF),
            Compression::Bgzf
        );
        assert_eq!(
            Compression::detect(&[0x04, 0x22, 0x4d, 0x18, 0]),
            Compression::Lz4
        );
        assert_eq!(
            output_compression_for_path(Path::new("a.pairs.gz")),
            Compression::Bgzf
        );
        assert_eq!(
            output_compression_for_path(Path::new("a.pairs.LZ4")),
            Compression::Lz4
        );
        assert_eq!(
            output_compression_for_path(Path::new("a.pairs")),
            Compression::None
        );
    }

    fn roundtrip_file(ext: &str, threads: usize) {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join(format!("t.pairs{ext}"));
        let payload: Vec<u8> = (0..500_000u32)
            .flat_map(|i| format!("line {i}\n").into_bytes())
            .collect();
        let mut w = open_output(
            Some(&p),
            OutputOptions {
                compression: None,
                level: 4,
                threads,
            },
        )
        .unwrap();
        w.write_all(&payload).unwrap();
        w.finish().unwrap();
        let mut src = open_input(Some(&p), threads).unwrap();
        let mut back = Vec::new();
        src.reader.read_to_end(&mut back).unwrap();
        assert_eq!(back, payload);
        let expected = output_compression_for_path(&p);
        assert_eq!(src.compression, expected);
    }

    #[test]
    fn roundtrips_all_formats() {
        roundtrip_file("", 1);
        roundtrip_file(".gz", 1);
        roundtrip_file(".gz", 4);
        roundtrip_file(".lz4", 2);
    }

    #[test]
    fn reads_plain_gzip() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.gz");
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(b"hello\nworld\n").unwrap();
        std::fs::write(&p, enc.finish().unwrap()).unwrap();
        let mut src = open_input(Some(&p), 2).unwrap();
        assert_eq!(src.compression, Compression::Gzip);
        let mut back = String::new();
        src.reader.read_to_string(&mut back).unwrap();
        assert_eq!(back, "hello\nworld\n");
    }

    #[test]
    fn truncated_gzip_is_an_error() {
        let mut enc = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        enc.write_all(&vec![1u8; 100_000]).unwrap();
        let data = enc.finish().unwrap();
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.gz");
        std::fs::write(&p, &data[..data.len() / 2]).unwrap();
        let mut src = open_input(Some(&p), 1).unwrap();
        let mut back = Vec::new();
        assert!(src.reader.read_to_end(&mut back).is_err());
    }
}
