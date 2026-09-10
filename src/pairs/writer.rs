//! Buffered `.pairs` writer with extension-driven compression.

use std::io::Write;
use std::path::Path;

use crate::error::{KiraError, Result};
use crate::io::compression::{FinishWrite, OutputOptions, open_output};
use crate::pairs::header::Header;

const FLUSH_THRESHOLD: usize = 1 << 20;

/// Writes header and body lines to a (possibly compressed) destination.
pub struct PairsWriter {
    out: Option<Box<dyn FinishWrite>>,
    buf: Vec<u8>,
    bytes_written: u64,
    lines_written: u64,
    name: String,
}

impl PairsWriter {
    /// Create a writer for `path` (`None`/`-` = stdout).
    pub fn create(path: Option<&Path>, opts: OutputOptions) -> Result<Self> {
        let out = open_output(path, opts)?;
        let name = path
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| "-".into());
        Ok(Self::from_boxed(out, name))
    }

    /// Wrap an existing sink.
    pub fn from_boxed(out: Box<dyn FinishWrite>, name: String) -> Self {
        Self {
            out: Some(out),
            buf: Vec::with_capacity(FLUSH_THRESHOLD + 4096),
            bytes_written: 0,
            lines_written: 0,
            name,
        }
    }

    /// Output name for diagnostics.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Write all header lines.
    pub fn write_header(&mut self, header: &Header) -> Result<()> {
        for l in header.lines() {
            self.buf.extend_from_slice(l.as_bytes());
            self.buf.push(b'\n');
        }
        self.maybe_flush()
    }

    /// Write one body line (without trailing newline; one is appended).
    #[inline]
    pub fn write_line(&mut self, line: &[u8]) -> Result<()> {
        self.buf.extend_from_slice(line);
        self.buf.push(b'\n');
        self.lines_written += 1;
        if self.buf.len() >= FLUSH_THRESHOLD {
            self.flush_buf()?;
        }
        Ok(())
    }

    /// Write raw bytes (must already contain newlines as needed).
    pub fn write_raw(&mut self, bytes: &[u8]) -> Result<()> {
        self.buf.extend_from_slice(bytes);
        self.maybe_flush()
    }

    #[inline]
    fn maybe_flush(&mut self) -> Result<()> {
        if self.buf.len() >= FLUSH_THRESHOLD {
            self.flush_buf()?;
        }
        Ok(())
    }

    fn flush_buf(&mut self) -> Result<()> {
        if self.buf.is_empty() {
            return Ok(());
        }
        let out = self
            .out
            .as_mut()
            .ok_or_else(|| KiraError::IoPlain(std::io::Error::other("writer already finished")))?;
        out.write_all(&self.buf)
            .map_err(|e| KiraError::io(Path::new(&self.name), e))?;
        self.bytes_written += self.buf.len() as u64;
        self.buf.clear();
        Ok(())
    }

    /// Bytes handed to the sink so far (uncompressed).
    pub fn bytes_written(&self) -> u64 {
        self.bytes_written + self.buf.len() as u64
    }

    /// Body lines written so far.
    pub fn lines_written(&self) -> u64 {
        self.lines_written
    }

    /// Flush buffered data to the sink (no trailer).
    pub fn flush(&mut self) -> Result<()> {
        self.flush_buf()?;
        if let Some(out) = self.out.as_mut() {
            out.flush()
                .map_err(|e| KiraError::io(Path::new(&self.name), e))?;
        }
        Ok(())
    }

    /// Flush and close; returns the number of uncompressed bytes written.
    pub fn finish(mut self) -> Result<u64> {
        self.flush_buf()?;
        if let Some(out) = self.out.take() {
            out.finish()
                .map_err(|e| KiraError::io(Path::new(&self.name), e))?;
        }
        Ok(self.bytes_written)
    }
}

impl Drop for PairsWriter {
    fn drop(&mut self) {
        if self.out.is_some() {
            let _ = self.flush_buf();
            if let Some(out) = self.out.take() {
                let _ = out.finish();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pairs::reader::PairsReader;
    use std::io::Read;

    #[test]
    fn roundtrip_through_gzip_file() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("out.pairs.gz");
        let header = Header::from_lines([
            "## pairs format v1.0.0",
            "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type",
        ])
        .unwrap();
        let mut w = PairsWriter::create(
            Some(&p),
            OutputOptions {
                threads: 2,
                ..Default::default()
            },
        )
        .unwrap();
        w.write_header(&header).unwrap();
        for i in 0..100_000u32 {
            w.write_line(format!("r{i}\tchr1\t{}\tchr1\t{}\t+\t-\tUU", i + 1, i + 2).as_bytes())
                .unwrap();
        }
        assert_eq!(w.lines_written(), 100_000);
        w.finish().unwrap();
        let r = PairsReader::open(Some(&p), 2).unwrap();
        assert_eq!(r.header(), &header);
        let mut n = 0u64;
        let mut body = r.into_body();
        while let Some(b) = body.next_block().unwrap() {
            n += b.n_lines;
        }
        assert_eq!(n, 100_000);
        let mut raw = Vec::new();
        std::fs::File::open(&p)
            .unwrap()
            .read_to_end(&mut raw)
            .unwrap();
        assert!(crate::io::bgzf::is_bgzf_header(&raw));
    }
}
