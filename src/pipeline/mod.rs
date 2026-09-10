//! The fused `process` pipeline: parse -> flip -> sort -> dedup ->
//! stats / bin / pairs output without text intermediates on disk.

use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use crate::binning::{BinConfig, Binner};
use crate::chroms::{ChromDict, ChromSizes};
use crate::dedup::{DedupConfig, Deduper, Emitted, Outcome};
use crate::error::Result;
use crate::pairs::columns::ColumnMap;
use crate::pairs::header::Header;
use crate::pairs::record::PairKey;
use crate::parse::{AlignmentSource, HicParser, ParseConfig, drive};
use crate::sort::external::{ExternalSorter, SortConfig, SortMetrics};
use crate::sort::key::{ParsedChunk, SortEntry};
use crate::stats::StatsAccumulator;

/// Chunk size (bytes of lines) handed to the sorter at once.
const CHUNK_BYTES: usize = 4 * 1024 * 1024;

/// Accumulates parsed lines into chunks for the sorter.
pub struct ChunkBuilder {
    chunk: ParsedChunk,
}

impl ChunkBuilder {
    /// Empty builder.
    pub fn new() -> Self {
        Self {
            chunk: ParsedChunk {
                data: Vec::with_capacity(CHUNK_BYTES + 4096),
                entries: Vec::with_capacity(CHUNK_BYTES / 80),
            },
        }
    }

    /// Append one record; returns a full chunk when the threshold is hit.
    #[inline]
    pub fn push(&mut self, key: &PairKey, line: &[u8]) -> Option<ParsedChunk> {
        let off = self.chunk.data.len() as u32;
        self.chunk.data.extend_from_slice(line);
        self.chunk.data.push(b'\n');
        self.chunk.entries.push(SortEntry {
            key: *key,
            chunk: 0,
            off,
            len: line.len() as u32,
        });
        if self.chunk.data.len() >= CHUNK_BYTES {
            Some(self.take())
        } else {
            None
        }
    }

    /// Take the current chunk.
    pub fn take(&mut self) -> ParsedChunk {
        std::mem::replace(
            &mut self.chunk,
            ParsedChunk {
                data: Vec::with_capacity(CHUNK_BYTES + 4096),
                entries: Vec::with_capacity(CHUNK_BYTES / 80),
            },
        )
    }

    /// Whether anything is buffered.
    pub fn is_empty(&self) -> bool {
        self.chunk.entries.is_empty()
    }
}

impl Default for ChunkBuilder {
    fn default() -> Self {
        Self::new()
    }
}

/// Stage timings and counters of a fused run.
#[derive(Debug, Default, Clone)]
pub struct ProcessMetrics {
    /// Alignment records read.
    pub alignments_read: u64,
    /// Pairs produced by the parser.
    pub pairs_parsed: u64,
    /// Seconds spent parsing (overlaps with run generation).
    pub parse_seconds: f64,
    /// Sorter metrics.
    pub sort: SortMetrics,
    /// Seconds spent in merge + dedup + outputs.
    pub dedup_seconds: f64,
    /// Duplicates found.
    pub duplicates: u64,
}

/// Run parse + sort over an alignment source, returning the sorted stream.
///
/// The caller supplies the dedup/output stage through [`consume_sorted`].
pub fn parse_and_sort(
    source: AlignmentSource,
    parser: &mut HicParser,
    name: &str,
    sort_cfg: SortConfig,
    stats: Option<&mut StatsAccumulator>,
) -> Result<(crate::sort::external::SortedStream, ProcessMetrics)> {
    let cols = ColumnMap::from_names(parser.header().columns())?;
    let dict = Arc::clone(parser.dict());
    let mut sorter = ExternalSorter::new(sort_cfg, cols, dict)?;
    let mut builder = ChunkBuilder::new();
    let t0 = Instant::now();
    let mut err: Option<crate::error::KiraError> = None;
    drive(source, parser, stats, name, |key, line| {
        if let Some(chunk) = builder.push(key, line)
            && let Err(e) = sorter.push_chunk(chunk)
        {
            err = Some(e);
            return Err(crate::error::KiraError::IoPlain(std::io::Error::other(
                "sorter failed",
            )));
        }
        Ok(())
    })
    .map_err(|e| err.take().unwrap_or(e))?;
    if !builder.is_empty() {
        sorter.push_chunk(builder.take())?;
    }
    let parse_seconds = t0.elapsed().as_secs_f64();
    let (stream, sm) = sorter.finish()?;
    let metrics = ProcessMetrics {
        alignments_read: parser.records_in,
        pairs_parsed: parser.pairs_out,
        parse_seconds,
        sort: sm,
        dedup_seconds: 0.0,
        duplicates: 0,
    };
    Ok((stream, metrics))
}

/// Dedup a sorted stream, routing records through `outputs` and feeding an
/// optional binner with kept pairs.
pub fn consume_sorted<F>(
    stream: &mut crate::sort::external::SortedStream,
    dedup_cfg: DedupConfig,
    dict: &ChromDict,
    mut outputs: F,
    mut binner: Option<&mut Binner>,
) -> Result<crate::dedup::DedupMetrics>
where
    F: FnMut(Emitted<'_>) -> Result<()>,
{
    let mut deduper = Deduper::new(dedup_cfg, dict);
    let mut sink = |e: Emitted<'_>| {
        if let Some(b) = binner.as_deref_mut()
            && e.outcome == Outcome::Unique
        {
            b.observe(e.key, None)?;
        }
        outputs(e)
    };
    while let Some((key, line)) = stream.next_record()? {
        deduper.push(key, line, key.seq + 1, &mut sink)?;
    }
    deduper.finish(&mut sink)?;
    Ok(deduper.metrics().clone())
}

/// Build the standard header for fused output.
pub fn output_header(parser: &HicParser) -> Header {
    parser.header().clone()
}

/// Convenience: build a binner config from process options.
pub fn bin_config(
    resolutions: Vec<u64>,
    chromsizes: &ChromSizes,
    memory_bytes: u64,
    tmpdir: Option<&Path>,
) -> BinConfig {
    BinConfig {
        resolutions,
        chromsizes: chromsizes.clone(),
        min_mapq: None,
        pair_types: None,
        zero_based: false,
        memory_bytes,
        tmpdir: tmpdir.map(Path::to_path_buf),
    }
}

/// Re-export for CLI convenience.
pub type ParseOptions = ParseConfig;
