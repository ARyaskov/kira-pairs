//! The external sorter: parallel parsing, bounded run generation, parallel
//! run serialisation and (multi-pass) k-way merging.

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;
use std::time::Instant;

use crossbeam_channel::{Sender, bounded};
use rayon::prelude::*;

use crate::chroms::ChromDict;
use crate::error::{KiraError, Location, Result};
use crate::io::buffered::LineBlock;
use crate::io::temp::TempManager;
use crate::memory::MemoryBudget;
use crate::pairs::columns::ColumnMap;
use crate::pairs::record::{PairKey, parse_line};
use crate::sort::key::{ParsedChunk, SortEntry, SortKeyContext};
use crate::sort::merge::{MemorySource, Merger, RunSource, Source};
use crate::sort::run::{BlockEncoder, DEFAULT_RUN_BLOCK, RunReader, RunWriter, compress_payload};

/// Sorter configuration.
#[derive(Debug, Clone)]
pub struct SortConfig {
    /// Worker threads for parsing, sorting and run compression.
    pub threads: usize,
    /// Memory budget.
    pub budget: MemoryBudget,
    /// Base directory for temporary runs (`None` = system default).
    pub tmpdir: Option<PathBuf>,
    /// LZ4-compress run blocks.
    pub compress_runs: bool,
    /// Maximum number of runs merged at once.
    pub max_fan_in: usize,
    /// Extra sort columns compared after `pair_type` (0-based indices).
    pub extra_cols: Vec<usize>,
    /// Raw block size used in run files.
    pub run_block_size: usize,
    /// Input name for error messages.
    pub input_name: String,
}

impl SortConfig {
    /// Configuration with sensible defaults for the given resources.
    pub fn new(threads: usize, budget: MemoryBudget) -> Self {
        Self {
            threads: threads.max(1),
            budget,
            tmpdir: None,
            compress_runs: true,
            max_fan_in: 64,
            extra_cols: Vec::new(),
            run_block_size: DEFAULT_RUN_BLOCK,
            input_name: "-".into(),
        }
    }
}

/// Counters collected by the sorter.
#[derive(Debug, Default, Clone)]
pub struct SortMetrics {
    /// Records sorted.
    pub records: u64,
    /// Runs written to disk.
    pub runs: u64,
    /// Bytes written to temporary storage.
    pub temp_bytes: u64,
    /// Seconds spent sorting and writing runs (run-writer thread).
    pub run_seconds: f64,
    /// Seconds spent in intermediate merge passes.
    pub merge_pass_seconds: f64,
    /// Number of intermediate merge passes.
    pub merge_passes: u64,
    /// Largest number of records buffered in memory at once.
    pub peak_records_buffered: u64,
}

struct Shared {
    error: Mutex<Option<KiraError>>,
    records: AtomicU64,
    peak_buffered: AtomicU64,
}

impl Shared {
    fn set_error(&self, e: KiraError) {
        let mut g = self.error.lock().unwrap_or_else(|p| p.into_inner());
        if g.is_none() {
            *g = Some(e);
        }
    }

    fn take_error(&self) -> Option<KiraError> {
        self.error.lock().unwrap_or_else(|p| p.into_inner()).take()
    }
}

struct RunOutput {
    paths: Vec<PathBuf>,
    temp: Option<Arc<TempManager>>,
    run_seconds: f64,
    temp_bytes: u64,
}

struct Collected {
    remaining: Vec<ParsedChunk>,
}

/// Parallel external merge sorter for `.pairs` records.
pub struct ExternalSorter {
    cfg: SortConfig,
    dict: Arc<ChromDict>,
    cols: Arc<ColumnMap>,
    shared: Arc<Shared>,
    block_tx: Option<Sender<LineBlock>>,
    chunk_tx: Option<Sender<ParsedChunk>>,
    parsers: Vec<JoinHandle<()>>,
    collector: Option<JoinHandle<Collected>>,
    run_writer: Option<JoinHandle<Result<RunOutput>>>,
    run_capacity: usize,
    pushed_records: u64,
}

impl ExternalSorter {
    /// Start a sorter. Parser and run-writer threads are spawned eagerly.
    pub fn new(cfg: SortConfig, cols: ColumnMap, dict: Arc<ChromDict>) -> Result<Self> {
        let threads = cfg.threads.max(1);
        let cols = Arc::new(cols);
        let shared = Arc::new(Shared {
            error: Mutex::new(None),
            records: AtomicU64::new(0),
            peak_buffered: AtomicU64::new(0),
        });
        // Two runs may be resident: the one being filled and the one being
        // sorted and written.
        let run_capacity = (cfg.budget.records / 2).max(32 * 1024 * 1024) as usize;

        let (block_tx, block_rx) = bounded::<LineBlock>(threads * 2);
        let (chunk_tx, chunk_rx) = bounded::<ParsedChunk>(threads * 2);
        let (run_tx, run_rx) = bounded::<Vec<ParsedChunk>>(0);

        let mut parsers = Vec::with_capacity(threads);
        for _ in 0..threads {
            let rx = block_rx.clone();
            let tx = chunk_tx.clone();
            let cols = Arc::clone(&cols);
            let dict = Arc::clone(&dict);
            let shared = Arc::clone(&shared);
            let name = cfg.input_name.clone();
            parsers.push(std::thread::spawn(move || {
                let mut ends: Vec<u32> = Vec::with_capacity(32);
                while let Ok(block) = rx.recv() {
                    match parse_block(block, &cols, &dict, &name, &mut ends) {
                        Ok(chunk) => {
                            if tx.send(chunk).is_err() {
                                break;
                            }
                        }
                        Err(e) => {
                            shared.set_error(e);
                            break;
                        }
                    }
                }
            }));
        }
        drop(block_rx);

        let shared_c = Arc::clone(&shared);
        let collector = std::thread::spawn(move || {
            let mut current: Vec<ParsedChunk> = Vec::new();
            let mut bytes = 0usize;
            let mut buffered = 0u64;
            let mut failed = false;
            while let Ok(chunk) = chunk_rx.recv() {
                if failed {
                    continue;
                }
                bytes += chunk.memory_size();
                buffered += chunk.entries.len() as u64;
                shared_c
                    .records
                    .fetch_add(chunk.entries.len() as u64, Ordering::Relaxed);
                current.push(chunk);
                if bytes >= run_capacity {
                    shared_c
                        .peak_buffered
                        .fetch_max(buffered, Ordering::Relaxed);
                    let run = std::mem::take(&mut current);
                    bytes = 0;
                    buffered = 0;
                    if run_tx.send(run).is_err() {
                        failed = true;
                    }
                }
            }
            shared_c
                .peak_buffered
                .fetch_max(buffered, Ordering::Relaxed);
            drop(run_tx);
            Collected { remaining: current }
        });

        let cfg_w = cfg.clone();
        let dict_w = Arc::clone(&dict);
        let cols_w = Arc::clone(&cols);
        let shared_w = Arc::clone(&shared);
        let run_writer = std::thread::spawn(move || -> Result<RunOutput> {
            let mut out = RunOutput {
                paths: Vec::new(),
                temp: None,
                run_seconds: 0.0,
                temp_bytes: 0,
            };
            while let Ok(chunks) = run_rx.recv() {
                let t0 = Instant::now();
                let temp = match &out.temp {
                    Some(t) => Arc::clone(t),
                    None => {
                        let t = Arc::new(TempManager::new(
                            cfg_w.tmpdir.as_deref(),
                            "kira-pairs-sort",
                        )?);
                        out.temp = Some(Arc::clone(&t));
                        t
                    }
                };
                let path = temp.next_path("run");
                let ctx = make_ctx(&dict_w, &cfg_w, &cols_w);
                let entries = sort_chunks(&chunks, &ctx, cfg_w.threads);
                match write_run(&path, &chunks, &entries, &cfg_w) {
                    Ok(bytes) => {
                        temp.add_bytes(bytes);
                        out.temp_bytes += bytes;
                        out.paths.push(path);
                    }
                    Err(e) => {
                        shared_w.set_error(e);
                        break;
                    }
                }
                out.run_seconds += t0.elapsed().as_secs_f64();
                log::info!("wrote run {} ({} records)", out.paths.len(), entries.len());
            }
            Ok(out)
        });

        Ok(Self {
            cfg,
            dict,
            cols,
            shared,
            block_tx: Some(block_tx),
            chunk_tx: Some(chunk_tx),
            parsers,
            collector: Some(collector),
            run_writer: Some(run_writer),
            run_capacity,
            pushed_records: 0,
        })
    }

    /// Column map used for parsing.
    pub fn columns(&self) -> &ColumnMap {
        &self.cols
    }

    /// Chromosome dictionary shared with callers.
    pub fn dict(&self) -> &Arc<ChromDict> {
        &self.dict
    }

    /// Bytes of records per run.
    pub fn run_capacity(&self) -> usize {
        self.run_capacity
    }

    fn check_error(&self) -> Result<()> {
        match self.shared.take_error() {
            Some(e) => Err(e),
            None => Ok(()),
        }
    }

    /// Queue a block of text lines for parsing and sorting.
    pub fn push_block(&mut self, block: LineBlock) -> Result<()> {
        let Some(tx) = &self.block_tx else {
            return Err(KiraError::IoPlain(std::io::Error::other(
                "sorter already finished",
            )));
        };
        if tx.send(block).is_err() {
            self.check_error()?;
            return Err(KiraError::IoPlain(std::io::Error::other(
                "sorter parser threads stopped",
            )));
        }
        self.check_error()
    }

    /// Queue an already-parsed chunk (entries must reference `chunk.data`).
    pub fn push_chunk(&mut self, chunk: ParsedChunk) -> Result<()> {
        self.pushed_records += chunk.entries.len() as u64;
        let Some(tx) = &self.chunk_tx else {
            return Err(KiraError::IoPlain(std::io::Error::other(
                "sorter already finished",
            )));
        };
        if tx.send(chunk).is_err() {
            self.check_error()?;
            return Err(KiraError::IoPlain(std::io::Error::other(
                "sorter collector stopped",
            )));
        }
        self.check_error()
    }

    /// Finish input and produce the sorted stream.
    pub fn finish(mut self) -> Result<(SortedStream, SortMetrics)> {
        drop(self.block_tx.take());
        for p in self.parsers.drain(..) {
            let _ = p.join();
        }
        drop(self.chunk_tx.take());
        let collected = self
            .collector
            .take()
            .and_then(|h| h.join().ok())
            .ok_or_else(|| {
                KiraError::IoPlain(std::io::Error::other("sorter collector thread failed"))
            })?;
        let run_out = self
            .run_writer
            .take()
            .map(|h| {
                h.join().map_err(|_| {
                    KiraError::IoPlain(std::io::Error::other("run writer thread panicked"))
                })
            })
            .transpose()?
            .transpose()?
            .ok_or_else(|| {
                KiraError::IoPlain(std::io::Error::other("run writer thread missing"))
            })?;
        self.check_error()?;

        let mut metrics = SortMetrics {
            records: self.shared.records.load(Ordering::Relaxed),
            runs: run_out.paths.len() as u64,
            temp_bytes: run_out.temp_bytes,
            run_seconds: run_out.run_seconds,
            merge_pass_seconds: 0.0,
            merge_passes: 0,
            peak_records_buffered: self.shared.peak_buffered.load(Ordering::Relaxed),
        };

        let ctx = make_ctx(&self.dict, &self.cfg, &self.cols);
        let t0 = Instant::now();
        let entries = sort_chunks(&collected.remaining, &ctx, self.cfg.threads);
        let memory = MemorySource::new(collected.remaining, entries);
        metrics.run_seconds += t0.elapsed().as_secs_f64();

        if run_out.paths.is_empty() {
            return Ok((
                SortedStream {
                    inner: StreamInner::Memory(memory),
                    _temp: None,
                },
                metrics,
            ));
        }

        let temp = run_out.temp.clone();
        let fan_in = self.fan_in();
        let mut paths = run_out.paths;
        let mut mem_source = Some(memory);
        // Intermediate passes while too many runs remain.
        while paths.len() + usize::from(mem_source.as_ref().is_some_and(|m| !m.is_empty())) > fan_in
        {
            let t1 = Instant::now();
            let temp_ref = temp.as_ref().ok_or_else(|| {
                KiraError::IoPlain(std::io::Error::other("temporary directory missing"))
            })?;
            // Spill the in-memory run so that every source is a file.
            if let Some(m) = mem_source.take()
                && !m.is_empty()
            {
                let path = temp_ref.next_path("run");
                let bytes = write_source_to_run(Source::Memory(m), &path, &self.cfg)?;
                temp_ref.add_bytes(bytes);
                metrics.temp_bytes += bytes;
                paths.push(path);
            }
            let groups: Vec<Vec<PathBuf>> = paths.chunks(fan_in).map(|g| g.to_vec()).collect();
            log::info!(
                "intermediate merge pass: {} runs -> {} runs (fan-in {fan_in})",
                paths.len(),
                groups.len()
            );
            let cfg = &self.cfg;
            let ctx_ref = &ctx;
            let new_paths: Vec<Result<(PathBuf, u64)>> = groups
                .into_par_iter()
                .map(|group| {
                    let out_path = temp_ref.next_path("merge");
                    let mut sources = Vec::with_capacity(group.len());
                    for p in &group {
                        sources.push(Source::File(RunReader::open(p)?));
                    }
                    let mut merger = Merger::new(sources, ctx_ref.clone());
                    let bytes = write_merger_to_run(&mut merger, &out_path, cfg)?;
                    for p in &group {
                        temp_ref.remove(p);
                    }
                    Ok((out_path, bytes))
                })
                .collect();
            paths = Vec::new();
            for r in new_paths {
                let (p, bytes) = r?;
                temp_ref.add_bytes(bytes);
                metrics.temp_bytes += bytes;
                paths.push(p);
            }
            metrics.merge_passes += 1;
            metrics.merge_pass_seconds += t1.elapsed().as_secs_f64();
        }
        let mut sources: Vec<Source> = Vec::with_capacity(paths.len() + 1);
        for p in &paths {
            sources.push(Source::File(RunReader::open(p)?));
        }
        if let Some(m) = mem_source.take()
            && !m.is_empty()
        {
            sources.push(Source::Memory(m));
        }
        log::info!("final merge of {} sources", sources.len());
        let merger = Merger::new(sources, ctx);
        Ok((
            SortedStream {
                inner: StreamInner::Merge(merger),
                _temp: temp,
            },
            metrics,
        ))
    }

    fn fan_in(&self) -> usize {
        let per_source = (self.cfg.run_block_size * 2 + 64 * 1024) as u64;
        let by_memory = (self.cfg.budget.merge / per_source) as usize;
        by_memory.clamp(2, self.cfg.max_fan_in.max(2))
    }
}

impl Drop for ExternalSorter {
    fn drop(&mut self) {
        drop(self.block_tx.take());
        drop(self.chunk_tx.take());
        for p in self.parsers.drain(..) {
            let _ = p.join();
        }
        if let Some(c) = self.collector.take() {
            let _ = c.join();
        }
        if let Some(w) = self.run_writer.take() {
            let _ = w.join();
        }
    }
}

fn make_ctx(dict: &ChromDict, cfg: &SortConfig, cols: &ColumnMap) -> SortKeyContext {
    SortKeyContext::new(Arc::new(dict.ranks()), Arc::new(cfg.extra_cols.clone()))
        .with_pair_type_col(cols.pair_type)
}

/// Parse one text block into a chunk.
pub fn parse_block(
    block: LineBlock,
    cols: &ColumnMap,
    dict: &ChromDict,
    name: &str,
    ends: &mut Vec<u32>,
) -> Result<ParsedChunk> {
    let data = block.data;
    let mut entries = Vec::with_capacity(block.n_lines as usize);
    let mut start = 0usize;
    for (i, p) in memchr::memchr_iter(b'\n', &data).enumerate() {
        let mut line = &data[start..p];
        if line.last() == Some(&b'\r') {
            line = &line[..line.len() - 1];
        }
        let lineno = block.first_line + i as u64;
        if !line.is_empty() {
            let key = parse_line(line, cols, dict, lineno, ends, || {
                Location::file(name).at_line(lineno)
            })?;
            entries.push(SortEntry {
                key,
                chunk: 0,
                off: start as u32,
                len: line.len() as u32,
            });
        }
        start = p + 1;
    }
    Ok(ParsedChunk { data, entries })
}

/// Concatenate and sort the entries of a run.
fn sort_chunks(chunks: &[ParsedChunk], ctx: &SortKeyContext, threads: usize) -> Vec<SortEntry> {
    let total: usize = chunks.iter().map(|c| c.entries.len()).sum();
    let mut entries = Vec::with_capacity(total);
    for (ci, c) in chunks.iter().enumerate() {
        for e in &c.entries {
            entries.push(SortEntry {
                chunk: ci as u32,
                ..*e
            });
        }
    }
    let cmp = |a: &SortEntry, b: &SortEntry| {
        ctx.cmp_full(
            &a.key,
            chunks[a.chunk as usize].line(a),
            &b.key,
            chunks[b.chunk as usize].line(b),
        )
    };
    if threads > 1 && entries.len() > 50_000 {
        entries.par_sort_unstable_by(cmp);
    } else {
        entries.sort_unstable_by(cmp);
    }
    entries
}

/// Serialise a sorted run to disk with parallel block compression.
fn write_run(
    path: &std::path::Path,
    chunks: &[ParsedChunk],
    entries: &[SortEntry],
    cfg: &SortConfig,
) -> Result<u64> {
    let mut writer = RunWriter::create(path, cfg.compress_runs)?;
    let threads = cfg.threads.max(1);
    if !cfg.compress_runs || threads == 1 {
        let mut enc = BlockEncoder::new(cfg.run_block_size);
        for e in entries {
            if enc.push(&e.key, chunks[e.chunk as usize].line(e)) {
                let (raw, n) = enc.take();
                writer.write_raw_block(&raw, n)?;
            }
        }
        if !enc.is_empty() {
            let (raw, n) = enc.take();
            writer.write_raw_block(&raw, n)?;
        }
        let (_, _, bytes) = writer.finish()?;
        return Ok(bytes);
    }
    let (raw_tx, raw_rx) = bounded::<(u64, Vec<u8>, u32)>(threads * 2);
    let (comp_tx, comp_rx) = bounded::<(u64, Vec<u8>, u32, u32)>(threads * 2);
    let result: Result<u64> = std::thread::scope(|s| {
        let mut workers = Vec::new();
        for _ in 0..threads {
            let rx = raw_rx.clone();
            let tx = comp_tx.clone();
            workers.push(s.spawn(move || {
                while let Ok((i, raw, n)) = rx.recv() {
                    let c = compress_payload(&raw);
                    if tx.send((i, c, raw.len() as u32, n)).is_err() {
                        break;
                    }
                }
            }));
        }
        drop(raw_rx);
        drop(comp_tx);
        let sink = s.spawn(move || -> Result<u64> {
            let mut pending: BTreeMap<u64, (Vec<u8>, u32, u32)> = BTreeMap::new();
            let mut next = 0u64;
            while let Ok((i, c, raw_len, n)) = comp_rx.recv() {
                pending.insert(i, (c, raw_len, n));
                while let Some((c, raw_len, n)) = pending.remove(&next) {
                    writer.write_block(&c, raw_len, n)?;
                    next += 1;
                }
            }
            let (_, _, bytes) = writer.finish()?;
            Ok(bytes)
        });
        let mut enc = BlockEncoder::new(cfg.run_block_size);
        let mut idx = 0u64;
        let mut send_err = false;
        for e in entries {
            if enc.push(&e.key, chunks[e.chunk as usize].line(e)) {
                let (raw, n) = enc.take();
                if raw_tx.send((idx, raw, n)).is_err() {
                    send_err = true;
                    break;
                }
                idx += 1;
            }
        }
        if !send_err && !enc.is_empty() {
            let (raw, n) = enc.take();
            let _ = raw_tx.send((idx, raw, n));
        }
        drop(raw_tx);
        for w in workers {
            let _ = w.join();
        }
        sink.join()
            .map_err(|_| KiraError::IoPlain(std::io::Error::other("run sink thread panicked")))?
    });
    result
}

fn write_merger_to_run<S: RunSource>(
    merger: &mut Merger<S>,
    path: &std::path::Path,
    cfg: &SortConfig,
) -> Result<u64> {
    let mut writer = RunWriter::create(path, cfg.compress_runs)?;
    let mut enc = BlockEncoder::new(cfg.run_block_size);
    merger.for_each(|k, l| {
        if enc.push(k, l) {
            let (raw, n) = enc.take();
            writer.write_raw_block(&raw, n)?;
        }
        Ok(())
    })?;
    if !enc.is_empty() {
        let (raw, n) = enc.take();
        writer.write_raw_block(&raw, n)?;
    }
    let (_, _, bytes) = writer.finish()?;
    Ok(bytes)
}

fn write_source_to_run(
    mut source: Source,
    path: &std::path::Path,
    cfg: &SortConfig,
) -> Result<u64> {
    let mut writer = RunWriter::create(path, cfg.compress_runs)?;
    let mut enc = BlockEncoder::new(cfg.run_block_size);
    while !source.is_done() {
        if enc.push(source.key(), source.line()) {
            let (raw, n) = enc.take();
            writer.write_raw_block(&raw, n)?;
        }
        source.advance()?;
    }
    if !enc.is_empty() {
        let (raw, n) = enc.take();
        writer.write_raw_block(&raw, n)?;
    }
    let (_, _, bytes) = writer.finish()?;
    Ok(bytes)
}

enum StreamInner {
    Memory(MemorySource),
    Merge(Merger<Source>),
}

/// The sorted output: a lending iterator over `(key, line)`.
pub struct SortedStream {
    inner: StreamInner,
    _temp: Option<Arc<TempManager>>,
}

impl SortedStream {
    /// Next record in sorted order.
    #[inline]
    pub fn next_record(&mut self) -> Result<Option<(&PairKey, &[u8])>> {
        match &mut self.inner {
            StreamInner::Memory(m) => m_next(m),
            StreamInner::Merge(m) => m.next_record(),
        }
    }

    /// Consume all records through a callback.
    pub fn for_each<F: FnMut(&PairKey, &[u8]) -> Result<()>>(&mut self, mut f: F) -> Result<()> {
        while let Some((k, l)) = self.next_record()? {
            f(k, l)?;
        }
        Ok(())
    }
}

#[inline]
fn m_next(m: &mut MemorySource) -> Result<Option<(&PairKey, &[u8])>> {
    // MemorySource has no pending-advance protocol, so emulate it with an
    // explicit cursor step: advance first unless this is the first call.
    if m.started() {
        m.advance()?;
    } else {
        m.mark_started();
    }
    if m.is_done() {
        return Ok(None);
    }
    Ok(Some((m.key(), m.line())))
}
