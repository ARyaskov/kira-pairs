//! Ordered parallel parsing of a `.pairs` body: blocks are parsed on worker
//! threads and handed back in input order.

use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};
use std::thread::JoinHandle;

use crossbeam_channel::{Receiver, bounded};

use crate::chroms::ChromDict;
use crate::error::{KiraError, Result};
use crate::pairs::columns::ColumnMap;
use crate::pairs::reader::BodyBlocks;
use crate::sort::external::parse_block;
use crate::sort::key::ParsedChunk;

/// Yields parsed chunks in input order while parsing in parallel.
pub struct OrderedParser {
    rx: Receiver<(u64, ParsedChunk)>,
    pending: BTreeMap<u64, ParsedChunk>,
    next: u64,
    error: Arc<Mutex<Option<KiraError>>>,
    threads: Vec<JoinHandle<()>>,
    bytes_read: Arc<std::sync::atomic::AtomicU64>,
    done: bool,
}

impl OrderedParser {
    /// Start reading `body` with `threads` parser workers.
    pub fn new(
        body: BodyBlocks,
        cols: Arc<ColumnMap>,
        dict: Arc<ChromDict>,
        threads: usize,
    ) -> Self {
        Self::with_depth(body, cols, dict, threads, threads.max(1) * 2)
    }

    /// Like [`OrderedParser::new`] with an explicit channel depth (see
    /// [`crate::memory::MemoryBudget::channel_depth`]).
    pub fn with_depth(
        mut body: BodyBlocks,
        cols: Arc<ColumnMap>,
        dict: Arc<ChromDict>,
        threads: usize,
        depth: usize,
    ) -> Self {
        let threads = threads.max(1);
        let depth = depth.max(1);
        let (block_tx, block_rx) = bounded(depth);
        let (out_tx, out_rx) = bounded(depth);
        let error = Arc::new(Mutex::new(None));
        let bytes_read = Arc::new(std::sync::atomic::AtomicU64::new(0));
        let mut handles = Vec::new();
        let name = body.name().to_string();
        let err_r = Arc::clone(&error);
        let bytes_r = Arc::clone(&bytes_read);
        handles.push(std::thread::spawn(move || {
            loop {
                match body.next_block() {
                    Ok(Some(b)) => {
                        bytes_r.store(body.bytes_read(), std::sync::atomic::Ordering::Relaxed);
                        if block_tx.send(b).is_err() {
                            break;
                        }
                    }
                    Ok(None) => break,
                    Err(e) => {
                        let mut g = err_r.lock().unwrap_or_else(|p| p.into_inner());
                        if g.is_none() {
                            *g = Some(e);
                        }
                        break;
                    }
                }
            }
        }));
        for _ in 0..threads {
            let rx = block_rx.clone();
            let tx = out_tx.clone();
            let cols = Arc::clone(&cols);
            let dict = Arc::clone(&dict);
            let err = Arc::clone(&error);
            let name = name.clone();
            handles.push(std::thread::spawn(move || {
                let mut ends = Vec::with_capacity(32);
                while let Ok(block) = rx.recv() {
                    let idx = block.index;
                    match parse_block(block, &cols, &dict, &name, &mut ends) {
                        Ok(chunk) => {
                            if tx.send((idx, chunk)).is_err() {
                                break;
                            }
                        }
                        Err(e) => {
                            let mut g = err.lock().unwrap_or_else(|p| p.into_inner());
                            if g.is_none() {
                                *g = Some(e);
                            }
                            break;
                        }
                    }
                }
            }));
        }
        drop(out_tx);
        Self {
            rx: out_rx,
            pending: BTreeMap::new(),
            next: 0,
            error,
            threads: handles,
            bytes_read,
            done: false,
        }
    }

    /// Bytes of decompressed input consumed so far.
    pub fn bytes_read(&self) -> u64 {
        self.bytes_read.load(std::sync::atomic::Ordering::Relaxed)
    }

    fn take_error(&self) -> Option<KiraError> {
        self.error.lock().unwrap_or_else(|p| p.into_inner()).take()
    }

    /// Next chunk in input order.
    pub fn next_chunk(&mut self) -> Result<Option<ParsedChunk>> {
        loop {
            if let Some(c) = self.pending.remove(&self.next) {
                self.next += 1;
                return Ok(Some(c));
            }
            if self.done {
                if let Some(e) = self.take_error() {
                    return Err(e);
                }
                return Ok(None);
            }
            match self.rx.recv() {
                Ok((i, c)) => {
                    self.pending.insert(i, c);
                }
                Err(_) => {
                    self.done = true;
                }
            }
        }
    }
}

impl Drop for OrderedParser {
    fn drop(&mut self) {
        let (_tx, rx) = bounded(0);
        self.rx = rx;
        for h in self.threads.drain(..) {
            let _ = h.join();
        }
    }
}
