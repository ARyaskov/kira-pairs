//! Temporary run-file management with predictable cleanup.

use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::error::{KiraError, Result};

/// A scratch directory holding private run files for one job.
///
/// Files are created with collision-safe names and removed when the manager
/// is dropped (or earlier via [`TempManager::remove`]).
#[derive(Debug)]
pub struct TempManager {
    dir: tempfile::TempDir,
    counter: AtomicU64,
    live: Mutex<Vec<PathBuf>>,
    bytes_written: AtomicU64,
}

impl TempManager {
    /// Create a job directory under `base` (or the system temp dir).
    pub fn new(base: Option<&Path>, prefix: &str) -> Result<Self> {
        let mut builder = tempfile::Builder::new();
        let prefix = format!("{prefix}-{}-", std::process::id());
        builder.prefix(&prefix);
        let dir = match base {
            Some(b) => {
                std::fs::create_dir_all(b).map_err(|e| KiraError::io(b, e))?;
                builder.tempdir_in(b).map_err(|e| KiraError::io(b, e))?
            }
            None => builder
                .tempdir()
                .map_err(|e| KiraError::io(std::env::temp_dir(), e))?,
        };
        log::debug!("temporary directory: {}", dir.path().display());
        Ok(Self {
            dir,
            counter: AtomicU64::new(0),
            live: Mutex::new(Vec::new()),
            bytes_written: AtomicU64::new(0),
        })
    }

    /// Path of the job directory.
    pub fn path(&self) -> &Path {
        self.dir.path()
    }

    /// Reserve a fresh unique file path with the given tag.
    pub fn next_path(&self, tag: &str) -> PathBuf {
        let n = self.counter.fetch_add(1, Ordering::Relaxed);
        let p = self.dir.path().join(format!("{tag}-{n:06}.kprun"));
        self.live
            .lock()
            .unwrap_or_else(|e| e.into_inner())
            .push(p.clone());
        p
    }

    /// Record bytes written to temporary storage (for metrics).
    pub fn add_bytes(&self, n: u64) {
        self.bytes_written.fetch_add(n, Ordering::Relaxed);
    }

    /// Total bytes written to temporary storage.
    pub fn bytes_written(&self) -> u64 {
        self.bytes_written.load(Ordering::Relaxed)
    }

    /// Remove one temporary file now (ignores missing files).
    pub fn remove(&self, path: &Path) {
        if let Err(e) = std::fs::remove_file(path)
            && e.kind() != std::io::ErrorKind::NotFound
        {
            log::warn!("could not remove temporary file {}: {e}", path.display());
        }
        let mut live = self.live.lock().unwrap_or_else(|e| e.into_inner());
        live.retain(|p| p != path);
    }

    /// Number of temporary files currently registered.
    pub fn live_count(&self) -> usize {
        self.live.lock().unwrap_or_else(|e| e.into_inner()).len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creates_and_cleans_up() {
        let dir;
        {
            let tm = TempManager::new(None, "kira-test").unwrap();
            dir = tm.path().to_path_buf();
            let p = tm.next_path("run");
            std::fs::write(&p, b"x").unwrap();
            assert_eq!(tm.live_count(), 1);
            tm.remove(&p);
            assert_eq!(tm.live_count(), 0);
            assert!(!p.exists());
            std::fs::write(tm.next_path("run"), b"y").unwrap();
        }
        assert!(!dir.exists());
    }
}
