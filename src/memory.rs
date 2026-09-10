//! Parsing and budgeting of the `--memory` option.

use crate::error::{KiraError, Result};

/// Parse a human-readable size such as `512M`, `8G`, `1024`, `2GiB`, `1.5G`.
///
/// Suffixes are case-insensitive and interpreted as binary multiples
/// (K = 1024, M = 1024², G = 1024³, T = 1024⁴), mirroring GNU `sort -S`.
pub fn parse_size(s: &str) -> Result<u64> {
    let s = s.trim();
    if s.is_empty() {
        return Err(KiraError::arg("empty size"));
    }
    let split = s
        .find(|c: char| !(c.is_ascii_digit() || c == '.'))
        .unwrap_or(s.len());
    let (num, suffix) = s.split_at(split);
    let value: f64 = num
        .parse()
        .map_err(|_| KiraError::arg(format!("invalid size {s:?}")))?;
    let suffix = suffix.trim().to_ascii_lowercase();
    let suffix = suffix
        .strip_suffix("ib")
        .or_else(|| suffix.strip_suffix('b'))
        .unwrap_or(&suffix);
    let mult: f64 = match suffix {
        "" => 1.0,
        "k" => 1024.0,
        "m" => 1024.0 * 1024.0,
        "g" => 1024.0 * 1024.0 * 1024.0,
        "t" => 1024.0 * 1024.0 * 1024.0 * 1024.0,
        _ => return Err(KiraError::arg(format!("unknown size suffix in {s:?}"))),
    };
    let bytes = value * mult;
    if !bytes.is_finite() || bytes < 0.0 || bytes > u64::MAX as f64 {
        return Err(KiraError::arg(format!("size out of range {s:?}")));
    }
    Ok(bytes as u64)
}

/// Format a byte count with a binary suffix for diagnostics.
pub fn format_size(bytes: u64) -> String {
    const UNITS: [&str; 5] = ["B", "KiB", "MiB", "GiB", "TiB"];
    let mut value = bytes as f64;
    let mut unit = 0;
    while value >= 1024.0 && unit < UNITS.len() - 1 {
        value /= 1024.0;
        unit += 1;
    }
    if unit == 0 {
        format!("{bytes} B")
    } else {
        format!("{value:.1} {}", UNITS[unit])
    }
}

/// How a total memory budget is split between the components of a job.
///
/// The split is deliberately conservative: the bulk goes to sortable
/// record buffers, and fixed reserves are kept for channel blocks,
/// compression buffers and merge read-ahead so that the process stays
/// close to the requested limit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MemoryBudget {
    /// Total budget requested by the user.
    pub total: u64,
    /// Bytes available for buffered records (sort runs, bin tables).
    pub records: u64,
    /// Bytes reserved for I/O blocks in flight through channels.
    pub io_buffers: u64,
    /// Bytes reserved for compression/decompression work buffers.
    pub compression: u64,
    /// Bytes reserved for merge-phase read-ahead buffers.
    pub merge: u64,
}

impl MemoryBudget {
    /// Smallest budget accepted; below this the sorter cannot make progress.
    pub const MIN_TOTAL: u64 = 64 * 1024 * 1024;

    /// Input block size for parallel parsing: shrinks below the 4 MiB default
    /// when the I/O reserve cannot hold one block per worker plus the channel
    /// depth (each block is resident about twice while it is parsed).
    pub fn block_size(&self, threads: usize) -> usize {
        let default = crate::io::buffered::DEFAULT_BLOCK_SIZE as u64;
        let slots = (threads as u64 + 2) * 3;
        let per_slot = self.io_buffers / slots.max(1);
        (per_slot.clamp(256 * 1024, default)) as usize
    }

    /// Depth of the bounded reader -> parser -> consumer channels so that the
    /// blocks in flight (two channels plus one block per worker, each block
    /// resident twice while it is parsed) fit into the I/O buffer reserve.
    pub fn channel_depth(&self, threads: usize) -> usize {
        let block = self.block_size(threads) as u64 * 2;
        let in_flight = self.io_buffers / block;
        let depth = in_flight.saturating_sub(threads as u64) / 2;
        (depth as usize).clamp(1, threads.max(1) * 2)
    }

    /// Derive a budget from a total and the thread counts that scale buffers.
    pub fn new(total: u64, threads: usize, io_threads: usize) -> Result<Self> {
        if total < Self::MIN_TOTAL {
            return Err(KiraError::Memory(format!(
                "at least {} is required, got {}",
                format_size(Self::MIN_TOTAL),
                format_size(total)
            )));
        }
        // Blocks in flight: reader -> parsers -> collector. Each stage may hold
        // a couple of blocks per thread.
        let block = crate::io::buffered::DEFAULT_BLOCK_SIZE as u64;
        let io_buffers = (block * (threads as u64 * 3 + io_threads as u64 * 2)).min(total / 4);
        // Compression workers each keep an input and an output block.
        let compression = (block * 2 * (threads as u64 + 1)).min(total / 8);
        // Merge read-ahead: fan-in × block, bounded to a quarter of the budget.
        let merge = (total / 4).max(16 * 1024 * 1024);
        let reserved = io_buffers + compression + merge;
        let records = total.saturating_sub(reserved);
        // Never let the record budget collapse below a workable minimum.
        let records = records.max(total / 3);
        Ok(Self {
            total,
            records,
            io_buffers,
            compression,
            merge,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_sizes() {
        assert_eq!(parse_size("1024").unwrap(), 1024);
        assert_eq!(parse_size("1K").unwrap(), 1024);
        assert_eq!(parse_size("8G").unwrap(), 8 << 30);
        assert_eq!(parse_size("512M").unwrap(), 512 << 20);
        assert_eq!(parse_size("2GiB").unwrap(), 2 << 30);
        assert_eq!(parse_size("1.5g").unwrap(), 3 << 29);
        assert_eq!(parse_size("2gb").unwrap(), 2 << 30);
        assert!(parse_size("x").is_err());
        assert!(parse_size("1X").is_err());
        assert!(parse_size("").is_err());
    }

    #[test]
    fn budget_is_split_sensibly() {
        let b = MemoryBudget::new(1 << 30, 8, 2).unwrap();
        assert!(b.records > b.total / 3);
        assert!(b.records + b.io_buffers + b.compression + b.merge <= b.total + b.total / 3);
        assert!(MemoryBudget::new(1 << 20, 1, 1).is_err());
    }

    #[test]
    fn formats_sizes() {
        assert_eq!(format_size(512), "512 B");
        assert_eq!(format_size(1536), "1.5 KiB");
    }
}
