//! Fixed-resolution contact binning with bounded memory.
//!
//! Each accepted pair is mapped to `(bin1, bin2)` in cooler's global bin
//! numbering (chromosomes in `chrom.sizes` order, bins of `resolution` bp,
//! 1-based `.pairs` coordinates converted to 0-based). Counts are aggregated
//! in a hash table that is spilled to sorted temporary runs when it grows
//! past its budget; runs are merge-reduced at the end.

pub mod aggregate;

use std::collections::HashSet;
use std::path::PathBuf;
use std::sync::Arc;

use crate::chroms::{ChromDict, ChromSizes};
use crate::error::{KiraError, Result};
use crate::io::temp::TempManager;
use crate::pairs::record::PairKey;
use aggregate::{Aggregator, BinTable};

/// Binning configuration.
#[derive(Debug, Clone)]
pub struct BinConfig {
    /// One or more resolutions (bp).
    pub resolutions: Vec<u64>,
    /// Chromosome sizes defining the bin table.
    pub chromsizes: ChromSizes,
    /// Minimal MAPQ on both sides (requires `mapq1`/`mapq2` columns).
    pub min_mapq: Option<u64>,
    /// Accepted pair types (`None` = all mapped pairs).
    pub pair_types: Option<HashSet<Vec<u8>>>,
    /// Treat positions as 0-based (default: 1-based `.pairs`).
    pub zero_based: bool,
    /// Memory budget for in-memory aggregation (all resolutions together).
    pub memory_bytes: u64,
    /// Temporary directory for spilled runs.
    pub tmpdir: Option<PathBuf>,
}

/// Bin offsets for one resolution.
#[derive(Debug, Clone)]
pub struct BinLayout {
    /// Resolution in bp.
    pub resolution: u64,
    /// First global bin id of each chromosome (chromsizes order), plus a
    /// final entry equal to the total number of bins.
    pub offsets: Vec<u64>,
}

impl BinLayout {
    /// Layout for a resolution.
    pub fn new(cs: &ChromSizes, resolution: u64) -> Self {
        let mut offsets = Vec::with_capacity(cs.len() + 1);
        let mut acc = 0u64;
        for size in cs.sizes() {
            offsets.push(acc);
            acc += size.div_ceil(resolution);
        }
        offsets.push(acc);
        Self {
            resolution,
            offsets,
        }
    }

    /// Total number of bins.
    pub fn n_bins(&self) -> u64 {
        *self.offsets.last().unwrap_or(&0)
    }

    /// Number of bins of chromosome `ci`.
    pub fn chrom_bins(&self, ci: usize) -> u64 {
        self.offsets[ci + 1] - self.offsets[ci]
    }

    /// Map a global bin id back to `(chrom index, start, end)`.
    pub fn locate(&self, bin: u64, cs: &ChromSizes) -> Option<(usize, u64, u64)> {
        let ci = self.offsets.partition_point(|&o| o <= bin).checked_sub(1)?;
        if ci + 1 >= self.offsets.len() {
            return None;
        }
        let local = bin - self.offsets[ci];
        let start = local * self.resolution;
        let end = (start + self.resolution).min(cs.sizes()[ci]);
        Some((ci, start, end))
    }
}

/// Counters collected by the binner.
#[derive(Debug, Default, Clone)]
pub struct BinMetrics {
    /// Pairs accepted for binning.
    pub accepted: u64,
    /// Pairs rejected by filters (pair type, MAPQ, unmapped).
    pub rejected: u64,
    /// Pairs on chromosomes absent from the chromosome sizes.
    pub unknown_chrom: u64,
    /// Pairs with positions beyond the chromosome length (clamped).
    pub out_of_range: u64,
}

/// Streaming multi-resolution binner.
pub struct Binner {
    cfg: BinConfig,
    dict: Arc<ChromDict>,
    chrom_index: Vec<Option<u32>>,
    layouts: Vec<BinLayout>,
    aggs: Vec<Aggregator>,
    unmapped_id: Option<u32>,
    metrics: BinMetrics,
    _temp: Option<Arc<TempManager>>,
}

impl Binner {
    /// Create a binner; spills go to a job directory under `cfg.tmpdir`.
    pub fn new(cfg: BinConfig, dict: Arc<ChromDict>) -> Result<Self> {
        if cfg.resolutions.is_empty() {
            return Err(KiraError::arg("at least one --resolution is required"));
        }
        for r in &cfg.resolutions {
            if *r == 0 {
                return Err(KiraError::arg("resolution must be positive"));
            }
        }
        let temp = Arc::new(TempManager::new(cfg.tmpdir.as_deref(), "kira-pairs-bin")?);
        let per_res = (cfg.memory_bytes / cfg.resolutions.len() as u64).max(8 << 20);
        let layouts: Vec<BinLayout> = cfg
            .resolutions
            .iter()
            .map(|r| BinLayout::new(&cfg.chromsizes, *r))
            .collect();
        let aggs = layouts
            .iter()
            .map(|l| Aggregator::new(l.resolution, per_res, Arc::clone(&temp)))
            .collect();
        let unmapped_id = dict.get(crate::chroms::UNMAPPED_CHROM);
        Ok(Self {
            cfg,
            dict,
            chrom_index: Vec::new(),
            layouts,
            aggs,
            unmapped_id,
            metrics: BinMetrics::default(),
            _temp: Some(temp),
        })
    }

    /// Bin layouts per resolution.
    pub fn layouts(&self) -> &[BinLayout] {
        &self.layouts
    }

    /// Metrics so far.
    pub fn metrics(&self) -> &BinMetrics {
        &self.metrics
    }

    #[inline]
    fn chrom_of(&mut self, id: u32) -> Option<u32> {
        let i = id as usize;
        if i >= self.chrom_index.len() {
            let n = self.dict.len().max(i + 1);
            for j in self.chrom_index.len()..n {
                let name = self.dict.name(j as u32);
                let name = String::from_utf8_lossy(&name);
                self.chrom_index
                    .push(self.cfg.chromsizes.index_of(&name).map(|x| x as u32));
            }
            if self.unmapped_id.is_none() {
                self.unmapped_id = self.dict.get(crate::chroms::UNMAPPED_CHROM);
            }
        }
        self.chrom_index[i]
    }

    /// Observe a pair. `mapq` carries `(mapq1, mapq2)` when available.
    #[inline]
    pub fn observe(&mut self, key: &PairKey, mapq: Option<(u64, u64)>) -> Result<()> {
        if let Some(u) = self.unmapped_id
            && (key.chrom1 == u || key.chrom2 == u)
        {
            self.metrics.rejected += 1;
            return Ok(());
        }
        if let Some(pts) = &self.cfg.pair_types
            && !pts.contains(key.pair_type_bytes())
        {
            self.metrics.rejected += 1;
            return Ok(());
        }
        if let Some(min) = self.cfg.min_mapq {
            match mapq {
                Some((a, b)) if a >= min && b >= min => {}
                _ => {
                    self.metrics.rejected += 1;
                    return Ok(());
                }
            }
        }
        let (Some(c1), Some(c2)) = (self.chrom_of(key.chrom1), self.chrom_of(key.chrom2)) else {
            self.metrics.unknown_chrom += 1;
            return Ok(());
        };
        let (p1, p2) = if self.cfg.zero_based {
            (key.pos1, key.pos2)
        } else {
            if key.pos1 == 0 || key.pos2 == 0 {
                self.metrics.rejected += 1;
                return Ok(());
            }
            (key.pos1 - 1, key.pos2 - 1)
        };
        let mut oor = false;
        for (li, layout) in self.layouts.iter().enumerate() {
            let res = layout.resolution;
            let mut local1 = p1 / res;
            let mut local2 = p2 / res;
            let n1 = layout.chrom_bins(c1 as usize);
            let n2 = layout.chrom_bins(c2 as usize);
            if local1 >= n1 {
                local1 = n1.saturating_sub(1);
                oor = true;
            }
            if local2 >= n2 {
                local2 = n2.saturating_sub(1);
                oor = true;
            }
            let mut b1 = layout.offsets[c1 as usize] + local1;
            let mut b2 = layout.offsets[c2 as usize] + local2;
            if b1 > b2 {
                std::mem::swap(&mut b1, &mut b2);
            }
            self.aggs[li].add(b1, b2)?;
        }
        if oor {
            self.metrics.out_of_range += 1;
        }
        self.metrics.accepted += 1;
        Ok(())
    }

    /// Finish aggregation and return one sorted table per resolution.
    pub fn finish(self) -> Result<(Vec<(BinLayout, BinTable)>, BinMetrics)> {
        let mut out = Vec::with_capacity(self.aggs.len());
        for (layout, agg) in self.layouts.into_iter().zip(self.aggs) {
            out.push((layout, agg.finish()?));
        }
        if self.metrics.out_of_range > 0 {
            log::warn!(
                "{} pairs had positions beyond the chromosome length and were clamped to the last bin",
                self.metrics.out_of_range
            );
        }
        if self.metrics.unknown_chrom > 0 {
            log::warn!(
                "{} pairs on chromosomes absent from the chromosome sizes were skipped",
                self.metrics.unknown_chrom
            );
        }
        Ok((out, self.metrics))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn layout_offsets() {
        let cs = ChromSizes::parse("chr1\t25\nchr2\t10\n", "t").unwrap();
        let l = BinLayout::new(&cs, 10);
        assert_eq!(l.offsets, vec![0, 3, 4]);
        assert_eq!(l.n_bins(), 4);
        assert_eq!(l.locate(2, &cs), Some((0, 20, 25)));
        assert_eq!(l.locate(3, &cs), Some((1, 0, 10)));
        assert_eq!(l.locate(4, &cs), None);
    }

    fn key(c1: u32, p1: u64, c2: u32, p2: u64, pt: &[u8; 2]) -> PairKey {
        let mut k = PairKey {
            seq: 0,
            pos1: p1,
            pos2: p2,
            chrom1: c1,
            chrom2: c2,
            pair_type: [0; 8],
            strand1: b'+',
            strand2: b'+',
            flags: 0,
            pair_type_len: 2,
        };
        k.pair_type[..2].copy_from_slice(pt);
        k
    }

    #[test]
    fn bins_and_conserves_counts() {
        let cs = ChromSizes::parse("chr1\t25\nchr2\t10\n", "t").unwrap();
        let dict = Arc::new(ChromDict::with_names(["!", "chr1", "chr2", "chrU"]));
        let cfg = BinConfig {
            resolutions: vec![10, 5],
            chromsizes: cs.clone(),
            min_mapq: None,
            pair_types: Some([b"UU".to_vec()].into_iter().collect()),
            zero_based: false,
            memory_bytes: 64 << 20,
            tmpdir: None,
        };
        let mut b = Binner::new(cfg, dict).unwrap();
        b.observe(&key(1, 1, 1, 11, b"UU"), None).unwrap(); // bins (0,1) @10
        b.observe(&key(1, 11, 1, 1, b"UU"), None).unwrap(); // flipped -> (0,1)
        b.observe(&key(1, 21, 2, 3, b"UU"), None).unwrap(); // (2,3)
        b.observe(&key(1, 21, 2, 3, b"UR"), None).unwrap(); // rejected
        b.observe(&key(0, 0, 2, 3, b"UU"), None).unwrap(); // unmapped
        b.observe(&key(3, 1, 2, 3, b"UU"), None).unwrap(); // unknown chrom
        b.observe(&key(1, 1, 1, 30, b"UU"), None).unwrap(); // pos beyond -> clamped (0,2)
        let (tables, m) = b.finish().unwrap();
        assert_eq!(m.accepted, 4);
        assert_eq!(m.rejected, 2);
        assert_eq!(m.unknown_chrom, 1);
        assert_eq!(m.out_of_range, 1);
        let mut t10 = tables[0].1.clone_for_test();
        let mut rows = Vec::new();
        while let Some(r) = t10.next_row().unwrap() {
            rows.push(r);
        }
        assert_eq!(rows, vec![(0, 1, 2), (0, 2, 1), (2, 3, 1)]);
        let total: u64 = rows.iter().map(|r| r.2).sum();
        assert_eq!(total, m.accepted);
    }
}
