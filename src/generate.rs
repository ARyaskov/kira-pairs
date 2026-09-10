//! Deterministic synthetic `.pairs` generator for benchmarks and tests.

use std::io::Write;

use crate::error::Result;
use crate::pairs::header::Header;

/// Generator parameters.
#[derive(Debug, Clone)]
pub struct GenerateConfig {
    /// Number of records.
    pub records: u64,
    /// RNG seed.
    pub seed: u64,
    /// Number of chromosomes (`chr1..chrN`).
    pub chromosomes: usize,
    /// Length of each chromosome.
    pub chrom_length: u64,
    /// Fraction of cis pairs.
    pub cis_fraction: f64,
    /// Fraction of records that are near-duplicates of a previous record.
    pub duplicate_rate: f64,
    /// Maximum per-side offset of a duplicate (bp).
    pub duplicate_radius: u64,
    /// Fraction of records with at least one unmapped side.
    pub unmapped_fraction: f64,
    /// Number of extra integer columns.
    pub extra_columns: usize,
    /// Read ID length.
    pub readid_length: usize,
    /// Emit records in block-sorted, upper-triangular order.
    pub sorted: bool,
}

impl Default for GenerateConfig {
    fn default() -> Self {
        Self {
            records: 1_000_000,
            seed: 42,
            chromosomes: 24,
            chrom_length: 100_000_000,
            cis_fraction: 0.75,
            duplicate_rate: 0.1,
            duplicate_radius: 2,
            unmapped_fraction: 0.02,
            extra_columns: 0,
            readid_length: 24,
            sorted: false,
        }
    }
}

/// SplitMix64 / xoshiro256** PRNG (no external dependency, stable output).
#[derive(Debug, Clone)]
pub struct Rng {
    s: [u64; 4],
}

impl Rng {
    /// Seeded generator.
    pub fn new(seed: u64) -> Self {
        let mut x = seed;
        let mut next = || {
            x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
            let mut z = x;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
            z ^ (z >> 31)
        };
        Self {
            s: [next(), next(), next(), next()],
        }
    }

    /// Next 64 random bits.
    #[inline]
    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// Uniform float in `[0, 1)`.
    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform integer in `[0, n)`.
    #[inline]
    pub fn below(&mut self, n: u64) -> u64 {
        if n == 0 { 0 } else { self.next_u64() % n }
    }
}

const READID_ALPHABET: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

/// Column names produced by the generator.
pub fn columns(extra: usize) -> Vec<String> {
    let mut c: Vec<String> = crate::pairs::header::STANDARD_COLUMNS
        .iter()
        .map(|s| s.to_string())
        .collect();
    for i in 0..extra {
        c.push(format!("extra{}", i + 1));
    }
    c
}

/// Header for generated data.
pub fn header(cfg: &GenerateConfig) -> Header {
    let mut cs = crate::chroms::ChromSizes::default();
    for i in 0..cfg.chromosomes {
        cs.push(&format!("chr{}", i + 1), cfg.chrom_length);
    }
    let mut h = Header::standard(
        Some("synthetic"),
        Some(&cs),
        &columns(cfg.extra_columns),
        if cfg.sorted {
            "upper triangle"
        } else {
            "whole matrix"
        },
    );
    if cfg.sorted {
        let _ = h.mark_sorted();
    }
    h
}

/// Write a synthetic dataset.
pub fn generate<W: Write>(cfg: &GenerateConfig, out: &mut W) -> Result<u64> {
    let h = header(cfg);
    h.write_to(out)?;
    let mut rng = Rng::new(cfg.seed);
    let n_chrom = cfg.chromosomes.max(1) as u64;
    let mut recent: Vec<(u64, u64, u64, u64, u8, u8)> = Vec::with_capacity(1024);
    let mut records: Vec<Vec<u8>> = Vec::new();
    let mut line = Vec::with_capacity(256);
    let mut written = 0u64;
    let mut readid = vec![b'A'; cfg.readid_length.max(1)];
    for i in 0..cfg.records {
        line.clear();
        // Read ID: deterministic pseudo-random string.
        for b in readid.iter_mut() {
            *b = READID_ALPHABET[rng.below(READID_ALPHABET.len() as u64) as usize];
        }
        line.extend_from_slice(&readid);
        let r = rng.next_f64();
        let (c1, p1, c2, p2, s1, s2, pt): (u64, u64, u64, u64, u8, u8, &[u8]);
        if r < cfg.unmapped_fraction {
            if rng.next_f64() < 0.5 {
                c1 = u64::MAX;
                p1 = 0;
                c2 = u64::MAX;
                p2 = 0;
                s1 = b'-';
                s2 = b'-';
                pt = b"NN";
            } else {
                c1 = u64::MAX;
                p1 = 0;
                c2 = rng.below(n_chrom);
                p2 = 1 + rng.below(cfg.chrom_length);
                s1 = b'-';
                s2 = if rng.next_f64() < 0.5 { b'+' } else { b'-' };
                pt = b"NU";
            }
        } else if !recent.is_empty() && rng.next_f64() < cfg.duplicate_rate {
            let src = recent[rng.below(recent.len() as u64) as usize];
            let jitter = |rng: &mut Rng, p: u64| {
                let d =
                    rng.below(2 * cfg.duplicate_radius + 1) as i64 - cfg.duplicate_radius as i64;
                (p as i64 + d).max(1) as u64
            };
            c1 = src.0;
            p1 = jitter(&mut rng, src.1);
            c2 = src.2;
            p2 = jitter(&mut rng, src.3);
            s1 = src.4;
            s2 = src.5;
            pt = b"UU";
        } else {
            let a = rng.below(n_chrom);
            let pa = 1 + rng.below(cfg.chrom_length);
            let (b, pb) = if rng.next_f64() < cfg.cis_fraction {
                // Cis with a heavy-tailed distance distribution.
                let d = (10f64.powf(rng.next_f64() * 6.5)) as u64;
                (a, (pa + d).min(cfg.chrom_length))
            } else {
                (rng.below(n_chrom), 1 + rng.below(cfg.chrom_length))
            };
            let (c1_, p1_, c2_, p2_) = if (a, pa) <= (b, pb) {
                (a, pa, b, pb)
            } else {
                (b, pb, a, pa)
            };
            c1 = c1_;
            p1 = p1_;
            c2 = c2_;
            p2 = p2_;
            s1 = if rng.next_f64() < 0.5 { b'+' } else { b'-' };
            s2 = if rng.next_f64() < 0.5 { b'+' } else { b'-' };
            pt = if rng.next_f64() < 0.9 { b"UU" } else { b"UR" };
            if recent.len() < 1024 {
                recent.push((c1, p1, c2, p2, s1, s2));
            } else {
                let slot = rng.below(1024) as usize;
                recent[slot] = (c1, p1, c2, p2, s1, s2);
            }
        }
        let push_chrom = |line: &mut Vec<u8>, c: u64| {
            if c == u64::MAX {
                line.push(b'!');
            } else {
                line.extend_from_slice(b"chr");
                crate::util::int::write_u64(line, c + 1);
            }
        };
        line.push(b'\t');
        push_chrom(&mut line, c1);
        line.push(b'\t');
        crate::util::int::write_u64(&mut line, p1);
        line.push(b'\t');
        push_chrom(&mut line, c2);
        line.push(b'\t');
        crate::util::int::write_u64(&mut line, p2);
        line.push(b'\t');
        line.push(s1);
        line.push(b'\t');
        line.push(s2);
        line.push(b'\t');
        line.extend_from_slice(pt);
        for _ in 0..cfg.extra_columns {
            line.push(b'\t');
            crate::util::int::write_u64(&mut line, rng.below(61));
        }
        line.push(b'\n');
        if cfg.sorted {
            records.push(line.clone());
        } else {
            out.write_all(&line)?;
            written += 1;
        }
        let _ = i;
    }
    if cfg.sorted {
        // Sort with pairtools semantics (lexicographic chrom, numeric pos,
        // pair type, then stable).
        let dict = crate::chroms::ChromDict::new();
        let cols = crate::pairs::columns::ColumnMap::from_names(columns(cfg.extra_columns))?;
        let mut ends = Vec::new();
        let mut keyed: Vec<(crate::pairs::record::PairKey, usize)> =
            Vec::with_capacity(records.len());
        for (i, r) in records.iter().enumerate() {
            let l = &r[..r.len() - 1];
            let k = crate::pairs::record::parse_line(
                l,
                &cols,
                &dict,
                i as u64,
                &mut ends,
                Default::default,
            )?;
            keyed.push((k, i));
        }
        let ctx = crate::sort::key::SortKeyContext::new(
            std::sync::Arc::new(dict.ranks()),
            std::sync::Arc::new(Vec::new()),
        );
        keyed.sort_by(|a, b| ctx.cmp_full(&a.0, &records[a.1], &b.0, &records[b.1]));
        for (_, i) in keyed {
            out.write_all(&records[i])?;
            written += 1;
        }
    }
    Ok(written)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic() {
        let cfg = GenerateConfig {
            records: 2000,
            extra_columns: 2,
            ..Default::default()
        };
        let mut a = Vec::new();
        let mut b = Vec::new();
        generate(&cfg, &mut a).unwrap();
        generate(&cfg, &mut b).unwrap();
        assert_eq!(a, b);
        let text = String::from_utf8(a).unwrap();
        let body: Vec<&str> = text.lines().filter(|l| !l.starts_with('#')).collect();
        assert_eq!(body.len(), 2000);
        assert!(body.iter().all(|l| l.split('\t').count() == 10));
        assert!(body.iter().any(|l| l.contains("\t!\t")));
    }

    #[test]
    fn sorted_output_is_sorted() {
        let cfg = GenerateConfig {
            records: 500,
            sorted: true,
            ..Default::default()
        };
        let mut a = Vec::new();
        generate(&cfg, &mut a).unwrap();
        let text = String::from_utf8(a).unwrap();
        assert!(text.contains("#sorted: chr1-chr2-pos1-pos2"));
        let keys: Vec<(String, String, u64, u64)> = text
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                (
                    f[1].to_string(),
                    f[3].to_string(),
                    f[2].parse().unwrap(),
                    f[4].parse().unwrap(),
                )
            })
            .collect();
        assert!(keys.windows(2).all(|w| w[0] <= w[1]));
    }
}
