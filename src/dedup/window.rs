//! Positional window of recent records with a `pos2` bucket index for dense
//! regions.

use std::collections::{HashMap, VecDeque};

/// A record in the window.
#[derive(Debug, Clone)]
pub struct Entry {
    /// Absolute record index.
    pub idx: u64,
    /// `pos1` of the record.
    pub pos1: u64,
    /// `pos2` of the record.
    pub pos2: u64,
    /// Signature of the non-positional match fields.
    pub sig: u64,
    /// Whether the record was classified as a duplicate (greedy mode).
    pub is_dup: bool,
    /// Extra-column values (only when extra column pairs are configured).
    pub extra: Option<Box<[u8]>>,
    /// Read ID (only when parent ids are requested in greedy mode).
    pub readid: Option<Box<[u8]>>,
}

/// A candidate neighbour returned by [`Window::candidates`].
#[derive(Debug, Clone, Copy)]
pub struct Neighbor {
    /// Absolute record index.
    pub idx: u64,
    /// `pos1` of the candidate.
    pub pos1: u64,
    /// `pos2` of the candidate.
    pub pos2: u64,
    /// Duplicate flag of the candidate.
    pub is_dup: bool,
    /// Slot in the window (valid until the next mutation).
    pub slot: usize,
}

const INDEX_THRESHOLD: usize = 48;

/// Records whose `pos1` is within `radius` of the current position.
#[derive(Debug)]
pub struct Window {
    radius: u64,
    cell: u64,
    entries: VecDeque<Entry>,
    buckets: HashMap<u64, Vec<u64>>,
    indexed: bool,
    stale: usize,
}

impl Window {
    /// Window with the given mismatch radius.
    pub fn new(radius: u64) -> Self {
        Self {
            radius,
            cell: radius.saturating_add(1).max(1),
            entries: VecDeque::new(),
            buckets: HashMap::new(),
            indexed: false,
            stale: 0,
        }
    }

    /// Number of records in the window.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True when empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Entry at a slot returned by [`Window::candidates`].
    #[inline]
    pub fn entry(&self, slot: usize) -> &Entry {
        &self.entries[slot]
    }

    /// Drop everything (new chromosome block).
    pub fn clear(&mut self) {
        self.entries.clear();
        self.buckets.clear();
        self.indexed = false;
        self.stale = 0;
    }

    /// Drop records with `pos1 < cur - radius`.
    pub fn expire(&mut self, cur_pos1: u64) {
        let min = cur_pos1.saturating_sub(self.radius);
        while let Some(front) = self.entries.front() {
            if front.pos1 < min {
                self.entries.pop_front();
                self.stale += 1;
            } else {
                break;
            }
        }
        if self.entries.is_empty() && self.indexed {
            self.buckets.clear();
            self.indexed = false;
            self.stale = 0;
        }
    }

    /// Add a record.
    pub fn insert(&mut self, entry: Entry) {
        let (pos2, idx) = (entry.pos2, entry.idx);
        self.entries.push_back(entry);
        if self.indexed {
            self.buckets.entry(pos2 / self.cell).or_default().push(idx);
        } else if self.entries.len() > INDEX_THRESHOLD {
            self.build_index();
        }
    }

    fn build_index(&mut self) {
        self.buckets.clear();
        for e in &self.entries {
            self.buckets
                .entry(e.pos2 / self.cell)
                .or_default()
                .push(e.idx);
        }
        self.indexed = true;
        self.stale = 0;
    }

    /// Collect candidates with matching signature and `|pos2 - pos2'| <= radius`,
    /// in increasing index order.
    pub fn candidates(&mut self, pos2: u64, sig: u64, out: &mut Vec<Neighbor>) {
        out.clear();
        if self.entries.is_empty() {
            return;
        }
        if !self.indexed {
            for (slot, e) in self.entries.iter().enumerate() {
                if e.sig == sig && e.pos2.abs_diff(pos2) <= self.radius {
                    out.push(Neighbor {
                        idx: e.idx,
                        pos1: e.pos1,
                        pos2: e.pos2,
                        is_dup: e.is_dup,
                        slot,
                    });
                }
            }
            return;
        }
        if self.stale > self.entries.len() * 4 + 1024 {
            self.build_index();
        }
        let b = pos2 / self.cell;
        let front = self.entries.front().map(|e| e.idx).unwrap_or(0);
        for bucket in b.saturating_sub(1)..=b.saturating_add(1) {
            let Some(list) = self.buckets.get_mut(&bucket) else {
                continue;
            };
            list.retain(|i| *i >= front);
            let entries = &self.entries;
            for &idx in list.iter() {
                let slot = (idx - front) as usize;
                let e = &entries[slot];
                if e.sig == sig && e.pos2.abs_diff(pos2) <= self.radius {
                    out.push(Neighbor {
                        idx: e.idx,
                        pos1: e.pos1,
                        pos2: e.pos2,
                        is_dup: e.is_dup,
                        slot,
                    });
                }
            }
        }
        out.sort_unstable_by_key(|n| n.idx);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn e(idx: u64, pos1: u64, pos2: u64, sig: u64) -> Entry {
        Entry {
            idx,
            pos1,
            pos2,
            sig,
            is_dup: false,
            extra: None,
            readid: None,
        }
    }

    #[test]
    fn linear_and_indexed_agree() {
        let mut w = Window::new(3);
        let mut out = Vec::new();
        for i in 0..200u64 {
            w.insert(e(i, 100, i * 2, 7));
        }
        w.candidates(50, 7, &mut out);
        let ids: Vec<u64> = out.iter().map(|n| n.idx).collect();
        assert_eq!(ids, vec![24, 25, 26]);
        w.candidates(50, 8, &mut out);
        assert!(out.is_empty());
        w.expire(104);
        assert!(w.is_empty());
        let mut w = Window::new(3);
        for i in 0..10u64 {
            w.insert(e(i, 100 + i, i * 2, 7));
        }
        w.candidates(6, 7, &mut out);
        let ids: Vec<u64> = out.iter().map(|n| n.idx).collect();
        assert_eq!(ids, vec![2, 3, 4]);
        w.expire(106);
        assert_eq!(w.len(), 7);
    }
}
