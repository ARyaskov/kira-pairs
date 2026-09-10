//! Union-find over a sliding range of record indices.

/// Disjoint-set forest whose elements are consecutive record indices. Roots
/// are always the smallest index of their set, which makes "first record in
/// input order" the natural cluster representative.
#[derive(Debug, Default)]
pub struct UnionFind {
    base: u64,
    parent: Vec<u64>,
    last_pos1: Vec<u64>,
    block: Vec<u64>,
    size: Vec<u64>,
}

impl UnionFind {
    /// Empty forest.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add element `idx` (must be `base + len`). Returns `idx`.
    pub fn push(&mut self, idx: u64, pos1: u64, block: u64) -> u64 {
        debug_assert_eq!(idx, self.base + self.parent.len() as u64);
        self.parent.push(idx);
        self.last_pos1.push(pos1);
        self.block.push(block);
        self.size.push(1);
        idx
    }

    #[inline]
    fn slot(&self, idx: u64) -> usize {
        (idx - self.base) as usize
    }

    /// Root of `idx` with path compression.
    pub fn find(&mut self, idx: u64) -> u64 {
        let mut root = idx;
        while self.parent[self.slot(root)] != root {
            root = self.parent[self.slot(root)];
        }
        let mut cur = idx;
        while self.parent[self.slot(cur)] != root {
            let s = self.slot(cur);
            let next = self.parent[s];
            self.parent[s] = root;
            cur = next;
        }
        root
    }

    /// Merge the sets of `a` and `b`; the smaller root wins.
    pub fn union(&mut self, a: u64, b: u64) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        let (root, child) = if ra < rb { (ra, rb) } else { (rb, ra) };
        let (rs, cs) = (self.slot(root), self.slot(child));
        self.parent[cs] = root;
        self.last_pos1[rs] = self.last_pos1[rs].max(self.last_pos1[cs]);
        self.size[rs] += self.size[cs];
    }

    /// Largest `pos1` in the set rooted at `root`.
    #[inline]
    pub fn last_pos1(&self, root: u64) -> u64 {
        self.last_pos1[self.slot(root)]
    }

    /// Block id of the set rooted at `root`.
    #[inline]
    pub fn block_of(&self, root: u64) -> u64 {
        self.block[self.slot(root)]
    }

    /// Number of members of the set rooted at `root`.
    #[inline]
    pub fn size(&self, root: u64) -> u64 {
        self.size[self.slot(root)]
    }

    /// Forget all elements below `new_base` (they must all be emitted and no
    /// remaining element may point at them).
    pub fn truncate_before(&mut self, new_base: u64) {
        if new_base <= self.base {
            return;
        }
        let n = (new_base - self.base) as usize;
        if n >= self.parent.len() {
            self.parent.clear();
            self.last_pos1.clear();
            self.block.clear();
            self.size.clear();
        } else {
            // Only safe when no remaining element references a dropped
            // root; callers guarantee this by draining whole clusters.
            let keep_ok = self.parent[n..].iter().all(|p| *p >= new_base);
            if !keep_ok {
                return;
            }
            self.parent.drain(..n);
            self.last_pos1.drain(..n);
            self.block.drain(..n);
            self.size.drain(..n);
        }
        self.base = new_base;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unions_and_roots() {
        let mut uf = UnionFind::new();
        for i in 0..6 {
            uf.push(i, i * 10, 1);
        }
        uf.union(1, 3);
        uf.union(3, 5);
        assert_eq!(uf.find(5), 1);
        assert_eq!(uf.size(1), 3);
        assert_eq!(uf.last_pos1(1), 50);
        uf.union(0, 5);
        assert_eq!(uf.find(3), 0);
        assert_eq!(uf.size(0), 4);
        assert_eq!(uf.find(2), 2);
        uf.truncate_before(6);
        uf.push(6, 60, 2);
        assert_eq!(uf.find(6), 6);
        assert_eq!(uf.block_of(6), 2);
    }
}
