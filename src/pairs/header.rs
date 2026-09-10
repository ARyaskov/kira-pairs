//! `.pairs` header model with pairtools-compatible editing operations.

use std::io::Write;

use crate::chroms::ChromSizes;
use crate::error::{KiraError, Result};

/// Version string written into new headers.
pub const PAIRS_FORMAT_VERSION: &str = "1.0.0";
/// Value pairtools writes into `#sorted:` after block sorting.
pub const SORTED_VALUE: &str = "chr1-chr2-pos1-pos2";
/// Standard `.pairs` column names in order.
pub const STANDARD_COLUMNS: [&str; 8] = [
    "readID",
    "chrom1",
    "pos1",
    "chrom2",
    "pos2",
    "strand1",
    "strand2",
    "pair_type",
];

/// A `.pairs` header: an ordered list of `#`-prefixed lines.
///
/// Lines are stored verbatim (without the trailing newline) so unknown
/// fields and comments round-trip unchanged.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Header {
    lines: Vec<String>,
}

impl Header {
    /// Header from raw lines (each must start with `#`).
    pub fn from_lines<I, S>(lines: I) -> Result<Self>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let mut out = Vec::new();
        for l in lines {
            let l: String = l.into();
            let l = l.trim_end_matches(['\n', '\r']).to_string();
            if !l.starts_with('#') {
                return Err(KiraError::header(format!(
                    "header line does not start with '#': {l:?}"
                )));
            }
            out.push(l);
        }
        Ok(Self { lines: out })
    }

    /// A pairtools-style standard header.
    pub fn standard(
        assembly: Option<&str>,
        chromsizes: Option<&ChromSizes>,
        columns: &[String],
        shape: &str,
    ) -> Self {
        let mut lines = vec![
            format!("## pairs format v{PAIRS_FORMAT_VERSION}"),
            format!("#shape: {shape}"),
            format!("#genome_assembly: {}", assembly.unwrap_or("unknown")),
        ];
        if let Some(cs) = chromsizes {
            for (name, size) in cs.iter() {
                lines.push(format!("#chromsize: {name} {size}"));
            }
        }
        lines.push(format!("#columns: {}", columns.join(" ")));
        Self { lines }
    }

    /// Raw header lines.
    pub fn lines(&self) -> &[String] {
        &self.lines
    }

    /// True when there are no header lines at all.
    pub fn is_empty(&self) -> bool {
        self.lines.is_empty()
    }

    /// pairtools `is_empty_header`: no lines or no leading `##` line.
    pub fn is_valid_pairs_header(&self) -> bool {
        self.lines.first().is_some_and(|l| l.starts_with("##"))
    }

    /// Values of all occurrences of a field, e.g. `field = "chromsize"`.
    pub fn fields(&self, field: &str) -> Vec<&str> {
        let prefix = format!("{field}:");
        self.lines
            .iter()
            .filter_map(|l| {
                let body = l.trim_start_matches('#');
                body.strip_prefix(&prefix).map(str::trim_start)
            })
            .collect()
    }

    /// First value of a field.
    pub fn field(&self, field: &str) -> Option<&str> {
        self.fields(field).into_iter().next()
    }

    /// Column names from `#columns:`; empty when absent.
    pub fn columns(&self) -> Vec<String> {
        self.field("columns")
            .map(|v| {
                v.split(' ')
                    .filter(|s| !s.is_empty())
                    .map(String::from)
                    .collect()
            })
            .unwrap_or_default()
    }

    /// Replace the `#columns:` line (or append one).
    pub fn set_columns(&mut self, columns: &[String]) {
        let new = format!("#columns: {}", columns.join(" "));
        let mut found = false;
        for l in &mut self.lines {
            if l.starts_with("#columns:") {
                *l = new.clone();
                found = true;
            }
        }
        if !found {
            self.lines.push(new);
        }
    }

    /// Append columns to `#columns:` (pairtools `append_columns`).
    pub fn append_columns(&mut self, extra: &[String]) {
        for l in &mut self.lines {
            if l.starts_with("#columns: ") {
                for c in extra {
                    l.push(' ');
                    l.push_str(c);
                }
            }
        }
    }

    /// `#chromsize:` entries in header order.
    pub fn chromsizes(&self) -> Result<ChromSizes> {
        let mut cs = ChromSizes::default();
        for v in self.fields("chromsize") {
            let mut it = v.split_whitespace();
            let (Some(name), Some(size)) = (it.next(), it.next()) else {
                return Err(KiraError::header(format!(
                    "malformed #chromsize line: {v:?}"
                )));
            };
            let size: u64 = size
                .parse()
                .map_err(|_| KiraError::header(format!("invalid chromosome size: {v:?}")))?;
            if !cs.push(name, size) {
                return Err(KiraError::header(format!(
                    "duplicate #chromsize entry {name}"
                )));
            }
        }
        Ok(cs)
    }

    /// True when a `#sorted:` line is present.
    pub fn is_sorted(&self) -> bool {
        self.lines.iter().any(|l| l.starts_with("#sorted"))
    }

    /// pairtools `mark_header_as_sorted`.
    pub fn mark_sorted(&mut self) -> Result<()> {
        if !self.is_valid_pairs_header() {
            return Err(KiraError::header(
                "input file is not valid .pairs: header is empty or lacks '## pairs format'",
            ));
        }
        if !self.is_sorted() {
            let line = format!("#sorted: {SORTED_VALUE}");
            if self.lines[0].starts_with("##") {
                self.lines.insert(1, line);
            } else {
                self.lines.insert(0, line);
            }
        }
        for l in &mut self.lines {
            if let Some(rest) = l.strip_prefix("#chromosomes:") {
                // pairtools slices the line at a fixed offset and keeps a
                // stray ':' token; we sort the actual names instead.
                let mut chroms: Vec<&str> = rest.split_whitespace().collect();
                chroms.sort_unstable();
                *l = format!("#chromosomes: {}", chroms.join(" "));
            }
        }
        Ok(())
    }

    /// pairtools `_update_header_entry`: replace or insert `#field: value`.
    pub fn set_field(&mut self, field: &str, value: &str) {
        let newline = format!("#{field}: {value}");
        let prefix = format!("#{field}");
        let mut found = false;
        for l in &mut self.lines {
            if l.starts_with(&prefix) {
                *l = newline.clone();
                found = true;
            }
        }
        if !found {
            if self.lines.last().is_some_and(|l| l.starts_with("#columns")) {
                let n = self.lines.len();
                self.lines.insert(n - 1, newline);
            } else {
                self.lines.push(newline);
            }
        }
    }

    /// `#samheader:` payloads in order.
    pub fn samheader(&self) -> Vec<String> {
        self.fields("samheader")
            .into_iter()
            .map(String::from)
            .collect()
    }

    /// pairtools `append_new_pg`: add a `@PG` record to every `@PG` chain of
    /// the embedded SAM header, re-inserting all `#samheader:` lines right
    /// before `#columns:` (this reorders the header the same way pairtools
    /// does).
    ///
    /// When the SAM header has no `@PG` lines, nothing is added (pairtools
    /// behaviour). Chains whose parent cannot be resolved are treated as new
    /// chains instead of failing.
    pub fn append_pg(&mut self, id: &str, pn: &str, cl: &str, vn: &str) -> Result<()> {
        if !self.is_valid_pairs_header() {
            return Err(KiraError::header(
                "input file is not valid .pairs: header is empty or lacks '## pairs format'",
            ));
        }
        let samheader = self.samheader();
        let other: Vec<String> = self
            .lines
            .iter()
            .filter(|l| !l.trim_start_matches('#').starts_with("samheader:"))
            .cloned()
            .collect();
        let new_sam = add_pg_to_samheader(&samheader, id, pn, cl, vn);
        let mut lines: Vec<String> = other
            .iter()
            .filter(|l| !l.starts_with("#columns"))
            .cloned()
            .collect();
        lines.extend(new_sam.iter().map(|l| format!("#samheader: {l}")));
        lines.extend(other.iter().filter(|l| l.starts_with("#columns")).cloned());
        self.lines = lines;
        Ok(())
    }

    /// pairtools `insert_samheader`: place SAM header lines before `#columns`.
    pub fn insert_samheader(&mut self, samheader: &[String]) {
        let mut lines: Vec<String> = self
            .lines
            .iter()
            .filter(|l| !l.starts_with("#columns"))
            .cloned()
            .collect();
        lines.extend(samheader.iter().map(|l| format!("#samheader: {l}")));
        lines.extend(
            self.lines
                .iter()
                .filter(|l| l.starts_with("#columns"))
                .cloned(),
        );
        self.lines = lines;
    }

    /// pairtools `subset_chroms_in_pairsheader`.
    pub fn subset_chromosomes(&mut self, keep: &[String]) {
        let keep: std::collections::HashSet<&str> = keep.iter().map(String::as_str).collect();
        let mut out = Vec::with_capacity(self.lines.len());
        for l in &self.lines {
            if l.starts_with("#chromsize:") {
                let name = l.split_whitespace().nth(1).unwrap_or("");
                if keep.contains(name) {
                    out.push(l.clone());
                }
            } else if l.starts_with("#chromosomes:") {
                let names: Vec<&str> = l
                    .split_whitespace()
                    .skip(1)
                    .filter(|c| keep.contains(c))
                    .collect();
                out.push(format!("#chromosomes: {}", names.join(" ")));
            } else {
                out.push(l.clone());
            }
        }
        self.lines = out;
    }

    /// Serialise to a writer, one line each with `\n`.
    pub fn write_to<W: Write>(&self, w: &mut W) -> std::io::Result<()> {
        for l in &self.lines {
            w.write_all(l.as_bytes())?;
            w.write_all(b"\n")?;
        }
        Ok(())
    }

    /// Serialise to bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut v = Vec::new();
        // Writing to a Vec cannot fail.
        let _ = self.write_to(&mut v);
        v
    }
}

#[derive(Debug, Clone)]
struct PgRecord {
    id: String,
    pp: Option<String>,
    raw: String,
}

fn parse_pg(line: &str) -> Option<PgRecord> {
    let mut id = None;
    let mut pp = None;
    for tv in line.split('\t').skip(1) {
        let (tag, value) = tv.split_once(':')?;
        match tag {
            "ID" => id = Some(value.to_string()),
            "PP" => pp = Some(value.to_string()),
            _ => {}
        }
    }
    Some(PgRecord {
        id: id?,
        pp,
        raw: line.to_string(),
    })
}

fn parse_pg_chains(samheader: &[String]) -> Vec<Vec<PgRecord>> {
    let mut pending: Vec<PgRecord> = samheader
        .iter()
        .filter(|l| l.starts_with("@PG"))
        .filter_map(|l| parse_pg(l))
        .collect();
    let mut chains: Vec<Vec<PgRecord>> = Vec::new();
    while !pending.is_empty() {
        let mut placed = false;
        for i in 0..pending.len() {
            match &pending[i].pp {
                None => {
                    let pg = pending.remove(i);
                    chains.push(vec![pg]);
                    placed = true;
                    break;
                }
                Some(pp) => {
                    if let Some(chain) = chains
                        .iter_mut()
                        .find(|c| c.last().is_some_and(|l| &l.id == pp))
                    {
                        let pg = pending.remove(i);
                        chain.push(pg);
                        placed = true;
                        break;
                    }
                }
            }
        }
        if !placed {
            // Parent unresolvable: start a new chain (pairtools `force=True`).
            let pg = pending.remove(0);
            chains.push(vec![pg]);
        }
    }
    chains
}

fn format_pg(id: &str, pn: &str, cl: &str, pp: &str, vn: &str) -> String {
    format!("@PG\tID:{id}\tPN:{pn}\tCL:{cl}\tPP:{pp}\tVN:{vn}")
}

fn add_pg_to_samheader(
    samheader: &[String],
    id: &str,
    pn: &str,
    cl: &str,
    vn: &str,
) -> Vec<String> {
    let is_pre = |l: &String| l.starts_with("@HD") || l.starts_with("@SQ") || l.starts_with("@RG");
    let pre: Vec<String> = samheader
        .iter()
        .filter(|l| is_pre(l))
        .map(|l| l.trim().to_string())
        .collect();
    let post: Vec<String> = samheader
        .iter()
        .filter(|l| !is_pre(l) && !l.starts_with("@PG"))
        .map(|l| l.trim().to_string())
        .collect();
    let mut chains = parse_pg_chains(samheader);
    let n_chains = chains.len();
    for (i, chain) in chains.iter_mut().enumerate() {
        let pp = chain.last().map(|l| l.id.clone()).unwrap_or_default();
        let new_id = if n_chains > 1 {
            format!("{id}-{}.{}", i + 1, chain.len() + 1)
        } else {
            id.to_string()
        };
        let raw = format_pg(&new_id, pn, cl, &pp, vn);
        chain.push(PgRecord {
            id: new_id,
            pp: Some(pp),
            raw,
        });
    }
    let mut out = pre;
    out.extend(chains.into_iter().flatten().map(|pg| pg.raw));
    out.extend(post);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Header {
        Header::from_lines([
            "## pairs format v1.0.0",
            "#shape: upper triangle",
            "#genome_assembly: hg38",
            "#chromsize: chr1 100",
            "#chromsize: chr2 50",
            "#chromosomes: chr2 chr1",
            "#samheader: @SQ\tSN:chr1\tLN:100",
            "#samheader: @PG\tID:bwa\tPN:bwa\tVN:1",
            "#custom_field: keep me",
            "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2",
        ])
        .unwrap()
    }

    #[test]
    fn extracts_fields() {
        let h = sample();
        assert_eq!(h.columns().len(), 10);
        assert_eq!(h.field("shape"), Some("upper triangle"));
        assert_eq!(h.chromsizes().unwrap().size_of("chr2"), Some(50));
        assert_eq!(h.field("custom_field"), Some("keep me"));
    }

    #[test]
    fn marks_sorted_like_pairtools() {
        let mut h = sample();
        h.mark_sorted().unwrap();
        assert_eq!(h.lines()[1], "#sorted: chr1-chr2-pos1-pos2");
        assert!(h.lines().iter().any(|l| l == "#chromosomes: chr1 chr2"));
        let n = h.lines().len();
        h.mark_sorted().unwrap();
        assert_eq!(h.lines().len(), n);
    }

    #[test]
    fn appends_pg_and_reorders_samheader() {
        let mut h = sample();
        h.append_pg("kira-pairs_sort", "kira-pairs", "kira-pairs sort", "0.1.0")
            .unwrap();
        let lines = h.lines();
        // columns is last, samheader lines right before it, custom field kept.
        assert!(lines.last().unwrap().starts_with("#columns:"));
        let n = lines.len();
        assert!(lines[n - 2].starts_with("#samheader: @PG\tID:kira-pairs_sort\tPN:kira-pairs\tCL:kira-pairs sort\tPP:bwa\tVN:0.1.0"));
        assert!(lines.iter().any(|l| l == "#custom_field: keep me"));
        assert_eq!(h.columns().len(), 10);
    }

    #[test]
    fn no_pg_added_without_existing_pg() {
        let mut h =
            Header::from_lines(["## pairs format v1.0.0", "#columns: readID chrom1"]).unwrap();
        h.append_pg("x", "x", "x", "1").unwrap();
        assert_eq!(h.lines().len(), 2);
    }

    #[test]
    fn set_field_inserts_before_columns() {
        let mut h = Header::from_lines(["## pairs format v1.0.0", "#columns: a b"]).unwrap();
        h.set_field("shape", "upper triangle");
        assert_eq!(h.lines()[1], "#shape: upper triangle");
        h.set_field("shape", "whole matrix");
        assert_eq!(h.lines().len(), 3);
        assert_eq!(h.field("shape"), Some("whole matrix"));
    }

    #[test]
    fn rejects_non_comment_lines() {
        assert!(Header::from_lines(["not a header"]).is_err());
    }
}
