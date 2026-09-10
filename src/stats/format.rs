//! Stats serialisation (pairtools TSV and YAML, JSON extension) and parsing
//! of stats files for merging.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use crate::error::{KiraError, Result};
use crate::stats::accumulator::{CIS_KB, DIRS, DistBins, StatsSnapshot};
use crate::util::pyfloat::format_py_float;

/// A scalar in the stats tree.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    /// Integer counter.
    Int(u64),
    /// Floating point summary.
    Float(f64),
    /// String.
    Str(String),
}

impl Value {
    fn render(&self) -> String {
        match self {
            Value::Int(i) => i.to_string(),
            Value::Float(f) => format_py_float(*f),
            Value::Str(s) => s.clone(),
        }
    }
}

/// Output format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum StatsFormat {
    /// pairtools `key<TAB>value` (flattened with `/`).
    #[default]
    Tsv,
    /// pairtools YAML (`--yaml`).
    Yaml,
    /// JSON (kira-pairs extension).
    Json,
}

/// A node of the nested stats tree (insertion ordered).
#[derive(Debug, Clone, PartialEq)]
pub enum Node {
    /// Leaf.
    Leaf(Value),
    /// Ordered map.
    Map(Vec<(String, Node)>),
}

impl Node {
    fn map() -> Self {
        Node::Map(Vec::new())
    }

    fn insert(&mut self, key: &str, node: Node) {
        if let Node::Map(m) = self {
            m.push((key.to_string(), node));
        }
    }

    /// Insert a `a/b/c` path.
    fn insert_path(&mut self, path: &str, value: Value) {
        let mut parts = path.split('/').peekable();
        let mut cur = self;
        while let Some(p) = parts.next() {
            if parts.peek().is_none() {
                cur.insert(p, Node::Leaf(value));
                return;
            }
            let Node::Map(m) = cur else { return };
            let pos = match m.iter().position(|(k, _)| k == p) {
                Some(i) => i,
                None => {
                    m.push((p.to_string(), Node::map()));
                    m.len() - 1
                }
            };
            cur = &mut m[pos].1;
        }
    }

    fn is_empty(&self) -> bool {
        match self {
            Node::Leaf(_) => false,
            Node::Map(m) => m.is_empty(),
        }
    }
}

/// Build the flat pairtools TSV rows for a snapshot.
pub fn flatten(s: &StatsSnapshot) -> Vec<(String, Value)> {
    let mut out: Vec<(String, Value)> = vec![
        ("total".into(), Value::Int(s.total)),
        ("total_unmapped".into(), Value::Int(s.total_unmapped)),
        (
            "total_single_sided_mapped".into(),
            Value::Int(s.total_single_sided_mapped),
        ),
        ("total_mapped".into(), Value::Int(s.total_mapped)),
        ("total_dups".into(), Value::Int(s.total_dups)),
        ("total_nodups".into(), Value::Int(s.total_nodups)),
        ("cis".into(), Value::Int(s.cis)),
        ("trans".into(), Value::Int(s.trans)),
    ];
    for (k, n) in &s.pair_types {
        out.push((format!("pair_types/{k}"), Value::Int(*n)));
    }
    for (i, kb) in CIS_KB.iter().enumerate() {
        out.push((format!("cis_{kb}kb+"), Value::Int(s.cis_kb[i])));
    }
    let sum = s.summary();
    for (k, v) in &sum.frac_cis {
        out.push((format!("summary/{k}"), v.clone()));
    }
    out.push(("summary/frac_dups".into(), sum.frac_dups.clone()));
    out.push((
        "summary/complexity_naive".into(),
        sum.complexity_naive.clone(),
    ));
    for (k, v) in &sum.convergence {
        out.push((format!("summary/dist_freq_convergence/{k}"), v.clone()));
    }
    for ((a, b), n) in &s.chrom_freq {
        out.push((format!("chrom_freq/{a}/{b}"), Value::Int(*n)));
    }
    let edges = s.bins.edges();
    for (i, &lo) in edges.iter().enumerate() {
        let range = if i + 1 < edges.len() {
            format!("{lo}-{}", edges[i + 1])
        } else {
            format!("{lo}+")
        };
        for (d, dir) in DIRS.iter().enumerate() {
            out.push((
                format!("dist_freq/{range}/{dir}"),
                Value::Int(s.dist_freq[d][i]),
            ));
        }
    }
    for (name, size) in &s.chromsizes {
        out.push((format!("chromsizes/{name}"), Value::Int(*size)));
    }
    out
}

/// Build the nested tree used for YAML/JSON output. When `pairtools_yaml`
/// is set, zero/empty values are omitted like pairtools' `format_yaml`.
pub fn tree(s: &StatsSnapshot, pairtools_yaml: bool) -> Node {
    let keep = |v: u64| !pairtools_yaml || v != 0;
    let mut root = Node::map();
    if !pairtools_yaml {
        root.insert("filter_expression", Node::Leaf(Value::Str(String::new())));
    }
    let scalars = [
        ("total", s.total),
        ("total_unmapped", s.total_unmapped),
        ("total_single_sided_mapped", s.total_single_sided_mapped),
        ("total_mapped", s.total_mapped),
        ("total_dups", s.total_dups),
        ("total_nodups", s.total_nodups),
        ("cis", s.cis),
        ("trans", s.trans),
    ];
    for (k, v) in scalars {
        if keep(v) {
            root.insert(k, Node::Leaf(Value::Int(v)));
        }
    }
    let mut pt = Node::map();
    for (k, n) in &s.pair_types {
        pt.insert(k, Node::Leaf(Value::Int(*n)));
    }
    if !pairtools_yaml || !pt.is_empty() {
        root.insert("pair_types", pt);
    }
    for (i, kb) in CIS_KB.iter().enumerate() {
        if keep(s.cis_kb[i]) {
            root.insert(&format!("cis_{kb}kb+"), Node::Leaf(Value::Int(s.cis_kb[i])));
        }
    }
    let sum = s.summary();
    let mut summary = Node::map();
    for (k, v) in &sum.frac_cis {
        summary.insert(k, Node::Leaf(v.clone()));
    }
    summary.insert("frac_dups", Node::Leaf(sum.frac_dups.clone()));
    summary.insert("complexity_naive", Node::Leaf(sum.complexity_naive.clone()));
    let mut conv = Node::map();
    for (k, v) in &sum.convergence {
        conv.insert_path(k, v.clone());
    }
    summary.insert("dist_freq_convergence", conv);
    root.insert("summary", summary);
    let mut cf = Node::map();
    for ((a, b), n) in &s.chrom_freq {
        cf.insert(&format!("{a}/{b}"), Node::Leaf(Value::Int(*n)));
    }
    if !pairtools_yaml || !cf.is_empty() {
        root.insert("chrom_freq", cf);
    }
    let mut df = Node::map();
    for (d, dir) in DIRS.iter().enumerate() {
        let mut m = Node::map();
        for (i, lo) in s.bins.edges().iter().enumerate() {
            m.insert(&lo.to_string(), Node::Leaf(Value::Int(s.dist_freq[d][i])));
        }
        df.insert(dir, m);
    }
    root.insert("dist_freq", df);
    let mut cs = Node::map();
    for (name, size) in &s.chromsizes {
        cs.insert(name, Node::Leaf(Value::Int(*size)));
    }
    if !pairtools_yaml || !cs.is_empty() {
        root.insert("chromsizes", cs);
    }
    root
}

/// Render a snapshot in the requested format.
pub fn render(s: &StatsSnapshot, format: StatsFormat) -> String {
    match format {
        StatsFormat::Tsv => {
            let mut out = String::new();
            for (k, v) in flatten(s) {
                let _ = writeln!(out, "{k}\t{}", v.render());
            }
            out
        }
        StatsFormat::Yaml => {
            let mut out = String::new();
            let mut root = Node::map();
            root.insert("no_filter", tree(s, true));
            write_yaml(&root, 0, &mut out, 0);
            out
        }
        StatsFormat::Json => {
            let mut root = Node::map();
            root.insert("no_filter", tree(s, false));
            let mut out = String::new();
            write_json(&root, 0, &mut out);
            out.push('\n');
            out
        }
    }
}

fn yaml_needs_quotes(s: &str) -> bool {
    if s.is_empty() {
        return true;
    }
    let lower = s.to_ascii_lowercase();
    if matches!(
        lower.as_str(),
        "true" | "false" | "null" | "~" | "yes" | "no" | "on" | "off" | ".nan" | ".inf" | "-.inf"
    ) {
        return true;
    }
    if s.parse::<f64>().is_ok() || s.parse::<i64>().is_ok() {
        return true;
    }
    let first = s.chars().next().unwrap_or(' ');
    if !(first.is_ascii_alphanumeric() || first == '_' || first == '/') {
        return true;
    }
    s.contains(": ")
        || s.contains(" #")
        || s.ends_with(':')
        || s.contains('\t')
        || s.contains('\n')
        || s.contains('\'')
        || s.contains('"')
}

fn yaml_scalar(v: &Value) -> String {
    match v {
        Value::Int(i) => i.to_string(),
        Value::Float(f) => {
            if f.is_nan() {
                ".nan".into()
            } else if f.is_infinite() {
                if *f > 0.0 {
                    ".inf".into()
                } else {
                    "-.inf".into()
                }
            } else {
                format_py_float(*f)
            }
        }
        Value::Str(s) => {
            if yaml_needs_quotes(s) {
                format!("'{}'", s.replace('\'', "''"))
            } else {
                s.clone()
            }
        }
    }
}

fn yaml_key(k: &str) -> String {
    if yaml_needs_quotes(k) {
        format!("'{}'", k.replace('\'', "''"))
    } else {
        k.to_string()
    }
}

/// `mode`: 0 = string keys, 1 = inside `dist_freq` (strand keys),
/// 2 = bin maps whose keys are integers.
fn write_yaml(node: &Node, indent: usize, out: &mut String, mode: u8) {
    let Node::Map(m) = node else { return };
    let pad = " ".repeat(indent);
    if m.is_empty() {
        let _ = writeln!(out, "{pad}{{}}");
        return;
    }
    for (k, v) in m {
        let key = if mode == 2 && k.parse::<u64>().is_ok() {
            k.clone()
        } else {
            yaml_key(k)
        };
        match v {
            Node::Leaf(val) => {
                let _ = writeln!(out, "{pad}{key}: {}", yaml_scalar(val));
            }
            Node::Map(inner) => {
                if inner.is_empty() {
                    let _ = writeln!(out, "{pad}{key}: {{}}");
                } else {
                    let _ = writeln!(out, "{pad}{key}:");
                    let child_mode = match mode {
                        0 if k == "dist_freq" => 1,
                        1 => 2,
                        _ => 0,
                    };
                    write_yaml(v, indent + 2, out, child_mode);
                }
            }
        }
    }
}

fn json_string(s: &str) -> String {
    serde_json::to_string(s).unwrap_or_else(|_| "\"\"".into())
}

fn write_json(node: &Node, indent: usize, out: &mut String) {
    match node {
        Node::Leaf(v) => match v {
            Value::Int(i) => out.push_str(&i.to_string()),
            Value::Float(f) => {
                if f.is_finite() {
                    out.push_str(&format_py_float(*f));
                } else {
                    out.push_str(&json_string(&format_py_float(*f)));
                }
            }
            Value::Str(s) => out.push_str(&json_string(s)),
        },
        Node::Map(m) => {
            if m.is_empty() {
                out.push_str("{}");
                return;
            }
            out.push_str("{\n");
            let pad = " ".repeat(indent + 2);
            for (i, (k, v)) in m.iter().enumerate() {
                out.push_str(&pad);
                out.push_str(&json_string(k));
                out.push_str(": ");
                write_json(v, indent + 2, out);
                if i + 1 < m.len() {
                    out.push(',');
                }
                out.push('\n');
            }
            out.push_str(&" ".repeat(indent));
            out.push('}');
        }
    }
}

/// Parse a pairtools TSV stats file into a snapshot (summary lines are
/// ignored; they are recomputed).
pub fn parse_tsv(text: &str, name: &str) -> Result<StatsSnapshot> {
    let mut flat: Vec<(String, String)> = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if line.trim().is_empty() {
            continue;
        }
        let mut it = line.splitn(2, '\t');
        let k = it.next().unwrap_or("").to_string();
        let v = it
            .next()
            .ok_or_else(|| {
                KiraError::format(
                    format!("{name} is not a valid stats file"),
                    crate::error::Location::file(name).at_line(i as u64 + 1),
                )
            })?
            .to_string();
        flat.push((k, v));
    }
    snapshot_from_flat(&flat, name)
}

/// Parse a YAML stats file (first filter only).
pub fn parse_yaml(text: &str, name: &str) -> Result<StatsSnapshot> {
    let doc: serde_yaml_ng::Value = serde_yaml_ng::from_str(text).map_err(|e| {
        KiraError::format(
            format!("{name}: invalid YAML stats: {e}"),
            crate::error::Location::file(name),
        )
    })?;
    let mapping = doc.as_mapping().ok_or_else(|| {
        KiraError::format(
            format!("{name}: YAML stats must be a mapping"),
            crate::error::Location::file(name),
        )
    })?;
    let (_, filter) = mapping.iter().next().ok_or_else(|| {
        KiraError::format(
            format!("{name}: empty YAML stats"),
            crate::error::Location::file(name),
        )
    })?;
    let mut flat = Vec::new();
    flatten_yaml(filter, String::new(), &mut flat);
    // dist_freq in YAML is strand -> bin -> count; TSV form is bin -> strand.
    let flat: Vec<(String, String)> = flat
        .into_iter()
        .map(|(k, v)| {
            let parts: Vec<&str> = k.split('/').collect();
            if parts.len() == 3 && parts[0] == "dist_freq" {
                (format!("dist_freq/{}/{}", parts[2], parts[1]), v)
            } else {
                (k, v)
            }
        })
        .collect();
    snapshot_from_flat(&flat, name)
}

fn flatten_yaml(v: &serde_yaml_ng::Value, prefix: String, out: &mut Vec<(String, String)>) {
    match v {
        serde_yaml_ng::Value::Mapping(m) => {
            for (k, v) in m {
                let ks = match k {
                    serde_yaml_ng::Value::String(s) => s.clone(),
                    serde_yaml_ng::Value::Number(n) => n.to_string(),
                    other => format!("{other:?}"),
                };
                let p = if prefix.is_empty() {
                    ks
                } else {
                    format!("{prefix}/{ks}")
                };
                flatten_yaml(v, p, out);
            }
        }
        serde_yaml_ng::Value::Number(n) => out.push((prefix, n.to_string())),
        serde_yaml_ng::Value::String(s) => out.push((prefix, s.clone())),
        serde_yaml_ng::Value::Bool(b) => out.push((prefix, b.to_string())),
        _ => {}
    }
}

fn snapshot_from_flat(flat: &[(String, String)], name: &str) -> Result<StatsSnapshot> {
    let bad = |k: &str, v: &str| {
        KiraError::format(
            format!("{name}: invalid value {v:?} for {k}"),
            crate::error::Location::file(name),
        )
    };
    let int = |k: &str, v: &str| -> Result<u64> {
        v.trim().parse::<u64>().or_else(|_| {
            v.trim()
                .parse::<f64>()
                .ok()
                .filter(|f| f.fract() == 0.0 && *f >= 0.0)
                .map(|f| f as u64)
                .ok_or_else(|| bad(k, v))
        })
    };
    // First pass: dist_freq edges.
    let mut edges: Vec<u64> = Vec::new();
    for (k, _) in flat {
        if let Some(rest) = k.strip_prefix("dist_freq/") {
            let range = rest.split('/').next().unwrap_or("");
            let lo = range.trim_end_matches('+').split('-').next().unwrap_or("");
            if let Ok(e) = lo.parse::<u64>() {
                edges.push(e);
            }
        }
    }
    let bins = if edges.is_empty() {
        DistBins::default()
    } else {
        DistBins::from_edges(edges)
    };
    let mut s = StatsSnapshot::empty(bins);
    let mut dist: BTreeMap<(usize, u64), u64> = BTreeMap::new();
    for (k, v) in flat {
        let parts: Vec<&str> = k.split('/').collect();
        match parts.as_slice() {
            ["total"] => s.total = int(k, v)?,
            ["total_unmapped"] => s.total_unmapped = int(k, v)?,
            ["total_single_sided_mapped"] => s.total_single_sided_mapped = int(k, v)?,
            ["total_mapped"] => s.total_mapped = int(k, v)?,
            ["total_dups"] => s.total_dups = int(k, v)?,
            ["total_nodups"] => s.total_nodups = int(k, v)?,
            ["cis"] => s.cis = int(k, v)?,
            ["trans"] => s.trans = int(k, v)?,
            ["pair_types", pt] => s.pair_types.push((pt.to_string(), int(k, v)?)),
            ["chrom_freq", a, b] => s
                .chrom_freq
                .push(((a.to_string(), b.to_string()), int(k, v)?)),
            ["chromsizes", c] => s.chromsizes.push((c.to_string(), int(k, v)?)),
            ["dist_freq", range, dir] => {
                let lo: u64 = range
                    .trim_end_matches('+')
                    .split('-')
                    .next()
                    .unwrap_or("")
                    .parse()
                    .map_err(|_| bad(k, v))?;
                let d = DIRS
                    .iter()
                    .position(|x| x == dir)
                    .ok_or_else(|| bad(k, v))?;
                dist.insert((d, lo), int(k, v)?);
            }
            [key] if key.starts_with("cis_") && key.ends_with("kb+") => {
                let kb: u64 = key[4..key.len() - 3].parse().map_err(|_| bad(k, v))?;
                if let Some(i) = CIS_KB.iter().position(|x| *x == kb) {
                    s.cis_kb[i] = int(k, v)?;
                }
            }
            _ => {}
        }
    }
    for ((d, lo), n) in dist {
        let i = s.bins.index(lo);
        s.dist_freq[d][i] = n;
    }
    s.pair_types.sort();
    s.chrom_freq.sort();
    Ok(s)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> StatsSnapshot {
        let mut s = StatsSnapshot::empty(DistBins::default());
        s.total = 10;
        s.total_mapped = 8;
        s.total_nodups = 6;
        s.total_dups = 2;
        s.total_unmapped = 1;
        s.total_single_sided_mapped = 1;
        s.cis = 4;
        s.trans = 2;
        s.cis_kb = [3, 2, 2, 1, 1, 0];
        s.pair_types = vec![("DD".into(), 2), ("UU".into(), 8)];
        s.chrom_freq = vec![
            (("chr1".into(), "chr1".into()), 4),
            (("chr1".into(), "chr2".into()), 2),
        ];
        s.dist_freq[0][5] = 1;
        s.dist_freq[3][40] = 3;
        s.chromsizes = vec![("chr1".into(), 1000), ("chr2".into(), 500)];
        s
    }

    #[test]
    fn tsv_roundtrip() {
        let s = sample();
        let text = render(&s, StatsFormat::Tsv);
        assert!(text.starts_with("total\t10\n"));
        assert!(text.contains("summary/frac_cis\t0.6666666666666666\n"));
        assert!(text.contains("dist_freq/0-1/+-\t0\n"));
        assert!(text.contains("dist_freq/1000000000+/++\t0\n"));
        assert!(text.contains("chromsizes/chr1\t1000\n"));
        let back = parse_tsv(&text, "t").unwrap();
        assert_eq!(back, s);
    }

    #[test]
    fn yaml_and_json_roundtrip() {
        let s = sample();
        let y = render(&s, StatsFormat::Yaml);
        assert!(y.starts_with("no_filter:\n  total: 10\n"));
        assert!(y.contains("  dist_freq:\n    '+-':\n      0: 0\n"));
        let back = parse_yaml(&y, "t").unwrap();
        assert_eq!(back, s);
        let j = render(&s, StatsFormat::Json);
        let v: serde_json::Value = serde_json::from_str(&j).unwrap();
        assert_eq!(v["no_filter"]["total"], 10);
        assert_eq!(v["no_filter"]["chrom_freq"]["chr1/chr2"], 2);
    }

    #[test]
    fn merge_sums() {
        let mut a = sample();
        let b = sample();
        a.merge(&b).unwrap();
        assert_eq!(a.total, 20);
        assert_eq!(a.pair_types, vec![("DD".into(), 4), ("UU".into(), 16)]);
        assert_eq!(a.dist_freq[3][40], 6);
        let mut c = sample();
        c.chromsizes = vec![("chrX".into(), 1)];
        assert!(a.merge(&c).is_err());
    }
}
