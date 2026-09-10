//! Integration tests: run the real CLI on committed fixtures and compare
//! against pairtools 1.1.3 golden outputs (normalising only the `@PG`
//! provenance record).

use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_kira-pairs"))
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

struct Out {
    status: i32,
    stdout: Vec<u8>,
    stderr: String,
}

fn run(args: &[&str], stdin: Option<&[u8]>) -> Out {
    let mut cmd = Command::new(bin());
    cmd.args(args).stdout(Stdio::piped()).stderr(Stdio::piped());
    cmd.stdin(if stdin.is_some() {
        Stdio::piped()
    } else {
        Stdio::null()
    });
    let mut child = cmd.spawn().expect("spawn kira-pairs");
    if let Some(data) = stdin {
        let mut si = child.stdin.take().unwrap();
        si.write_all(data).unwrap();
        drop(si);
    }
    let out = child.wait_with_output().unwrap();
    Out {
        status: out.status.code().unwrap_or(-1),
        stdout: out.stdout,
        stderr: String::from_utf8_lossy(&out.stderr).into_owned(),
    }
}

fn ok(args: &[&str]) -> Vec<u8> {
    let o = run(args, None);
    assert_eq!(o.status, 0, "command {args:?} failed: {}", o.stderr);
    o.stdout
}

/// Header lines minus provenance `@PG` records, and body lines.
fn normalise(data: &[u8]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let mut header = Vec::new();
    let mut body = Vec::new();
    for line in data.split(|b| *b == b'\n') {
        if line.is_empty() {
            continue;
        }
        if line.starts_with(b"#") {
            if line.starts_with(b"#samheader: @PG\tID:pairtools_")
                || line.starts_with(b"#samheader: @PG\tID:kira-pairs_")
            {
                continue;
            }
            // pairtools 1.1.3 `sort` slices "#chromosomes:" at a fixed offset and
            // emits a stray ':' token; kira-pairs writes the names only.
            let line = if let Some(rest) = line.strip_prefix(b"#chromosomes: : ") {
                [b"#chromosomes: ".as_slice(), rest].concat()
            } else {
                line.to_vec()
            };
            header.push(line);
        } else {
            body.push(line.to_vec());
        }
    }
    (header, body)
}

fn assert_matches_golden(actual: &[u8], golden_name: &str) {
    let expected =
        std::fs::read(golden(golden_name)).unwrap_or_else(|e| panic!("golden {golden_name}: {e}"));
    let (ha, ba) = normalise(actual);
    let (he, be) = normalise(&expected);
    assert_eq!(ba.len(), be.len(), "{golden_name}: body length differs");
    for (i, (a, e)) in ba.iter().zip(be.iter()).enumerate() {
        assert_eq!(
            String::from_utf8_lossy(a),
            String::from_utf8_lossy(e),
            "{golden_name}: body line {i} differs"
        );
    }
    assert_eq!(
        ha.iter()
            .map(|l| String::from_utf8_lossy(l).into_owned())
            .collect::<Vec<_>>(),
        he.iter()
            .map(|l| String::from_utf8_lossy(l).into_owned())
            .collect::<Vec<_>>(),
        "{golden_name}: header differs"
    );
}

fn stats_map(text: &str) -> std::collections::BTreeMap<String, String> {
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let (k, v) = l.split_once('\t').expect("stats line");
            (k.to_string(), v.to_string())
        })
        .collect()
}

fn assert_stats_match(actual: &str, golden_name: &str) {
    let expected = std::fs::read_to_string(golden(golden_name)).unwrap();
    let (a, e) = (stats_map(actual), stats_map(&expected));
    assert_eq!(
        a.keys().collect::<Vec<_>>(),
        e.keys().collect::<Vec<_>>(),
        "{golden_name}: keys"
    );
    for (k, ev) in &e {
        let av = &a[k];
        if av == ev {
            continue;
        }
        let (fa, fe): (f64, f64) = (av.parse().unwrap(), ev.parse().unwrap());
        if fa.is_nan() && fe.is_nan() {
            continue;
        }
        assert!(
            (fa - fe).abs() <= 1e-9 * fe.abs().max(1.0),
            "{golden_name}: {k}: {av} vs {ev}"
        );
    }
}

fn read_gz(path: &Path) -> Vec<u8> {
    let mut v = Vec::new();
    flate2::read::MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_end(&mut v)
        .unwrap();
    v
}

#[test]
fn help_and_version() {
    let o = run(&["--version"], None);
    assert_eq!(o.status, 0);
    assert!(String::from_utf8_lossy(&o.stdout).contains("kira-pairs"));
    for cmd in [
        "sort", "flip", "dedup", "stats", "select", "bin", "parse", "process", "generate",
    ] {
        let o = run(&[cmd, "--help"], None);
        assert_eq!(o.status, 0, "{cmd} --help");
    }
}

#[test]
fn sort_matches_pairtools() {
    let out = ok(&["sort", fixture("small.pairs").to_str().unwrap()]);
    assert_matches_golden(&out, "small.sorted.pairs");
    let out = ok(&["sort", fixture("edge.pairs").to_str().unwrap()]);
    assert_matches_golden(&out, "edge.sorted.pairs");
    // Unknown header fields and comments round-trip.
    let text = String::from_utf8_lossy(&out);
    assert!(text.contains("#custom_field: round trip me"));
    assert!(text.contains("# a free-form comment line"));
    assert!(text.contains("#sorted: chr1-chr2-pos1-pos2"));
}

#[test]
fn sort_external_runs_and_thread_counts_are_deterministic() {
    let dir = tempfile::tempdir().unwrap();
    let data = dir.path().join("gen.pairs");
    ok(&[
        "generate",
        "--records",
        "400000",
        "--seed",
        "3",
        "--chromosomes",
        "8",
        "-o",
        data.to_str().unwrap(),
    ]);
    let reference = ok(&["sort", data.to_str().unwrap(), "--threads", "1"]);
    // Tiny memory budget with a tiny fan-in forces several runs and an
    // intermediate merge pass.
    for threads in ["2", "8"] {
        let out = ok(&[
            "sort",
            data.to_str().unwrap(),
            "--threads",
            threads,
            "--memory",
            "64M",
            "--max-fan-in",
            "2",
            "--tmpdir",
            dir.path().to_str().unwrap(),
        ]);
        assert_eq!(out, reference, "threads={threads}");
    }
    // Sorted order and stability: equal keys keep input order.
    let (_, body) = normalise(&reference);
    assert_eq!(body.len(), 400_000);
    // No temp files left behind.
    let leftovers: Vec<_> = std::fs::read_dir(dir.path())
        .unwrap()
        .filter(|e| e.as_ref().unwrap().file_name() != "gen.pairs")
        .collect();
    assert!(leftovers.is_empty(), "temporary files left: {leftovers:?}");
}

#[test]
fn sort_ties_across_runs_keep_input_order() {
    // 50k records with only 4 distinct keys -> massive ties spanning runs.
    let mut data = b"## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n".to_vec();
    for i in 0..50_000u32 {
        let k = i % 4;
        data.extend_from_slice(
            format!(
                "r{i:06}\tchr{}\t{}\tchr1\t10\t+\t-\tUU\n",
                1 + k % 2,
                100 + k / 2
            )
            .as_bytes(),
        );
    }
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("ties.pairs");
    std::fs::write(&p, &data).unwrap();
    let out = ok(&[
        "sort",
        p.to_str().unwrap(),
        "--memory",
        "64M",
        "--threads",
        "4",
    ]);
    let (_, body) = normalise(&out);
    let mut last_id_per_key: std::collections::HashMap<(String, String), u32> =
        std::collections::HashMap::new();
    for line in &body {
        let f: Vec<&str> = std::str::from_utf8(line).unwrap().split('\t').collect();
        let id: u32 = f[0][1..].parse().unwrap();
        let key = (f[1].to_string(), f[2].to_string());
        if let Some(prev) = last_id_per_key.get(&key) {
            assert!(
                *prev < id,
                "input order lost for key {key:?}: {prev} then {id}"
            );
        }
        last_id_per_key.insert(key, id);
    }
}

#[test]
fn stdin_stdout_pipeline_and_gzip() {
    let input = std::fs::read(fixture("small.pairs")).unwrap();
    let sorted = run(&["sort"], Some(&input));
    assert_eq!(sorted.status, 0, "{}", sorted.stderr);
    assert_matches_golden(&sorted.stdout, "small.sorted.pairs");
    let dedup = run(&["dedup"], Some(&sorted.stdout));
    assert_eq!(dedup.status, 0, "{}", dedup.stderr);
    assert_matches_golden(&dedup.stdout, "small.dedup.pairs");
    let stats = run(&["stats"], Some(&sorted.stdout));
    assert_eq!(stats.status, 0);
    assert_stats_match(&String::from_utf8_lossy(&stats.stdout), "small.stats");
    // gzip input (BGZF) by magic, gzip output by extension.
    let gz_in = ok(&["sort", fixture("small.pairs.gz").to_str().unwrap()]);
    assert_matches_golden(&gz_in, "small.sorted.pairs");
    let dir = tempfile::tempdir().unwrap();
    let outp = dir.path().join("out.pairs.gz");
    ok(&[
        "sort",
        fixture("small.pairs").to_str().unwrap(),
        "-o",
        outp.to_str().unwrap(),
    ]);
    assert_matches_golden(&read_gz(&outp), "small.sorted.pairs");
    let lz4p = dir.path().join("out.pairs.lz4");
    ok(&[
        "sort",
        fixture("small.pairs").to_str().unwrap(),
        "-o",
        lz4p.to_str().unwrap(),
    ]);
    let back = ok(&["sort", lz4p.to_str().unwrap()]);
    assert_matches_golden(&back, "small.sorted.pairs");
}

#[test]
fn dedup_matches_pairtools() {
    let dir = tempfile::tempdir().unwrap();
    let dups = dir.path().join("dups.pairs");
    let stats = dir.path().join("stats.txt");
    let out = ok(&[
        "dedup",
        golden("small.sorted.pairs").to_str().unwrap(),
        "--output-dups",
        dups.to_str().unwrap(),
        "--output-stats",
        stats.to_str().unwrap(),
    ]);
    assert_matches_golden(&out, "small.dedup.pairs");
    assert_matches_golden(&std::fs::read(&dups).unwrap(), "small.dups.pairs");
    assert_stats_match(
        &std::fs::read_to_string(&stats).unwrap(),
        "small.dedup.stats",
    );
    let out = ok(&[
        "dedup",
        golden("small.sorted.pairs").to_str().unwrap(),
        "--method",
        "sum",
        "--max-mismatch",
        "5",
        "--output-dups",
        dups.to_str().unwrap(),
    ]);
    assert_matches_golden(&out, "small.dedup_sum5.pairs");
    assert_matches_golden(&std::fs::read(&dups).unwrap(), "small.dups_sum5.pairs");
    let out = ok(&[
        "dedup",
        golden("small.sorted.pairs").to_str().unwrap(),
        "--keep-parent-id",
        "--output-dups",
        "-",
    ]);
    assert_matches_golden(&out, "small.dedup_with_dups_parent.pairs");
}

#[test]
fn dedup_adversarial_cases_match_pairtools() {
    let cases = fixture("dedup_cases.pairs");
    // Greedy (cython) goldens are compared without parent ids: pairtools'
    // cython backend reports wrong parent_readID values (see docs).
    for (extra, name) in [
        (
            vec!["--method", "max", "--keep-parent-id"],
            "dedup_cases.max.pairs",
        ),
        (
            vec!["--method", "sum", "--keep-parent-id"],
            "dedup_cases.sum.pairs",
        ),
        (
            vec!["--method", "max", "--backend", "cython"],
            "dedup_cases.max.cython.pairs",
        ),
        (
            vec!["--method", "sum", "--clustering", "greedy"],
            "dedup_cases.sum.cython.pairs",
        ),
    ] {
        let mut args = vec!["dedup", cases.to_str().unwrap(), "--output-dups", "-"];
        args.extend(extra);
        let out = ok(&args);
        assert_matches_golden(&out, name);
    }
    let dir = tempfile::tempdir().unwrap();
    let stats = dir.path().join("s.txt");
    let out = ok(&[
        "dedup",
        cases.to_str().unwrap(),
        "--output-unmapped",
        "-",
        "--output-dups",
        "-",
        "--output-stats",
        stats.to_str().unwrap(),
    ]);
    assert_matches_golden(&out, "dedup_cases.all.pairs");
    assert_stats_match(
        &std::fs::read_to_string(&stats).unwrap(),
        "dedup_cases.stats",
    );
}

#[test]
fn dedup_groups_crossing_buffer_and_chunk_boundaries() {
    // Duplicate groups straddling 4 MiB parser blocks and the dedup pending
    // window: every 5000th record is duplicated by the next one.
    let mut data = b"## pairs format v1.0.0\n#sorted: chr1-chr2-pos1-pos2\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n".to_vec();
    let mut expected_dups = 0u64;
    let n = 300_000u64;
    let mut pos = 1_000u64;
    let mut i = 0u64;
    while i < n {
        let padding = "x".repeat(40);
        data.extend_from_slice(
            format!(
                "r{i}_{padding}\tchr1\t{pos}\tchr1\t{}\t+\t-\tUU\n",
                pos + 5000
            )
            .as_bytes(),
        );
        i += 1;
        if i.is_multiple_of(97) {
            // near duplicate within 3 bp on both sides
            data.extend_from_slice(
                format!(
                    "d{i}_{padding}\tchr1\t{}\tchr1\t{}\t+\t-\tUU\n",
                    pos + 2,
                    pos + 5003
                )
                .as_bytes(),
            );
            expected_dups += 1;
            i += 1;
        }
        pos += 20;
    }
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("big.pairs");
    std::fs::write(&p, &data).unwrap();
    let dups = dir.path().join("dups.pairs");
    for threads in ["1", "3"] {
        let out = ok(&[
            "dedup",
            p.to_str().unwrap(),
            "--threads",
            threads,
            "--output-dups",
            dups.to_str().unwrap(),
        ]);
        let (_, body) = normalise(&out);
        let (_, dbody) = normalise(&std::fs::read(&dups).unwrap());
        assert_eq!(dbody.len() as u64, expected_dups, "threads={threads}");
        assert!(dbody.iter().all(|l| l.starts_with(b"d")));
        assert_eq!(body.len() as u64 + expected_dups, i);
    }
}

#[test]
fn dedup_rejects_unsorted_input() {
    let data = b"## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\nr1\tchr1\t100\tchr1\t200\t+\t-\tUU\nr2\tchr1\t50\tchr1\t200\t+\t-\tUU\n";
    let o = run(&["dedup"], Some(data));
    assert_ne!(o.status, 0);
    assert!(o.stderr.contains("not sorted"), "{}", o.stderr);
    assert!(o.stderr.contains("line 4"), "{}", o.stderr);
}

#[test]
fn stats_matches_pairtools_and_merges() {
    let out = ok(&["stats", golden("small.sorted.pairs").to_str().unwrap()]);
    assert_stats_match(&String::from_utf8_lossy(&out), "small.stats");
    let yaml = ok(&[
        "stats",
        "--yaml",
        golden("small.sorted.pairs").to_str().unwrap(),
    ]);
    let parsed: serde_yaml_ng::Value = serde_yaml_ng::from_slice(&yaml).unwrap();
    let expected: serde_yaml_ng::Value =
        serde_yaml_ng::from_str(&std::fs::read_to_string(golden("small.stats.yaml")).unwrap())
            .unwrap();
    assert_eq!(parsed["no_filter"]["total"], expected["no_filter"]["total"]);
    assert_eq!(
        parsed["no_filter"]["chrom_freq"],
        expected["no_filter"]["chrom_freq"]
    );
    assert_eq!(
        parsed["no_filter"]["dist_freq"],
        expected["no_filter"]["dist_freq"]
    );
    assert_eq!(
        parsed["no_filter"]["pair_types"],
        expected["no_filter"]["pair_types"]
    );
    let json = ok(&[
        "stats",
        "--json",
        golden("small.sorted.pairs").to_str().unwrap(),
    ]);
    let j: serde_json::Value = serde_json::from_slice(&json).unwrap();
    assert_eq!(
        j["no_filter"]["total"],
        expected["no_filter"]["total"].as_u64().unwrap()
    );
    // merge: two copies double every counter
    let merged = ok(&[
        "stats",
        "--merge",
        golden("small.stats").to_str().unwrap(),
        golden("small.stats").to_str().unwrap(),
    ]);
    let m = stats_map(&String::from_utf8_lossy(&merged));
    let single = stats_map(&std::fs::read_to_string(golden("small.stats")).unwrap());
    assert_eq!(
        m["total"].parse::<u64>().unwrap(),
        2 * single["total"].parse::<u64>().unwrap()
    );
    assert_eq!(
        m["cis"].parse::<u64>().unwrap(),
        2 * single["cis"].parse::<u64>().unwrap()
    );
    assert_eq!(m["summary/frac_cis"], single["summary/frac_cis"]);
}

#[test]
fn flip_matches_pairtools() {
    let out = ok(&[
        "flip",
        fixture("small.unflipped.pairs").to_str().unwrap(),
        "-c",
        fixture("small.chrom.sizes").to_str().unwrap(),
    ]);
    assert_matches_golden(&out, "small.flipped.pairs");
}

#[test]
fn select_matches_pairtools() {
    let dir = tempfile::tempdir().unwrap();
    let rest = dir.path().join("rest.pairs");
    let out = ok(&[
        "select",
        "(pair_type == \"UU\") and (abs(pos1-pos2) < 1000)",
        fixture("small.pairs").to_str().unwrap(),
        "--output-rest",
        rest.to_str().unwrap(),
    ]);
    assert_matches_golden(&out, "small.select.pairs");
    assert_matches_golden(&std::fs::read(&rest).unwrap(), "small.select_rest.pairs");
    let o = run(
        &[
            "select",
            "nosuchcol == 1",
            fixture("small.pairs").to_str().unwrap(),
        ],
        None,
    );
    assert_ne!(o.status, 0);
    assert!(o.stderr.contains("unknown column"));
}

#[test]
fn bin_conserves_counts_and_supports_formats() {
    let dir = tempfile::tempdir().unwrap();
    let coo = dir.path().join("c.tsv");
    let bg2 = dir.path().join("c.bg2");
    let bins = dir.path().join("bins.tsv");
    let sorted = golden("small.dedup.pairs");
    let cs = fixture("small.chrom.sizes");
    ok(&[
        "bin",
        sorted.to_str().unwrap(),
        "-c",
        cs.to_str().unwrap(),
        "-r",
        "100000",
        "-o",
        coo.to_str().unwrap(),
        "--bins-out",
        bins.to_str().unwrap(),
    ]);
    ok(&[
        "bin",
        sorted.to_str().unwrap(),
        "-c",
        cs.to_str().unwrap(),
        "-r",
        "100000",
        "-o",
        bg2.to_str().unwrap(),
        "--format",
        "bg2",
    ]);
    let coo_text = std::fs::read_to_string(&coo).unwrap();
    let total: u64 = coo_text
        .lines()
        .map(|l| l.split('\t').nth(2).unwrap().parse::<u64>().unwrap())
        .sum();
    let mapped = std::fs::read_to_string(&sorted)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .filter(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            f[1] != "!" && f[3] != "!"
        })
        .count() as u64;
    assert_eq!(total, mapped);
    let mut prev = (0u64, 0u64);
    for l in coo_text.lines() {
        let f: Vec<u64> = l.split('\t').map(|x| x.parse().unwrap()).collect();
        assert!(f[0] <= f[1]);
        assert!((f[0], f[1]) > prev || prev == (0, 0));
        prev = (f[0], f[1]);
    }
    let bg2_total: u64 = std::fs::read_to_string(&bg2)
        .unwrap()
        .lines()
        .map(|l| l.split('\t').nth(6).unwrap().parse::<u64>().unwrap())
        .sum();
    assert_eq!(bg2_total, mapped);
    let n_bins = std::fs::read_to_string(&bins).unwrap().lines().count();
    let expected_bins: u64 = std::fs::read_to_string(&cs)
        .unwrap()
        .lines()
        .map(|l| {
            l.split('\t')
                .nth(1)
                .unwrap()
                .parse::<u64>()
                .unwrap()
                .div_ceil(100_000)
        })
        .sum();
    assert_eq!(n_bins as u64, expected_bins);
    // multi-resolution in one pass
    let tmpl = dir.path().join("multi.{res}.tsv");
    ok(&[
        "bin",
        sorted.to_str().unwrap(),
        "-c",
        cs.to_str().unwrap(),
        "-r",
        "50000",
        "-r",
        "500000",
        "-o",
        tmpl.to_str().unwrap(),
    ]);
    for r in ["50000", "500000"] {
        let t = std::fs::read_to_string(dir.path().join(format!("multi.{r}.tsv"))).unwrap();
        let s: u64 = t
            .lines()
            .map(|l| l.split('\t').nth(2).unwrap().parse::<u64>().unwrap())
            .sum();
        assert_eq!(s, mapped, "resolution {r}");
    }
}

#[test]
fn parse_matches_pairtools() {
    let bam = fixture("hic.bam");
    let cs = fixture("hic.chrom.sizes");
    let out = ok(&["parse", "-c", cs.to_str().unwrap(), bam.to_str().unwrap()]);
    assert_matches_golden(&out, "hic.parse.pairsam");
    let out = ok(&[
        "parse",
        "-c",
        cs.to_str().unwrap(),
        bam.to_str().unwrap(),
        "--drop-sam",
    ]);
    assert_matches_golden(&out, "hic.parse.pairs");
    let out = ok(&[
        "parse",
        "-c",
        cs.to_str().unwrap(),
        fixture("hic.sam").to_str().unwrap(),
        "--drop-sam",
    ]);
    let (_, body) = normalise(&out);
    let (_, gbody) = normalise(&std::fs::read(golden("hic.parse.pairs")).unwrap());
    assert_eq!(body, gbody);
    let out = ok(&[
        "parse",
        "-c",
        cs.to_str().unwrap(),
        bam.to_str().unwrap(),
        "--drop-sam",
        "--add-columns",
        "mapq,pos5,pos3,cigar,read_len,matched_bp,algn_ref_span,algn_read_span,dist_to_5,dist_to_3,read_side,algn_idx,same_side_algn_count,NM,AS,XS",
        "--add-pair-index",
    ]);
    assert_matches_golden(&out, "hic.parse.cols.pairs");
    let out = ok(&[
        "parse",
        "-c",
        cs.to_str().unwrap(),
        bam.to_str().unwrap(),
        "--drop-sam",
        "--walks-policy",
        "mask",
        "--min-mapq",
        "30",
    ]);
    assert_matches_golden(&out, "hic.parse.mask_mapq30.pairs");
    let dir = tempfile::tempdir().unwrap();
    let stats = dir.path().join("s.txt");
    ok(&[
        "parse",
        "-c",
        cs.to_str().unwrap(),
        bam.to_str().unwrap(),
        "--drop-sam",
        "--output-stats",
        stats.to_str().unwrap(),
        "-o",
        "/dev/null",
    ]);
    assert_stats_match(&std::fs::read_to_string(&stats).unwrap(), "hic.parse.stats");
    let o = run(
        &[
            "parse",
            "-c",
            cs.to_str().unwrap(),
            bam.to_str().unwrap(),
            "--walks-policy",
            "all",
        ],
        None,
    );
    assert_ne!(o.status, 0);
}

#[test]
fn process_matches_parse_sort_dedup() {
    let bam = fixture("hic.bam");
    let cs = fixture("hic.chrom.sizes");
    let dir = tempfile::tempdir().unwrap();
    let pairs = dir.path().join("p.pairs.gz");
    let dups = dir.path().join("d.pairs");
    let stats = dir.path().join("s.txt");
    let bins = dir.path().join("b.tsv");
    for threads in ["1", "6"] {
        ok(&[
            "process",
            bam.to_str().unwrap(),
            "-c",
            cs.to_str().unwrap(),
            "--drop-sam",
            "--threads",
            threads,
            "--output-pairs",
            pairs.to_str().unwrap(),
            "--output-dups",
            dups.to_str().unwrap(),
            "--output-stats",
            stats.to_str().unwrap(),
            "--output-bins",
            bins.to_str().unwrap(),
            "--resolution",
            "100000",
        ]);
        assert_matches_golden(&read_gz(&pairs), "hic.process.pairs");
        assert_matches_golden(&std::fs::read(&dups).unwrap(), "hic.process.dups.pairs");
        assert_stats_match(
            &std::fs::read_to_string(&stats).unwrap(),
            "hic.process.stats",
        );
        let total: u64 = std::fs::read_to_string(&bins)
            .unwrap()
            .lines()
            .map(|l| l.split('\t').nth(2).unwrap().parse::<u64>().unwrap())
            .sum();
        assert!(total > 0);
    }
}

#[test]
fn malformed_input_gives_clear_errors() {
    let bad_pos = b"## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\nr1\tchr1\t12x\tchr1\t200\t+\t-\tUU\n";
    let o = run(&["sort"], Some(bad_pos));
    assert_ne!(o.status, 0);
    assert!(
        o.stderr.contains("line 3") && o.stderr.contains("column 3") && o.stderr.contains("12x"),
        "{}",
        o.stderr
    );
    let bad_strand = b"## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\nr1\tchr1\t12\tchr1\t200\t*\t-\tUU\n";
    let o = run(&["stats"], Some(bad_strand));
    assert_ne!(o.status, 0);
    assert!(o.stderr.contains("strand"), "{}", o.stderr);
    let few = b"## pairs format v1.0.0\n#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\nr1\tchr1\t12\n";
    let o = run(&["dedup"], Some(few));
    assert_ne!(o.status, 0);
    assert!(o.stderr.contains("fields"), "{}", o.stderr);
    let o = run(&["sort", "/nonexistent/file.pairs"], None);
    assert_ne!(o.status, 0);
    assert!(o.stderr.contains("nonexistent"));
    let dir = tempfile::tempdir().unwrap();
    let trunc = dir.path().join("t.pairs.gz");
    let full = std::fs::read(fixture("small.pairs.gz")).unwrap();
    std::fs::write(&trunc, &full[..full.len() / 2]).unwrap();
    let o = run(&["sort", trunc.to_str().unwrap()], None);
    assert_ne!(o.status, 0, "truncated gzip must fail");
    let bad_cs = dir.path().join("bad.sizes");
    std::fs::write(&bad_cs, "chr1\tabc\n").unwrap();
    let o = run(
        &[
            "flip",
            fixture("small.pairs").to_str().unwrap(),
            "-c",
            bad_cs.to_str().unwrap(),
        ],
        None,
    );
    assert_ne!(o.status, 0);
    assert!(o.stderr.contains("chromosome sizes"), "{}", o.stderr);
}

#[test]
fn stdout_is_data_only_and_metrics_go_to_stderr() {
    let input = std::fs::read(fixture("small.pairs")).unwrap();
    let o = run(&["sort", "--metrics", "-v"], Some(&input));
    assert_eq!(o.status, 0);
    assert!(o.stderr.contains("records_read"));
    assert!(!String::from_utf8_lossy(&o.stdout).contains("records_read"));
    assert_matches_golden(&o.stdout, "small.sorted.pairs");
}

#[test]
fn memory_budget_is_respected() {
    // 64 MiB budget on a ~100 MB input must not balloon memory: check the
    // peak RSS reported by the process itself stays well under 8x budget.
    let dir = tempfile::tempdir().unwrap();
    let data = dir.path().join("gen.pairs");
    ok(&[
        "generate",
        "--records",
        "1500000",
        "--seed",
        "5",
        "-o",
        data.to_str().unwrap(),
    ]);
    let o = run(
        &[
            "sort",
            data.to_str().unwrap(),
            "--memory",
            "64M",
            "--threads",
            "4",
            "--metrics",
            "-o",
            "/dev/null",
            "--tmpdir",
            dir.path().to_str().unwrap(),
        ],
        None,
    );
    assert_eq!(o.status, 0, "{}", o.stderr);
    let rss_line = o
        .stderr
        .lines()
        .find(|l| l.starts_with("peak_rss_bytes"))
        .expect("peak_rss metric");
    let rss: u64 = rss_line.split('\t').nth(1).unwrap().parse().unwrap();
    let runs: u64 = o
        .stderr
        .lines()
        .find(|l| l.starts_with("number_of_runs"))
        .unwrap()
        .split('\t')
        .nth(1)
        .unwrap()
        .parse()
        .unwrap();
    assert!(runs >= 2, "expected external runs, got {runs}");
    assert!(
        rss < 512 * 1024 * 1024,
        "peak RSS {} exceeds 512 MiB for a 64M budget",
        rss
    );
}
