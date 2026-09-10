#!/usr/bin/env python3
"""Differential compatibility tests: run the same inputs through pairtools and
kira-pairs and compare normalised outputs.

Only provenance fields are normalised (the `@PG` record each tool appends to
`#samheader:`); everything else, including body bytes, must match exactly.
Stats files are compared as key/value maps with a relative tolerance for
floating point summaries.

Usage:
  compare_pairtools.py --kira PATH --pairtools PATH [--fixtures DIR] [--quick] [--big N]
"""
import argparse
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def run(cmd, stdin=None, stdout=None, check=True):
    r = subprocess.run(cmd, stdin=stdin, stdout=stdout or subprocess.PIPE, stderr=subprocess.PIPE, shell=isinstance(cmd, str))
    if check and r.returncode != 0:
        raise RuntimeError(f"command failed ({r.returncode}): {cmd}\n{r.stderr.decode(errors='replace')[-2000:]}")
    return r


def read_pairs(path):
    header, body = [], []
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b"#"):
                if line.startswith(b"#samheader: @PG\tID:pairtools_") or line.startswith(b"#samheader: @PG\tID:kira-pairs_"):
                    continue
                # pairtools 1.1.3 sort emits a stray ':' token in #chromosomes:
                if line.startswith(b"#chromosomes: : "):
                    line = b"#chromosomes: " + line[len(b"#chromosomes: : "):]
                header.append(line)
            else:
                body.append(line)
    return header, body


def read_stats(path):
    out = {}
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            k, v = line.rstrip("\n").split("\t", 1)
            out[k] = v
    return out


def stats_equal(a, b):
    if set(a) != set(b):
        return False, f"keys differ: {sorted(set(a) ^ set(b))[:10]}"
    for k in a:
        x, y = a[k], b[k]
        if x == y:
            continue
        try:
            fx, fy = float(x), float(y)
        except ValueError:
            return False, f"{k}: {x!r} != {y!r}"
        if math.isnan(fx) and math.isnan(fy):
            continue
        if not math.isclose(fx, fy, rel_tol=1e-9, abs_tol=1e-12):
            return False, f"{k}: {x} != {y}"
    return True, ""


class Suite:
    def __init__(self, kira, pairtools, tmp):
        self.kira, self.pt, self.tmp = kira, pairtools, tmp
        self.passed, self.failed = 0, 0

    def check_pairs(self, name, pt_path, k_path, compare_header=True):
        ha, ba = read_pairs(pt_path)
        hb, bb = read_pairs(k_path)
        ok = ba == bb and (not compare_header or ha == hb)
        self.report(name, ok, "" if ok else self.first_diff(ha, hb, ba, bb))

    def check_stats(self, name, pt_path, k_path):
        ok, msg = stats_equal(read_stats(pt_path), read_stats(k_path))
        self.report(name, ok, msg)

    def report(self, name, ok, msg):
        if ok:
            self.passed += 1
            print(f"PASS {name}")
        else:
            self.failed += 1
            print(f"FAIL {name}: {msg}")

    @staticmethod
    def first_diff(ha, hb, ba, bb):
        if ha != hb:
            for i, (x, y) in enumerate(zip(ha, hb)):
                if x != y:
                    return f"header line {i}: {x!r} vs {y!r}"
            return f"header length {len(ha)} vs {len(hb)}"
        for i, (x, y) in enumerate(zip(ba, bb)):
            if x != y:
                return f"body line {i}: {x!r} vs {y!r}"
        return f"body length {len(ba)} vs {len(bb)}"

    def path(self, name):
        return os.path.join(self.tmp, name)


def make_unflipped(src, dst, seed=3):
    rng = random.Random(seed)
    with open(src) as f, open(dst, "w") as out:
        for line in f:
            if line.startswith("#"):
                out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if rng.random() < 0.5:
                fields[1], fields[3] = fields[3], fields[1]
                fields[2], fields[4] = fields[4], fields[2]
                fields[5], fields[6] = fields[6], fields[5]
                fields[7] = fields[7][1] + fields[7][0]
            out.write("\t".join(fields) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kira", default=os.path.join(ROOT, "target", "release", "kira-pairs"))
    ap.add_argument("--pairtools", default="pairtools")
    ap.add_argument("--fixtures", default=os.path.join(ROOT, "tests", "fixtures"))
    ap.add_argument("--quick", action="store_true", help="skip the larger synthetic dataset")
    ap.add_argument("--big", type=int, default=300000, help="records in the synthetic dataset")
    ap.add_argument("--keep", action="store_true")
    a = ap.parse_args()
    tmp = tempfile.mkdtemp(prefix="kira-compat-")
    s = Suite(a.kira, a.pairtools, tmp)
    fx = a.fixtures
    try:
        datasets = [("small", os.path.join(fx, "small.pairs"), os.path.join(fx, "small.chrom.sizes"))]
        if not a.quick:
            gen = s.path("gen.pairs")
            run([a.kira, "generate", "--records", str(a.big), "--seed", "7", "--chromosomes", "30", "--extra-columns", "1", "-o", gen])
            cs = s.path("gen.chrom.sizes")
            with open(gen) as f, open(cs, "w") as out:
                for line in f:
                    if line.startswith("#chromsize:"):
                        _, c, l = line.split()
                        out.write(f"{c}\t{l}\n")
                    elif not line.startswith("#"):
                        break
            datasets.append(("gen", gen, cs))
        for label, pairs, chromsizes in datasets:
            pt_sorted, k_sorted = s.path(f"{label}.pt.sorted.pairs"), s.path(f"{label}.k.sorted.pairs")
            run([a.pairtools, "sort", pairs, "-o", pt_sorted])
            run([a.kira, "sort", pairs, "-o", k_sorted, "--memory", "64M", "--threads", "4"])
            s.check_pairs(f"{label}: sort", pt_sorted, k_sorted)
            # sort via stdin/stdout and gzip
            with open(pairs, "rb") as fin, open(s.path(f"{label}.k.sorted2.pairs"), "wb") as fout:
                run([a.kira, "sort"], stdin=fin, stdout=fout)
            s.check_pairs(f"{label}: sort stdin->stdout", pt_sorted, s.path(f"{label}.k.sorted2.pairs"))
            gz = s.path(f"{label}.k.sorted.pairs.gz")
            run([a.kira, "sort", pairs, "-o", gz])
            run(f"gzip -dc {gz} > {s.path(label + '.k.sorted3.pairs')}")
            s.check_pairs(f"{label}: sort gz output", pt_sorted, s.path(f"{label}.k.sorted3.pairs"))
            # dedup
            for extra, tag in [([], "max3"), (["--method", "sum", "--max-mismatch", "5"], "sum5"), (["--max-mismatch", "0"], "max0"), (["--backend", "cython"], "cython")]:
                pt_out, pt_dups, pt_stats = s.path(f"{label}.{tag}.pt.dedup.pairs"), s.path(f"{label}.{tag}.pt.dups.pairs"), s.path(f"{label}.{tag}.pt.stats")
                k_out, k_dups, k_stats = s.path(f"{label}.{tag}.k.dedup.pairs"), s.path(f"{label}.{tag}.k.dups.pairs"), s.path(f"{label}.{tag}.k.stats")
                run([a.pairtools, "dedup", pt_sorted, "-o", pt_out, "--output-dups", pt_dups, "--output-stats", pt_stats] + extra)
                run([a.kira, "dedup", k_sorted, "-o", k_out, "--output-dups", k_dups, "--output-stats", k_stats] + extra)
                s.check_pairs(f"{label}: dedup {tag}", pt_out, k_out)
                s.check_pairs(f"{label}: dedup {tag} dups", pt_dups, k_dups)
                s.check_stats(f"{label}: dedup {tag} stats", pt_stats, k_stats)
            # dedup with dups in the same output + no-mark-dups
            pt_out, k_out = s.path(f"{label}.pt.dedup_all.pairs"), s.path(f"{label}.k.dedup_all.pairs")
            run([a.pairtools, "dedup", pt_sorted, "-o", pt_out, "--output-dups", "-", "--no-mark-dups"])
            run([a.kira, "dedup", k_sorted, "-o", k_out, "--output-dups", "-", "--no-mark-dups"])
            s.check_pairs(f"{label}: dedup dups-in-output no-mark", pt_out, k_out)
            # stats
            pt_stats, k_stats = s.path(f"{label}.pt.stats"), s.path(f"{label}.k.stats")
            run([a.pairtools, "stats", pt_sorted, "-o", pt_stats])
            run([a.kira, "stats", k_sorted, "-o", k_stats])
            s.check_stats(f"{label}: stats", pt_stats, k_stats)
            # stats merge
            run([a.pairtools, "stats", "--merge", pt_stats, pt_stats, "-o", s.path(f"{label}.pt.merged")])
            run([a.kira, "stats", "--merge", k_stats, k_stats, "-o", s.path(f"{label}.k.merged")])
            s.check_stats(f"{label}: stats --merge", s.path(f"{label}.pt.merged"), s.path(f"{label}.k.merged"))
            # flip
            unflipped = s.path(f"{label}.unflipped.pairs")
            make_unflipped(pairs, unflipped)
            run([a.pairtools, "flip", unflipped, "-c", chromsizes, "-o", s.path(f"{label}.pt.flip.pairs")])
            run([a.kira, "flip", unflipped, "-c", chromsizes, "-o", s.path(f"{label}.k.flip.pairs")])
            s.check_pairs(f"{label}: flip", s.path(f"{label}.pt.flip.pairs"), s.path(f"{label}.k.flip.pairs"))
            # select
            for q in ['(pair_type == "UU") and (abs(pos1-pos2) < 1000)', "chrom1==chrom2", 'regex_match(chrom1, "chr1\\d") and (chrom2 != "!")', 'csv_match(pair_type, "UU,UR") or wildcard_match(chrom2, "chr1*")']:
                run([a.pairtools, "select", q, pairs, "-o", s.path("pt.sel.pairs"), "--output-rest", s.path("pt.rest.pairs")])
                run([a.kira, "select", q, pairs, "-o", s.path("k.sel.pairs"), "--output-rest", s.path("k.rest.pairs")])
                s.check_pairs(f"{label}: select {q}", s.path("pt.sel.pairs"), s.path("k.sel.pairs"))
                s.check_pairs(f"{label}: select rest {q}", s.path("pt.rest.pairs"), s.path("k.rest.pairs"))
        # dedup adversarial cases
        cases = os.path.join(fx, "dedup_cases.pairs")
        # cython backend: parent ids are not compared (pairtools bug, see docs).
        for extra, tag in [(["--method", "max", "--keep-parent-id"], "max"), (["--method", "sum", "--keep-parent-id"], "sum"), (["--backend", "cython"], "cython"), (["--backend", "cython", "--method", "sum"], "cython-sum")]:
            run([a.pairtools, "dedup", cases, "-o", s.path("pt.cases.pairs"), "--output-dups", "-", "--output-unmapped", "-"] + extra)
            run([a.kira, "dedup", cases, "-o", s.path("k.cases.pairs"), "--output-dups", "-", "--output-unmapped", "-"] + extra)
            s.check_pairs(f"dedup adversarial {tag}", s.path("pt.cases.pairs"), s.path("k.cases.pairs"))
        # parse
        bam, cs = os.path.join(fx, "hic.bam"), os.path.join(fx, "hic.chrom.sizes")
        variants = [[], ["--drop-sam"], ["--drop-sam", "--add-columns", "mapq,pos5,pos3,cigar,read_len,matched_bp,algn_ref_span,algn_read_span,dist_to_5,dist_to_3,read_side,algn_idx,same_side_algn_count,NM,AS,XS", "--add-pair-index"],
                    ["--drop-sam", "--walks-policy", "mask"], ["--drop-sam", "--walks-policy", "5any"], ["--drop-sam", "--walks-policy", "3any"], ["--drop-sam", "--walks-policy", "3unique"],
                    ["--drop-sam", "--min-mapq", "30"], ["--drop-sam", "--min-mapq", "0"], ["--drop-sam", "--no-flip"], ["--drop-sam", "--report-alignment-end", "3"], ["--drop-sam", "--max-inter-align-gap", "5"], ["--drop-sam", "--max-molecule-size", "100"], ["--drop-readid", "--drop-sam"]]
        for v in variants:
            run([a.pairtools, "parse", "-c", cs, bam, "-o", s.path("pt.parse.pairs")] + v)
            run([a.kira, "parse", "-c", cs, bam, "-o", s.path("k.parse.pairs")] + v)
            s.check_pairs(f"parse {' '.join(v) or '(defaults)'}", s.path("pt.parse.pairs"), s.path("k.parse.pairs"))
        run([a.pairtools, "parse", "-c", cs, os.path.join(fx, "hic.sam"), "-o", s.path("pt.parse.pairs"), "--drop-sam"])
        run([a.kira, "parse", "-c", cs, os.path.join(fx, "hic.sam"), "-o", s.path("k.parse.pairs"), "--drop-sam"])
        s.check_pairs("parse SAM input", s.path("pt.parse.pairs"), s.path("k.parse.pairs"))
        run([a.pairtools, "parse", "-c", cs, bam, "--drop-sam", "--output-stats", s.path("pt.parse.stats"), "-o", os.devnull])
        run([a.kira, "parse", "-c", cs, bam, "--drop-sam", "--output-stats", s.path("k.parse.stats"), "-o", os.devnull])
        s.check_stats("parse --output-stats", s.path("pt.parse.stats"), s.path("k.parse.stats"))
        # fused process vs pairtools chain
        run(f"{a.pairtools} parse -c {cs} {bam} --drop-sam | {a.pairtools} sort | {a.pairtools} dedup -o {s.path('pt.proc.pairs')} --output-dups {s.path('pt.proc.dups')} --output-stats {s.path('pt.proc.stats')}")
        run([a.kira, "process", bam, "-c", cs, "--drop-sam", "--output-pairs", s.path("k.proc.pairs"), "--output-dups", s.path("k.proc.dups"), "--output-stats", s.path("k.proc.stats"), "--output-bins", s.path("k.proc.bins"), "--resolution", "100000"])
        s.check_pairs("process vs parse|sort|dedup", s.path("pt.proc.pairs"), s.path("k.proc.pairs"))
        s.check_pairs("process dups", s.path("pt.proc.dups"), s.path("k.proc.dups"))
        s.check_stats("process stats", s.path("pt.proc.stats"), s.path("k.proc.stats"))
        # interop chains
        run(f"{a.kira} parse -c {cs} {bam} --drop-sam | {a.kira} sort | {a.pairtools} dedup -o {s.path('x1.pairs')}")
        s.check_pairs("kira parse|sort -> pairtools dedup", s.path("pt.proc.pairs"), s.path("x1.pairs"))
        run(f"{a.pairtools} parse -c {cs} {bam} --drop-sam | {a.kira} sort | {a.kira} dedup -o {s.path('x2.pairs')}")
        s.check_pairs("pairtools parse -> kira sort|dedup", s.path("pt.proc.pairs"), s.path("x2.pairs"))
    finally:
        if a.keep:
            print("kept", tmp)
        else:
            shutil.rmtree(tmp, ignore_errors=True)
    print(f"\n{s.passed} passed, {s.failed} failed")
    sys.exit(1 if s.failed else 0)


if __name__ == "__main__":
    main()
