#!/usr/bin/env python3
"""Generate deterministic SAM/BAM fixtures covering the alignment cases that
`kira-pairs parse` must handle identically to `pairtools parse`.

Usage: generate_fixture.py OUT_PREFIX [--seed N] [--reads N]

Writes OUT_PREFIX.sam, OUT_PREFIX.bam (if pysam is importable) and
OUT_PREFIX.chrom.sizes. No aligner is needed: records are synthesised with
realistic flags, CIGARs, MAPQs, SA tags and supplementary/secondary records.
"""
import argparse
import random
import sys

CHROMS = [("chr1", 5_000_000), ("chr2", 3_000_000), ("chr10", 2_000_000), ("chrX", 1_000_000), ("chrM", 16_569)]
READ_LEN = 100
BASES = "ACGT"


def seq(rng, n):
    return "".join(rng.choice(BASES) for _ in range(n))


def qual(n):
    return "I" * n


def sam_line(qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, s, q, tags):
    return "\t".join([qname, str(flag), rname, str(pos), str(mapq), cigar, rnext, str(pnext), str(tlen), s, q] + tags)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--reads", type=int, default=400)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    lines = ["@HD\tVN:1.6\tSO:queryname"]
    for c, l in CHROMS:
        lines.append(f"@SQ\tSN:{c}\tLN:{l}")
    lines.append("@RG\tID:rg1\tSM:sample1")
    lines.append("@PG\tID:bwa\tPN:bwa\tVN:0.7.17-r1188\tCL:bwa mem -SP5M ref.fa r1.fq r2.fq")
    lines.append("@CO\tsynthetic fixture for kira-pairs")
    records = []

    def add_pair(i, kind):
        name = f"read{i:05d}:{kind}"
        s1, s2 = seq(rng, READ_LEN), seq(rng, READ_LEN)
        q = qual(READ_LEN)
        c1, l1 = rng.choice(CHROMS)
        p1 = rng.randint(1, l1 - 2 * READ_LEN)
        cis = rng.random() < 0.7
        if cis:
            c2, l2 = c1, l1
            d = int(10 ** rng.uniform(2, 6))
            p2 = min(p1 + d, l2 - READ_LEN)
        else:
            c2, l2 = rng.choice(CHROMS)
            p2 = rng.randint(1, l2 - READ_LEN)
        f1 = 0x1 | 0x40 | (0x10 if rng.random() < 0.5 else 0)
        f2 = 0x1 | 0x80 | (0x10 if rng.random() < 0.5 else 0)
        mq = lambda: rng.choice([0, 1, 5, 20, 30, 40, 60])
        tags1 = [f"NM:i:{rng.randint(0,3)}", f"AS:i:{rng.randint(50,100)}", "XS:i:0", "RG:Z:rg1"]
        tags2 = [f"NM:i:{rng.randint(0,3)}", f"AS:i:{rng.randint(50,100)}", "XS:i:0", "RG:Z:rg1"]
        if kind == "UU":
            recs = [
                sam_line(name, f1, c1, p1, 60, f"{READ_LEN}M", "=" if c1 == c2 else c2, p2, 0, s1, q, tags1),
                sam_line(name, f2, c2, p2, 60, f"{READ_LEN}M", "=" if c1 == c2 else c1, p1, 0, s2, q, tags2),
            ]
        elif kind == "MQ":  # varying MAPQ on both sides
            recs = [
                sam_line(name, f1, c1, p1, mq(), f"{READ_LEN}M", "=", p2, 0, s1, q, tags1),
                sam_line(name, f2, c2, p2, mq(), f"{READ_LEN}M", "=", p1, 0, s2, q, tags2),
            ]
        elif kind == "SOFT":  # soft clips of various lengths (gap conversion)
            cl1 = rng.choice([5, 10, 19, 20, 21, 30, 45])
            cl2 = rng.choice([0, 10, 25])
            cig1 = f"{cl1}S{READ_LEN-cl1}M" if not (f1 & 0x10) else f"{READ_LEN-cl1}M{cl1}S"
            cig2 = f"{READ_LEN-cl2}M{cl2}S" if cl2 else f"{READ_LEN}M"
            recs = [
                sam_line(name, f1, c1, p1, 60, cig1, "=", p2, 0, s1, q, tags1),
                sam_line(name, f2, c2, p2, 60, cig2, "=", p1, 0, s2, q, tags2),
            ]
        elif kind == "UN":  # read 2 unmapped
            recs = [
                sam_line(name, f1 | 0x8, c1, p1, 60, f"{READ_LEN}M", "=", p1, 0, s1, q, tags1),
                sam_line(name, 0x1 | 0x80 | 0x4, c1, p1, 0, "*", "=", p1, 0, s2, q, ["RG:Z:rg1"]),
            ]
        elif kind == "NU":  # read 1 unmapped
            recs = [
                sam_line(name, 0x1 | 0x40 | 0x4, c2, p2, 0, "*", "=", p2, 0, s1, q, ["RG:Z:rg1"]),
                sam_line(name, f2 | 0x8, c2, p2, 60, f"{READ_LEN}M", "=", p2, 0, s2, q, tags2),
            ]
        elif kind == "NN":
            recs = [
                sam_line(name, 0x1 | 0x40 | 0x4 | 0x8, "*", 0, 0, "*", "*", 0, 0, s1, q, ["RG:Z:rg1"]),
                sam_line(name, 0x1 | 0x80 | 0x4 | 0x8, "*", 0, 0, "*", "*", 0, 0, s2, q, ["RG:Z:rg1"]),
            ]
        elif kind == "MM":  # multimappers (mapq 0) with secondary records
            recs = [
                sam_line(name, f1, c1, p1, 0, f"{READ_LEN}M", "=", p2, 0, s1, q, tags1 + ["XA:Z:chr2,+100,100M,0;"]),
                sam_line(name, f1 | 0x100, "chr2", 100, 0, f"{READ_LEN}M", "=", p2, 0, "*", "*", ["NM:i:0"]),
                sam_line(name, f2, c2, p2, 60, f"{READ_LEN}M", "=", p1, 0, s2, q, tags2),
            ]
        elif kind == "CHIM":  # simple chimeric read 1 with a supplementary record (SA tags)
            split = rng.choice([40, 50, 60])
            # 5' part: first `split` bases map at p1; 3' part maps elsewhere.
            strand1 = "-" if f1 & 0x10 else "+"
            c3, l3 = rng.choice(CHROMS)
            p3 = rng.randint(1, l3 - READ_LEN)
            rescue = rng.random() < 0.5
            if rescue:
                # 3' part on the mate's chromosome, pointing towards it, within 750 bp.
                c3 = c2
                if f2 & 0x10:
                    p3 = max(1, p2 - 200)
                    strand3 = "+"
                else:
                    p3 = p2 + 200
                    strand3 = "-"
            else:
                strand3 = rng.choice("+-")
            f3 = 0x1 | 0x40 | 0x800 | (0x10 if strand3 == "-" else 0)
            if strand1 == "+":
                cig_primary = f"{split}M{READ_LEN-split}S"
            else:
                cig_primary = f"{READ_LEN-split}S{split}M"
            if strand3 == "+":
                cig_supp = f"{split}H{READ_LEN-split}M"
            else:
                cig_supp = f"{READ_LEN-split}M{split}H"
            sa_primary = f"SA:Z:{c3},{p3},{strand3},{cig_supp.replace('H','S')},60,0;"
            sa_supp = f"SA:Z:{c1},{p1},{strand1},{cig_primary},60,0;"
            recs = [
                sam_line(name, f1, c1, p1, 60, cig_primary, "=", p2, 0, s1, q, tags1 + [sa_primary]),
                sam_line(name, f3, c3, p3, 60, cig_supp, "=", p2, 0, s1[split:] if strand3 == "+" else s1[: READ_LEN - split], q[split:], ["NM:i:0", sa_supp]),
                sam_line(name, f2, c2, p2, 60, f"{READ_LEN}M", "=", p1, 0, s2, q, tags2),
            ]
        elif kind == "CHIM3":  # three-part chimera on read 1 (never rescuable)
            recs = [
                sam_line(name, f1 & ~0x10, c1, p1, 60, "30M70S", "=", p2, 0, s1, q, tags1 + ["SA:Z:x"]),
                sam_line(name, 0x1 | 0x40 | 0x800, "chr2", 1000 + i, rng.choice([0, 60]), "30H40M30H", "=", p2, 0, "*", "*", ["SA:Z:x"]),
                sam_line(name, 0x1 | 0x40 | 0x800 | 0x10, "chr10", 2000 + i, 60, "30M70H", "=", p2, 0, "*", "*", ["SA:Z:x"]),
                sam_line(name, f2, c2, p2, 60, f"{READ_LEN}M", "=", p1, 0, s2, q, tags2),
            ]
        elif kind == "DUP":  # exact and near duplicates of a fixed location
            base1, base2 = 100_000 + (i % 7) * 3, 200_000 + (i % 5) * 2
            recs = [
                sam_line(name, 0x1 | 0x40 | 0x20, "chr1", base1, 60, f"{READ_LEN}M", "=", base2, 0, s1, q, tags1),
                sam_line(name, 0x1 | 0x80 | 0x10, "chr1", base2, 60, f"{READ_LEN}M", "=", base1, 0, s2, q, tags2),
            ]
        else:
            raise ValueError(kind)
        # Random record order within the read (BWA writes R1 first; shuffle sometimes).
        if rng.random() < 0.2:
            recs.reverse()
        records.extend(recs)

    kinds = ["UU"] * 10 + ["MQ"] * 3 + ["SOFT"] * 3 + ["UN", "NU", "NN", "MM"] + ["CHIM"] * 4 + ["CHIM3"] + ["DUP"] * 4
    for i in range(a.reads):
        add_pair(i, rng.choice(kinds))
    lines.extend(records)
    with open(a.prefix + ".sam", "w") as f:
        f.write("\n".join(lines) + "\n")
    with open(a.prefix + ".chrom.sizes", "w") as f:
        # Deliberately partial order with a chromosome missing (chrM) and one
        # absent from the SAM (chrY) to exercise pairtools' ordering rules.
        for c, l in [("chr1", 5_000_000), ("chr2", 3_000_000), ("chrX", 1_000_000), ("chrY", 500_000), ("chr10", 2_000_000)]:
            f.write(f"{c}\t{l}\n")
    try:
        import pysam
        pysam.view("-b", "-o", a.prefix + ".bam", a.prefix + ".sam", catch_stdout=False)
        print("wrote", a.prefix + ".bam")
    except ImportError:
        print("pysam not available; BAM not written", file=sys.stderr)
    print("wrote", a.prefix + ".sam", "with", len(records), "records")


if __name__ == "__main__":
    main()
