#!/usr/bin/env bash
# Regenerate committed fixtures (tests/fixtures) and pairtools 1.1.3 golden
# outputs (tests/golden). Requires `pairtools` 1.1.3 on PATH (or PAIRTOOLS=...)
# with pysam, and a release build of kira-pairs (for the synthetic generator).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PT="${PAIRTOOLS:-pairtools}"
PY="${PYTHON:-python3}"
KP="${KIRA_PAIRS:-$ROOT/target/release/kira-pairs}"
FX="$ROOT/tests/fixtures"
GD="$ROOT/tests/golden"
mkdir -p "$FX" "$GD"
cd "$FX"

echo "# synthetic pairs (unsorted, whole matrix orientation already upper-triangular for mapped pairs)"
"$KP" generate --records 3000 --seed 11 --chromosomes 12 --chrom-length 5000000 --extra-columns 1 --duplicate-rate 0.2 --unmapped-fraction 0.05 -o small.pairs
grep '^#chromsize:' small.pairs | awk '{print $2"\t"$3}' > small.chrom.sizes
"$PY" - <<'PY'
import random
random.seed(3)
with open('small.unflipped.pairs','w') as out:
    for l in open('small.pairs'):
        if l.startswith('#'):
            out.write(l); continue
        f=l.rstrip('\n').split('\t')
        if random.random()<0.5:
            f[1],f[3]=f[3],f[1]; f[2],f[4]=f[4],f[2]; f[5],f[6]=f[6],f[5]; f[7]=f[7][1]+f[7][0]
        out.write('\t'.join(f)+'\n')
PY
bgzip -c small.pairs > small.pairs.gz

echo "# hand-written edge cases"
cat > edge.pairs <<'PAIRS'
## pairs format v1.0.0
#shape: upper triangle
#genome_assembly: test
#chromsize: chr1 300000000
#chromsize: chr2 250000000
#chromsize: chr10 130000000
#chromosomes: chr2 chr10 chr1
#custom_field: round trip me
#samheader: @SQ	SN:chr1	LN:300000000
#samheader: @PG	ID:bwa	PN:bwa	VN:0.7.17
# a free-form comment line
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2 note
r1	chr10	5	chr10	9	+	+	UU	60	60	small
r2	chr1	10	chr1	20	+	-	UU	60	60	tie1
r3	chr1	10	chr1	20	+	-	UU	60	60	tie2
r4	chr1	10	chr1	20	+	-	UR	60	60	pairtype-order
r5	chr1	10	chr1	20	+	-	UU	60	60	tie3
r6	chr2	4294967297	chr2	18446744073709551615	-	-	UU	1	2	huge
r7	!	0	chr1	100	-	+	NU	0	30	single-sided
r8	!	0	!	0	-	-	NN	0	0	unmapped
r9	chr1	9	chr2	1	+	+	UU	10	20	trans
r10	chr1	100	chr1	100	-	-	UU	60	60	self
r11	chr1	100	chr10	5	+	-	UU	3	4	chr1-vs-chr10
r12	chr1	2	chr1	1000000	+	+	UU	60	60	x
PAIRS

echo "# adversarial dedup cases (already block-sorted and flipped)"
cat > dedup_cases.pairs <<'PAIRS'
## pairs format v1.0.0
#sorted: chr1-chr2-pos1-pos2
#shape: upper triangle
#chromsize: chr1 1000000
#chromsize: chr2 1000000
#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type
u1	!	0	!	0	-	-	NN
u2	!	0	chr1	500	-	+	NU
a	chr1	100	chr1	100	+	+	UU
b	chr1	103	chr1	100	+	+	UU
c	chr1	106	chr1	100	+	+	UU
d	chr1	109	chr1	100	+	+	UU
e	chr1	113	chr1	100	+	+	UU
s1	chr1	200	chr1	200	+	-	UU
s2	chr1	200	chr1	203	+	-	UU
s3	chr1	200	chr1	206	+	-	UU
s4	chr1	201	chr1	209	+	-	UU
t1	chr1	300	chr1	300	-	-	UU
t2	chr1	302	chr1	302	-	-	UU
t3	chr1	303	chr1	305	-	-	UU
t4	chr1	306	chr1	301	-	-	UU
m1	chr1	400	chr1	400	+	+	UU
m2	chr1	400	chr1	406	+	+	UU
m3	chr1	401	chr1	403	+	+	UU
m4	chr1	900	chr1	100	+	+	UU
x1	chr1	1000	chr1	1000	+	+	UU
x2	chr1	1000	chr1	1000	+	-	UU
x3	chr1	1000	chr1	1000	-	+	UU
x4	chr1	1000	chr1	1000	-	-	UU
x5	chr1	1000	chr1	1000	+	+	UR
y1	chr1	5000	chr2	5000	+	+	UU
y2	chr1	5003	chr2	5003	+	+	UU
y3	chr1	5003	chr2	5007	+	+	UU
z1	chr2	10	chr2	20	+	+	UU
z2	chr2	10	chr2	20	+	+	UU
z3	chr2	10	chr2	20	+	+	UU
PAIRS

echo "# SAM/BAM fixture"
"$PY" "$ROOT/scripts/generate_fixture.py" hic --reads 400 --seed 1 > /dev/null

echo "# golden outputs from $($PT --version 2>/dev/null | tail -1)"
strip_pg() { grep -v $'^#samheader: @PG\tID:pairtools_'; }
$PT sort small.pairs 2>/dev/null | strip_pg > "$GD/small.sorted.pairs"
$PT sort edge.pairs 2>/dev/null | strip_pg > "$GD/edge.sorted.pairs"
$PT dedup "$GD/small.sorted.pairs" --output-dups "$GD/small.dups.pairs" --output-stats "$GD/small.dedup.stats" 2>/dev/null | strip_pg > "$GD/small.dedup.pairs"
$PT dedup "$GD/small.sorted.pairs" --keep-parent-id --output-dups - 2>/dev/null | strip_pg > "$GD/small.dedup_with_dups_parent.pairs"
$PT dedup "$GD/small.sorted.pairs" --method sum --max-mismatch 5 --output-dups "$GD/small.dups_sum5.pairs" 2>/dev/null | strip_pg > "$GD/small.dedup_sum5.pairs"
$PT stats "$GD/small.sorted.pairs" 2>/dev/null > "$GD/small.stats"
$PT stats --yaml "$GD/small.sorted.pairs" 2>/dev/null > "$GD/small.stats.yaml"
$PT flip small.unflipped.pairs -c small.chrom.sizes 2>/dev/null | strip_pg > "$GD/small.flipped.pairs"
$PT select '(pair_type == "UU") and (abs(pos1-pos2) < 1000)' small.pairs --output-rest "$GD/small.select_rest.pairs" 2>/dev/null | strip_pg > "$GD/small.select.pairs"
sed -i $'/^#samheader: @PG\tID:pairtools_/d' "$GD/small.select_rest.pairs" "$GD/small.dups.pairs" "$GD/small.dups_sum5.pairs" 2>/dev/null || true
for m in max sum; do
  $PT dedup dedup_cases.pairs --method $m --keep-parent-id --output-dups - 2>/dev/null | strip_pg > "$GD/dedup_cases.$m.pairs"
  # The cython backend's parent_readID values are wrong after its internal
  # buffer shrinks (indices are not rebased), so its golden covers
  # classification only.
  $PT dedup dedup_cases.pairs --method $m --backend cython --output-dups - 2>/dev/null | strip_pg > "$GD/dedup_cases.$m.cython.pairs"
done
$PT dedup dedup_cases.pairs --output-unmapped - --output-dups - --output-stats "$GD/dedup_cases.stats" 2>/dev/null | strip_pg > "$GD/dedup_cases.all.pairs"
$PT parse -c hic.chrom.sizes hic.bam 2>/dev/null | strip_pg > "$GD/hic.parse.pairsam"
$PT parse -c hic.chrom.sizes hic.bam --drop-sam 2>/dev/null | strip_pg > "$GD/hic.parse.pairs"
$PT parse -c hic.chrom.sizes hic.bam --drop-sam --add-columns mapq,pos5,pos3,cigar,read_len,matched_bp,algn_ref_span,algn_read_span,dist_to_5,dist_to_3,read_side,algn_idx,same_side_algn_count,NM,AS,XS --add-pair-index 2>/dev/null | strip_pg > "$GD/hic.parse.cols.pairs"
$PT parse -c hic.chrom.sizes hic.bam --drop-sam --walks-policy mask --min-mapq 30 2>/dev/null | strip_pg > "$GD/hic.parse.mask_mapq30.pairs"
$PT parse -c hic.chrom.sizes hic.bam --drop-sam --output-stats "$GD/hic.parse.stats" -o /dev/null 2>/dev/null
$PT parse -c hic.chrom.sizes hic.bam --drop-sam 2>/dev/null | $PT sort 2>/dev/null | $PT dedup --output-dups "$GD/hic.process.dups.pairs" --output-stats "$GD/hic.process.stats" 2>/dev/null | strip_pg > "$GD/hic.process.pairs"
sed -i $'/^#samheader: @PG\tID:pairtools_/d' "$GD/hic.process.dups.pairs"
ls -la "$FX" "$GD"
