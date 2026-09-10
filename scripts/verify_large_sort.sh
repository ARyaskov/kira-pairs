#!/usr/bin/env bash
# Verify memory-bounded sorting on data much larger than the budget:
# generates N records, sorts with a small --memory and a small --max-fan-in
# (forcing several runs and intermediate merge passes), checks ordering and
# record conservation with a streaming awk check, and reports peak RSS.
#
# Usage: verify_large_sort.sh [--records N] [--memory SIZE] [--threads N] [--tmpdir DIR]
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
KP="${KIRA_PAIRS:-$ROOT/target/release/kira-pairs}"
RECORDS=30000000; MEMORY=256M; THREADS=$(nproc); TMP="${TMPDIR:-/tmp}"
while [ $# -gt 0 ]; do case "$1" in
  --records) RECORDS=$2; shift 2;; --memory) MEMORY=$2; shift 2;; --threads) THREADS=$2; shift 2;; --tmpdir) TMP=$2; shift 2;;
  *) echo "unknown option $1"; exit 2;; esac; done
WORK=$(mktemp -d -p "$TMP" kira-verify-XXXX); trap 'rm -rf "$WORK"' EXIT
echo "# generating $RECORDS records"
$KP generate --records "$RECORDS" --seed 9 --chromosomes 30 -o "$WORK/in.pairs" --threads "$THREADS"
ls -la "$WORK/in.pairs"
echo "# sorting with --memory $MEMORY --max-fan-in 4 --threads $THREADS"
/usr/bin/time -f "wall=%es user=%Us sys=%Ss maxrss=%MkB" \
  $KP sort "$WORK/in.pairs" -o "$WORK/out.pairs" --memory "$MEMORY" --max-fan-in 4 --threads "$THREADS" --tmpdir "$WORK" --metrics 2>&1 | grep -E "wall=|number_of_runs|merge_passes|temporary_bytes|peak_rss"
echo "# verifying order and conservation"
grep -vc '^#' "$WORK/in.pairs"
LC_ALL=C awk -F'\t' '!/^#/ { n++; key=$2 SUBSEP $4; if (n>1) { if ($2<pc1 || ($2==pc1 && ($4<pc2 || ($4==pc2 && ($3+0<pp1 || ($3+0==pp1 && ($5+0<pp2 || ($5+0==pp2 && $8<ppt)))))))) { print "ORDER VIOLATION at record " n; exit 1 } } pc1=$2; pc2=$4; pp1=$3+0; pp2=$5+0; ppt=$8 } END { print n " records in sorted order" }' "$WORK/out.pairs"
cmp <(grep -v '^#' "$WORK/in.pairs" | LC_ALL=C sort) <(grep -v '^#' "$WORK/out.pairs" | LC_ALL=C sort) && echo "# record multiset conserved"
echo "# leftover temp files: $(ls "$WORK" | grep -vc 'pairs$' || true)"
