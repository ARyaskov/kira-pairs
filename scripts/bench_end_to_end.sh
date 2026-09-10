#!/usr/bin/env bash
# End-to-end benchmarks of kira-pairs against pairtools 1.1.3.
#
# Measures wall time, CPU time, peak RSS and throughput for sort, dedup,
# stats, sort|dedup and the full parse|sort|dedup pipeline on synthetic
# datasets that are large enough not to fit in CPU caches. Results are
# appended as TSV to $OUT (default bench-results.tsv).
#
# Usage: bench_end_to_end.sh [--records N] [--threads N] [--memory SIZE]
#        [--pairtools PATH] [--kira PATH] [--tmpdir DIR] [--out FILE] [--bam FILE --chroms FILE]
set -euo pipefail
export LC_ALL=C
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RECORDS=10000000
THREADS=$(nproc)
MEMORY=4G
PT="${PAIRTOOLS:-pairtools}"
KP="${KIRA_PAIRS:-$ROOT/target/release/kira-pairs}"
TMPDIR_ARG="${TMPDIR:-/tmp}"
OUT="$ROOT/bench-results.tsv"
BAM=""
CHROMS=""
while [ $# -gt 0 ]; do
  case "$1" in
    --records) RECORDS=$2; shift 2;;
    --threads) THREADS=$2; shift 2;;
    --memory) MEMORY=$2; shift 2;;
    --pairtools) PT=$2; shift 2;;
    --kira) KP=$2; shift 2;;
    --tmpdir) TMPDIR_ARG=$2; shift 2;;
    --out) OUT=$2; shift 2;;
    --bam) BAM=$2; shift 2;;
    --chroms) CHROMS=$2; shift 2;;
    *) echo "unknown option $1"; exit 2;;
  esac
done
WORK=$(mktemp -d -p "$TMPDIR_ARG" kira-bench-XXXX)
trap 'rm -rf "$WORK"' EXIT
TIME=/usr/bin/time
[ -x $TIME ] || { echo "/usr/bin/time is required"; exit 1; }

measure() { # label cmd...
  local label=$1; shift
  local tf="$WORK/time.$$"
  sync; { echo 3 > /proc/sys/vm/drop_caches; } 2>/dev/null || true
  $TIME -f "%e\t%U\t%S\t%M" -o "$tf" bash -c "$*"
  read -r wall user sys rss < "$tf"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$(date -Iseconds)" "$label" "$RECORDS" "$wall" "$user" "$sys" "$rss" "$THREADS" >> "$OUT"
  awk -v l="$label" -v w="$wall" -v u="$user" -v s="$sys" -v r="$rss" 'BEGIN { printf "%-50s wall %8.2fs  cpu %8.2fs  peakRSS %8.1f MB\n", l, w, u + s, r / 1024 }'
}

[ -f "$OUT" ] || printf "timestamp\tlabel\trecords\twall_s\tuser_s\tsys_s\tmax_rss_kb\tthreads\n" > "$OUT"
echo "# host: $(uname -srm); cpu: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs); mem: $(free -g | awk '/Mem/ {print $2}') GB; rust: $(rustc --version); pairtools: $($PT --version 2>/dev/null | tail -1)"
echo "# records=$RECORDS threads=$THREADS memory=$MEMORY tmpdir=$WORK"

echo "# generating dataset"
$KP generate --records "$RECORDS" --seed 42 --duplicate-rate 0.15 --extra-columns 0 -o "$WORK/bench.pairs" --threads "$THREADS"
$KP generate --records "$RECORDS" --seed 42 --duplicate-rate 0.15 --extra-columns 0 -o "$WORK/bench.pairs.gz" --threads "$THREADS"
ls -la "$WORK"/bench.pairs*

measure "kira sort (plain in, plain out)" "$KP sort $WORK/bench.pairs -o $WORK/k.sorted.pairs --threads $THREADS --memory $MEMORY --tmpdir $WORK"
measure "kira sort (gz in, gz out)" "$KP sort $WORK/bench.pairs.gz -o $WORK/k.sorted.pairs.gz --threads $THREADS --memory $MEMORY --tmpdir $WORK"
measure "pairtools sort (plain in, plain out)" "$PT sort $WORK/bench.pairs -o $WORK/pt.sorted.pairs --nproc $THREADS --memory $MEMORY --tmpdir $WORK"
measure "pairtools sort (gz in, gz out)" "$PT sort $WORK/bench.pairs.gz -o $WORK/pt.sorted.pairs.gz --nproc $THREADS --memory $MEMORY --tmpdir $WORK"
cmp <(grep -v '^#' "$WORK/k.sorted.pairs") <(grep -v '^#' "$WORK/pt.sorted.pairs") && echo "# sort outputs identical"

measure "kira dedup (plain)" "$KP dedup $WORK/k.sorted.pairs -o $WORK/k.dedup.pairs --output-stats $WORK/k.dedup.stats --threads $THREADS"
measure "pairtools dedup (plain)" "$PT dedup $WORK/pt.sorted.pairs -o $WORK/pt.dedup.pairs --output-stats $WORK/pt.dedup.stats"
cmp <(grep -v '^#' "$WORK/k.dedup.pairs") <(grep -v '^#' "$WORK/pt.dedup.pairs") && echo "# dedup outputs identical"

measure "kira stats" "$KP stats $WORK/k.sorted.pairs -o $WORK/k.stats --threads $THREADS"
measure "pairtools stats" "$PT stats $WORK/pt.sorted.pairs -o $WORK/pt.stats"

measure "kira sort | dedup (gz in, gz out)" "$KP sort $WORK/bench.pairs.gz --threads $THREADS --memory $MEMORY --tmpdir $WORK | $KP dedup -o $WORK/k.chain.pairs.gz --output-stats $WORK/k.chain.stats --threads $THREADS"
measure "pairtools sort | dedup (gz in, gz out)" "$PT sort $WORK/bench.pairs.gz --nproc $THREADS --memory $MEMORY --tmpdir $WORK | $PT dedup -o $WORK/pt.chain.pairs.gz --output-stats $WORK/pt.chain.stats"

if [ -n "$BAM" ] && [ -n "$CHROMS" ]; then
  measure "kira process (bam -> pairs.gz + stats + bins)" "$KP process $BAM -c $CHROMS --drop-sam --output-pairs $WORK/k.proc.pairs.gz --output-stats $WORK/k.proc.stats --output-bins $WORK/k.proc.bins.gz --resolution 10000 --threads $THREADS --memory $MEMORY --tmpdir $WORK"
  measure "pairtools parse | sort | dedup (bam -> pairs.gz + stats)" "$PT parse -c $CHROMS $BAM --drop-sam --nproc-in $THREADS | $PT sort --nproc $THREADS --memory $MEMORY --tmpdir $WORK | $PT dedup -o $WORK/pt.proc.pairs.gz --output-stats $WORK/pt.proc.stats"
  measure "kira parse (bam -> pairs.gz)" "$KP parse $BAM -c $CHROMS --drop-sam -o $WORK/k.parse.pairs.gz --threads $THREADS"
  measure "pairtools parse (bam -> pairs.gz)" "$PT parse -c $CHROMS $BAM --drop-sam -o $WORK/pt.parse.pairs.gz --nproc-in $THREADS"
fi
echo "# results appended to $OUT"
