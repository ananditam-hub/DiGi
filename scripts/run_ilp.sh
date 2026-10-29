#!/bin/bash

# Infer a sample's two MHC haplotypes with the recombination-capped ILP
# (src/ilp/dip_main_recomb.cpp) from its read k-mers and the graph node k-mers.
# Usage: bash run_ilp.sh <graph_file> <node_kmers.tsv> <read_kmers.tsv> <sample> <out_dir>
#   graph_file and node_kmers.tsv come from build_graph.sh (debug_<graph_file> must sit
#   next to graph_file); read_kmers.tsv comes from process_reads.sh.
#   C1 (missing k-mer penalty, default 5.2), C2 (novel k-mer penalty, default 3.5) and
#   R (max recombinations per haplotype, default 10) can be set as environment variables.
#   Needs bin/dip_ilp (make ilp) and a Gurobi license (e.g. GRB_LICENSE_FILE or WLS variables).
# -> <out_dir>/<sample>_{1,2}_<weights>_<read_kmers name>.fa (+ _haplotype_regions.tsv), <out_dir>/<sample>_<weights>.log

set -euo pipefail

LABELS="$(realpath "$1")"
GK="$(realpath "$2")"
SK="$(realpath "$3")"
SAMPLE="$4"
OUTDIR="$5"
c1=${C1:-5.2}
c2=${C2:-3.5}
R=${R:-10}

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/bin/dip_ilp"
[ -x "$BIN" ] || { echo "missing $BIN, run make ilp in $ROOT first" >&2; exit 1; }

COORDS="$(dirname "$LABELS")/debug_$(basename "$LABELS")"

for f in "$LABELS" "$COORDS" "$GK" "$SK"; do
    [ -s "$f" ] || { echo "missing $f" >&2; exit 1; }
done

# H = haplotypes in the graph, r = junctions per haplotype = nodes per haplotype - 1
H=$(grep -c '^Processing ' "$COORDS")
nodes=$(wc -l < "$LABELS")
[ $((nodes % H)) -eq 0 ] || { echo "$nodes nodes not divisible by H=$H" >&2; exit 1; }
r=$((nodes / H - 1))

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

rel_weights="c1-${c1}_c2-${c2}_R-${R}"
cov=$(basename "$SK" .tsv)
echo "Running ILP on $SAMPLE (H=$H, r=$r, $rel_weights)..."

# Run from the output directory so Gurobi's dip_ilp.log lands there
cd "$OUTDIR"
"$BIN" --graph_kmers "$GK" \
  --sk "$SK" \
  --c1 "$c1" --c2 "$c2" \
  --R "$R" \
  --H "$H" --r "$r" \
  --labels "$LABELS" \
  --node_coords "$COORDS" \
  --sample "$SAMPLE" \
  --out_dir "$OUTDIR" \
  --cov "$cov" \
  --rel_weights "$rel_weights" \
  > "$OUTDIR/${SAMPLE}_${rel_weights}.log" 2>&1 \
  || { echo "ILP failed, see $OUTDIR/${SAMPLE}_${rel_weights}.log" >&2; exit 1; }

echo "Inferred haplotypes in $OUTDIR"
