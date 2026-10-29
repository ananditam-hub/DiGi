#!/bin/bash

# Build the junction-linked haplotype graph (node sequences) and its
# FracMinHash-sampled node k-mers.
# Usage: bash build_graph.sh <graph.yml> [samples_included.txt]
#   FRACTION (default 0.05) and K (default 31) can be set as environment variables.
# -> <graph_file>, debug_<graph_file>, node_kmers_<fraction>_k<k>.tsv (next to graph_file)

set -euo pipefail

CONFIG="$1"
SAMPLES="${2:-}"
fraction=${FRACTION:-0.05}
k=${K:-31}

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/bin/node_kmers_gen"
[ -x "$BIN" ] || { echo "missing $BIN, run make in $ROOT first" >&2; exit 1; }

# Graph output path, read from the config
GRAPH_FILE=$(python3 -c "import sys, yaml; print(yaml.safe_load(open(sys.argv[1]))['graph']['graph_file'])" "$CONFIG")
OUTDIR=$(dirname "$GRAPH_FILE")
DEBUG="$OUTDIR/debug_$(basename "$GRAPH_FILE")"
KMERS="$OUTDIR/node_kmers_${fraction}_k${k}.tsv"
mkdir -p "$OUTDIR"

# Build graph nodes
python3 "$ROOT/src/graph/graph_builder.py" --config "$CONFIG" ${SAMPLES:+--samples_included "$SAMPLES"} > "$DEBUG"

# Every haplotype must have the same number of nodes (junctions + 1)
n_haps=$(grep -c '^Processing ' "$DEBUG" || true)
grep -o 'wrote [0-9]* nodes' "$DEBUG" | sort | uniq -c
[ "$(grep -o 'wrote [0-9]* nodes' "$DEBUG" | sort -u | wc -l)" = 1 ] \
  && [ "$(grep -c 'wrote [0-9]* nodes' "$DEBUG")" = "$n_haps" ] \
  || { echo "unequal node counts / missing haplotypes, see $DEBUG" >&2; exit 1; }

# Sample node k-mers
"$BIN" --fraction "$fraction" --k "$k" --input "$GRAPH_FILE" --output "$KMERS.tmp" \
  && [ -s "$KMERS.tmp" ] && mv "$KMERS.tmp" "$KMERS" \
  || { echo "node k-mer sampling failed" >&2; exit 1; }

echo "Graph:      $GRAPH_FILE"
echo "Node k-mers: $KMERS"
