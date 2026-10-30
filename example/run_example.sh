#!/bin/bash

# Infer HG002's two MHC haplotypes from 10x MHC-filtered Illumina reads with an
# 18-haplotype HGSVC graph cut every 10 kb.
# Usage: bash example/run_example.sh
#   Needs the digi conda env and `make` in the repo root; the ILP step also needs
#   `make ilp` and a Gurobi license, and is skipped without bin/dip_ilp.
#   FRACTION, K, THREADS, C1, C2 and R are passed through to the pipeline scripts.
# -> example/work/{panel,liftover,ref_graph,reads,inferred}

set -euo pipefail

EX="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(dirname "$EX")"
cd "$EX"

export FRACTION=${FRACTION:-0.05}
export K=${K:-31}
export C1=${C1:-6.5127} C2=${C2:-4.0225} R=${R:-10}

mkdir -p work
[ -d work/panel ] || tar -C work -xJf data/HGSVC_MHC_panel_18hap_grid10kb.tar.xz

# 1. Graph and node k-mers
bash "$ROOT/scripts/build_graph.sh" graph.yml

# 2. Read k-mers
bash "$ROOT/scripts/process_reads.sh" data/HG002.MHC_R1.10x.fastq.gz data/HG002.MHC_R2.10x.fastq.gz HG002 work/reads

# 3. Diploid inference
if [ -x "$ROOT/bin/dip_ilp" ]; then
  bash "$ROOT/scripts/run_ilp.sh" work/ref_graph/ref_graph.txt \
    "work/ref_graph/node_kmers_${FRACTION}_k${K}.tsv" \
    "work/reads/HG002_read_kmers_${FRACTION}_k${K}.tsv" \
    HG002 work/inferred
else
  echo "bin/dip_ilp not built (make ilp, needs Gurobi), skipping the ILP step" >&2
fi
