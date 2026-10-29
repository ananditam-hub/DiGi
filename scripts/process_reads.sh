#!/bin/bash

# Quality-filter MHC-filtered paired-end reads and count their
# FracMinHash-sampled k-mers.
# Usage: bash process_reads.sh <R1.fastq[.gz]> <R2.fastq[.gz]> <sample> <out_dir>
#   FRACTION (default 0.05), K (default 31) and THREADS (default 8) can be set as
#   environment variables. Use the same FRACTION and K as the node k-mers.
# -> <out_dir>/<sample>_R{1,2}.clean.fastq, <sample>.fastp.{html,json},
#    <sample>_read_kmers_<fraction>_k<k>.tsv

set -euo pipefail

R1="$1"
R2="$2"
SAMPLE="$3"
OUTDIR="$4"
fraction=${FRACTION:-0.05}
k=${K:-31}
threads=${THREADS:-8}

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/bin/read_kmer_sampling"
[ -x "$BIN" ] || { echo "missing $BIN, run make in $ROOT first" >&2; exit 1; }
mkdir -p "$OUTDIR"

# Quality filtering (uncompressed output, since the sampler reads plain FASTQ)
CLEAN1="$OUTDIR/${SAMPLE}_R1.clean.fastq"
CLEAN2="$OUTDIR/${SAMPLE}_R2.clean.fastq"
fastp \
  -i "$R1" \
  -I "$R2" \
  -o "$CLEAN1" \
  -O "$CLEAN2" \
  --detect_adapter_for_pe \
  --cut_front \
  --cut_tail \
  --cut_mean_quality 20 \
  --length_required 30 \
  --thread "$threads" \
  --html "$OUTDIR/${SAMPLE}.fastp.html" \
  --json "$OUTDIR/${SAMPLE}.fastp.json"

# Sample read k-mers; R1 and R2 are streamed into one run
out="$OUTDIR/${SAMPLE}_read_kmers_${fraction}_k${k}.tsv"
"$BIN" "$fraction" "$k" <(cat "$CLEAN1" "$CLEAN2") "$out.tmp" \
  && [ -s "$out.tmp" ] && mv "$out.tmp" "$out" \
  || { echo "read k-mer sampling failed for $SAMPLE" >&2; exit 1; }

echo "Read k-mers: $out"
