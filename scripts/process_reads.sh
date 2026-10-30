#!/bin/bash

# Quality-filter MHC-filtered paired-end reads and count their
# FracMinHash-sampled k-mers.
# Usage: bash process_reads.sh <R1.fastq[.gz]> <R2.fastq[.gz]> <sample> <out_dir>
#   FRACTION (default 0.05), K (default 31), THREADS (default 8), KMC_MEM in GB
#   (default 16) and KMC_TMP (dir for KMC's bins/db, default <out_dir>) can be set
#   as environment variables. Use the same FRACTION and K as the node k-mers.
# -> <out_dir>/<sample>_R{1,2}.clean.fastq.gz, <sample>.fastp.{html,json},
#    <sample>_read_kmers_<fraction>_k<k>.tsv

set -euo pipefail

R1="$1"
R2="$2"
SAMPLE="$3"
OUTDIR="$4"
fraction=${FRACTION:-0.05}
k=${K:-31}
threads=${THREADS:-8}
kmc_mem=${KMC_MEM:-16}

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/bin/kmer_set_sampling"
[ -x "$BIN" ] || { echo "missing $BIN, run make in $ROOT first" >&2; exit 1; }
for tool in fastp kmc kmc_tools; do
  command -v "$tool" >/dev/null || { echo "missing $tool, activate the digi env" >&2; exit 1; }
done
mkdir -p "$OUTDIR"

# Quality filtering
CLEAN1="$OUTDIR/${SAMPLE}_R1.clean.fastq.gz"
CLEAN2="$OUTDIR/${SAMPLE}_R2.clean.fastq.gz"
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

# Count canonical k-mers of R1 + R2 with KMC, then sample them. KMC drops k-mers
# containing N and canonicalises as the lexicographic min, like node_kmers_gen.
# Sampling depends on the k-mer alone, so counting everything first gives the
# same counts as sampling while counting.
work=$(mktemp -d "${KMC_TMP:-$OUTDIR}/kmc.${SAMPLE}.XXXXXX")
trap 'rm -rf "$work"' EXIT
mkdir -p "$work/bins"
readlink -f "$CLEAN1" "$CLEAN2" > "$work/reads.lst"

# -ci1: keep singletons (KMC default -ci2 drops them); -cs: no 255 cap on counts
kmc -k"$k" -ci1 -cs4294967295 -t"$threads" -m"$kmc_mem" -fq @"$work/reads.lst" "$work/db" "$work/bins" \
  > "$work/kmc.log" 2>&1 || { cat "$work/kmc.log" >&2; echo "kmc failed for $SAMPLE" >&2; exit 1; }
grep -E "No. of unique counted k-mers|Total no. of k-mers|Total no. of reads" "$work/kmc.log" >&2

# Stream the dump straight into the sampler instead of writing it out
out="$OUTDIR/${SAMPLE}_read_kmers_${fraction}_k${k}.tsv"
kmc_tools -t"$threads" transform "$work/db" dump /dev/stdout 2> "$work/dump.log" \
  | "$BIN" "$fraction" /dev/stdin "$out.tmp" \
  && [ -s "$out.tmp" ] && mv "$out.tmp" "$out" \
  || { cat "$work/dump.log" >&2; echo "read k-mer sampling failed for $SAMPLE" >&2; rm -f "$out.tmp"; exit 1; }

echo "Read k-mers: $out"
