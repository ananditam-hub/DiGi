#!/bin/bash

# Lift junctions from a GRCh38 BED to every assembly in the panel.
# Usage: bash liftover.sh <GRCh38.fa> <junctions.bed> <panel_dir> <out_dir>

set -euo pipefail

# Reference genome (GRCh38)
REF="$1"

# BED file containing junction coordinates on GRCh38
JUNCTIONS="$2"

# Directory of panel assembly FASTAs
PANEL_DIR="$3"

# Main output directory
OUTBASE="$4"

mkdir -p "$OUTBASE"

# Logs
date

# Run offline, and pull the pipeline once up front
export NXF_OFFLINE=1
nextflow pull evotools/nf-LO

index=0

# Skip extensions with no matching files
shopt -s nullglob

for var in "$PANEL_DIR"/*.{fa,fasta,fna,fas}
do
    base=$(basename "${var%.*}")

    echo "Processing: $base"

    # Separate output directory per assembly
    SAMPLE_OUT="$OUTBASE/$base"
    mkdir -p "$SAMPLE_OUT"

    # Run nf-LO liftover
    nextflow run evotools/nf-LO \
        --source "$REF" \
        --target "$var" \
        --distance near \
        --aligner minimap2 \
        -profile conda \
        --outdir "$SAMPLE_OUT"

    # Chain file produced by nf-LO
    CHAIN="$SAMPLE_OUT/chainnet/liftover.chain"

    # Output lifted BED
    OUTPUT_BED="$SAMPLE_OUT/${base}_junctions.bed"

    # Run CrossMap
    CrossMap region \
        "$CHAIN" \
        "$JUNCTIONS" \
        > "$OUTPUT_BED" \
        2> "$SAMPLE_OUT/crossmap.log"

    index=$((index + 1))

    echo "Finished: $base"
done

echo "Processed $index assemblies."
