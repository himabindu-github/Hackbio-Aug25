#!/bin/bash

#========================================
# This script runs Kraken2 on all contig FASTA files in a folder.
# Outputs classification and report files for each sample.
#==========================================

# Path to downloaded Kraken2 database
DB="/Volumes/Crucial X6/wgs_analysis_hackbio/minikraken_8GB_20200312"

# Directory containing individual contig files
CONTIG_DIR="/Volumes/Crucial X6/wgs_analysis_hackbio/spades_contigs"

# Output directory
OUTPUT_DIR="/Volumes/Crucial X6/wgs_analysis_hackbio/kraken_results"
mkdir -p "$OUTPUT_DIR"

# Loop through all .fasta files in the contigs directory
for CONTIG in "$CONTIG_DIR"/*.fasta; do
  SAMPLE_NAME=$(basename "$CONTIG" .fasta)

  echo "Running Kraken2 on $SAMPLE_NAME..."

  kraken2 \
    --db "$DB" \
    --output "$OUTPUT_DIR/${SAMPLE_NAME}_kraken.txt" \
    --report "$OUTPUT_DIR/${SAMPLE_NAME}_report.txt" \
    --use-names \
    --threads 2 \
    "$CONTIG"

  echo "Finished $SAMPLE_NAME"
done

echo "All samples processed."

