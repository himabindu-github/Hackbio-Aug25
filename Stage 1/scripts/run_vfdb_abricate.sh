#!/bin/bash

#========================================
# Run Abricate with the vfdb database on SPAdes contigs (excluding SRR27013337).
# Results are appended to a combined TSV file with sample names prefixed.
#====================================

# Set input and outputs
INPUT_DIR="/home/maa/himabindu/results/SPAdes_output"
OUTPUT_FILE="/home/maa/himabindu/reports/abricate_resf_results.tsv"  # You can change this
DB="vfdb"
# Empty the output file first
> "$OUTPUT_FILE"

# loop through each subdirectory in the input folder
for subdir in "$INPUT_DIR"/*; do
  if [[ -d "$subdir" ]]; then
    CONTIG_FILE="$subdir/contigs.fasta"     #get the contigs.fasta file and assign it to a variable CONTIG_FILE
    SAMPLE_NAME=$(basename "$subdir")       # extract the sample name

    # Check if contigs.fasta exists and is not from SRR27013337
    if [[ -f "$CONTIG_FILE" && "$SAMPLE_NAME" != SRR27013337* ]]; then
      echo "Running abricate on: $SAMPLE_NAME"

      abricate --db "$DB" "$CONTIG_FILE" | sed "s/^/$SAMPLE_NAME\t/" >> "$OUTPUT_FILE"
    else
      echo "Skipping: $SAMPLE_NAME"
    fi
  fi
done
