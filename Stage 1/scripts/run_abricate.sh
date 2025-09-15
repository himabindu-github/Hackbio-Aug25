#!/bin/bash

INPUT_DIR="SPAdes_output"
OUTPUT_FILE="abricate_results.tsv"  # You can change this

# Empty the output file first (optional)
> "$OUTPUT_FILE"

for subdir in "$INPUT_DIR"/*; do
  if [[ -d "$subdir" ]]; then
    CONTIG_FILE="$subdir/contigs.fasta"
    SAMPLE_NAME=$(basename "$subdir")

    # Check if contigs.fasta exists and is not from SRR27013337
    if [[ -f "$CONTIG_FILE" && "$SAMPLE_NAME" != SRR27013337* ]]; then
      echo "Running abricate on: $SAMPLE_NAME"

      abricate "$CONTIG_FILE" | sed "s/^/$SAMPLE_NAME\t/" >> "$OUTPUT_FILE"
    else
      echo "Skipping: $SAMPLE_NAME"
    fi
  fi
done
