#!/bin/bash

#===============================
# Runs Abricate on contigs from each contigs.fasta file in SPAdes output folder (excluding SRR27013337).
# Appends results with sample names to a combined output TSV file.
#===============================

INPUT_DIR="/home/maa/himabindu/results/SPAdes_output"
OUTPUT_FILE="/home/maa/himabindu/reports/abricate_resf_results.tsv"  # You can change this
DB="card"

# Check if input directory exists
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: Input directory $INPUT_DIR does not exist. Exiting."
  exit 1
fi

# Check if abricate database exists
if ! abricate --list | grep -q "^$DB$"; then
  echo "Error: Abricate database '$DB' not found. Please download or install the database first."
  exit 1
fi

# Empty the output file first
> "$OUTPUT_FILE"

# Loop through each subdirectory in the input folder
for subdir in "$INPUT_DIR"/*; do
  if [[ -d "$subdir" ]]; then
    CONTIG_FILE="$subdir/contigs.fasta"     # Get the contigs.fasta file path
    SAMPLE_NAME=$(basename "$subdir")       # Extract the sample name

    # Check if contigs.fasta exists and sample is not SRR27013337 (excluded due to contamination)
    if [[ -f "$CONTIG_FILE" && "$SAMPLE_NAME" != SRR27013337* ]]; then
      echo "Running Abricate on: $SAMPLE_NAME"
      abricate --db "$DB" "$CONTIG_FILE" | sed "s/^/$SAMPLE_NAME\t/" >> "$OUTPUT_FILE"
    else
      echo "Skipping: $SAMPLE_NAME"
    fi
  fi
done

echo "Abricate screening complete. Results saved to $OUTPUT_FILE."
