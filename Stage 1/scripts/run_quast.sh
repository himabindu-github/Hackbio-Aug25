#!/bin/bash

# -------------------------------------------------------------
# This script runs QUAST on all SPAdes assemblies in a directory.
# For each sample (subdirectory), it finds contigs.fasta and runs
# QUAST, saving the results in a separate output folder.
# -------------------------------------------------------------

#----------------------------------------------------------
# Step 1: Specify input and output directories
#--------------------------------------------------------
INPUT_DIR="SPAdes_output"
OUTPUT_DIR="quast_output"

#------------------------------------------------------------------
# Step 2: Create the output directory if it doesn't exist
#-----------------------------------------------------------------
mkdir -p "$OUTPUT_DIR"

echo "Starting QUAST analysis for assemblies in $INPUT_DIR"

#-----------------------------------------------------------------
# Step 3: Loop over each subdirectory inside the SPAdes output
#------------------------------------------------------------------
for subdir in "$INPUT_DIR"/* ; do
  # Step 3a: Check if it's a directory (to skip any non-folder files)
  if [[ -d "$subdir" ]]; then

    # Step 3b: Define the path to contigs.fasta
    CONTIG_FILE="$subdir/contigs.fasta"

    # Step 3c: Check if contigs.fasta exists
    if [[ -f "$CONTIG_FILE" ]]; then
      sample=$(basename "$subdir")
      echo "Running QUAST on: $sample"

      # Step 3d: Run QUAST and output results in a sample-named folder
      quast.py "$CONTIG_FILE" -o "$OUTPUT_DIR/$sample"

      echo "QUAST completed for: $sample"
      

    else
      echo "Skipping $subdir — contigs.fasta not found."
    fi

  fi
done

echo "All QUAST analyses completed."

