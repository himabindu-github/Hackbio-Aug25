#!/bin/bash

#=============================
# This script performs SPAdes assembly (de novo assembly) on fastq reads that are processed
#==============================


#------------------------------------------------------
# Step 1: Set input and output directories
#-----------------------------------------------------

INPUT_DIR="fastp_output"
OUTPUT_DIR="spades_output"

#-------------------------------------------------------
# Step 2: Create output directory if it doesn't exist
#--------------------------------------------------------

mkdir -p "$OUTPUT_DIR"

#------------------------------------------------------------
# Step 3: Set SPAdes parameters
#----------------------------------------------------------

THREADS=4
MEMORY=16             # in GB

echo "Assembly using spades"

#--------------------------------------------------------
# Step 4: Loop through all trimmed paired-end fastq files
#------------------------------------------------------------

for f1 in "$INPUT_DIR"/*_1.trimmed.fastq.gz; do
# Get corresponding read 2
  f2="${f1/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz}"
  sample=$(basename "$f1" _1.trimmed.fastq.gz)

  #echo "Assembling: $sample"

#----------------------------------------------------------
# Step 5: Run SPAdes
#-----------------------------------------------------------
  spades.py \
    -1 "$f1" \
    -2 "$f2" \
    -o "$OUTPUT_DIR/$sample" \
    --threads "$THREADS" \
    --memory "$MEMORY"

  echo "Done assembling!"
done


