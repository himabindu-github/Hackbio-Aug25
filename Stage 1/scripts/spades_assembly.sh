#!/bin/bash

INPUT_DIR="fastp_output"
OUTPUT_DIR="SPAdes_output"

mkdir -p "$OUTPUT_DIR"

echo "Assembly using spades"

for f1 in "$INPUT_DIR"/*_1.trimmed.fastq.gz; do
  f2="${f1/_1.trimmed.fastq.gz/_2.trimmed.fastq.gz}"
  sample=$(basename "$f1" _1.trimmed.fastq.gz)

  #echo "Assembling: $sample"

  spades.py --only-assembler \
    -1 "$f1" \
    -2 "$f2" \
    -o "$OUTPUT_DIR/$sample" 
   # --threads 4 \
   # --memory 16

  echo "Done assembling!"
done
