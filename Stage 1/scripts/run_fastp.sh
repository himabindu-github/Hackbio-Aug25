#!/bin/bash

INPUT_DIR="raw_data_wgs"
OUTPUT_DIR="fastp_output"
FASTQC_RAW_DIR="fastqc_raw_output"
FASTQC_TRIMMED_DIR="fastqc_trimmed_output"
MULTIQC_RAW_DIR="multiqc_raw_output"
MULTIQC_TRIMMED_DIR="multiqc_trimmed_output"

mkdir -p "$OUTPUT_DIR" "$FASTQC_RAW_DIR" "$FASTQC_TRIMMED_DIR" "$MULTIQC_RAW_DIR" "$MULTIQC_TRIMMED_DIR"

echo "Running FastQC on raw FASTQ files..."

fastqc "$INPUT_DIR"/*.fastq.gz -o "$FASTQC_RAW_DIR"

echo "Performing MultiQC with fastqc_raw_outpu..."
multiqc "$FASTQC_RAW_DIR" -o "$MULTIQC_RAW_DIR"

# "Check the MultiQC report to decide trimming parameters and sample quality."

echo "Running Fastp for read trimming..."

# loop through all forward read files (-1.fastq.gz) in the input directory
for f1 in "$INPUT_DIR"/*_1.fastq.gz; do
  # Get the reverse read files by replacing -1 with -2
  f2="${f1/_1.fastq.gz/_2.fastq.gz}"

  # retrieves the sample name for naming the output
  name=$(basename "$f1" _1.fastq.gz)

  echo "🔄 Processing sample: $name"

  fastp \
    -i "$f1" \
    -I "$f2" \
    -o "$OUTPUT_DIR/${name}_1.trimmed.fastq.gz" \
    -O "$OUTPUT_DIR/${name}_2.trimmed.fastq.gz" \
    -h "$OUTPUT_DIR/${name}_fastp.html" \
    -j "$OUTPUT_DIR/${name}_fastp.json" \
    --thread 4

  echo "Done!"

done

echo "Running FastQC on trimmed FASTQ files..."

fastqc "$OUTPUT_DIR"/*_trimmed.fastq.gz -o "$FASTQC_TRIMMED_DIR"

echo "Performing MultiQC on trimmed data"
multiqc "$FASTQC_TRIMMED_DIR" "$OUTPUT_DIR" -o "$MULTIQC_TRIMMED_DIR"

echo "All processing completed!"
