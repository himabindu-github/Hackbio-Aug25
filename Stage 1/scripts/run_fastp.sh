#!/bin/bash


#==========================================
#This script performs quality control (QC) and read trimming for raw WGS FASTQ files using FastQC, fastp, and MultiQC.
#==========================================

# -----------------------------------------------------
# STEP 1: Define directory paths for inputs and outputs
# -----------------------------------------------------
INPUT_DIR="raw_data_wgs"                         # Raw FASTQ files
OUTPUT_DIR="fastp_output"                        # Output directory for trimmed files and fastp reports
FASTQC_RAW_DIR="fastqc_raw_output"              # FastQC results for raw reads
FASTQC_TRIMMED_DIR="fastqc_trimmed_output"      # FastQC results for trimmed reads
MULTIQC_RAW_DIR="multiqc_raw_output"            # MultiQC summary for raw reads
MULTIQC_TRIMMED_DIR="multiqc_trimmed_output"    # MultiQC summary for trimmed reads

# -----------------------------------------------------
# STEP 2: Create all necessary output directories
# -----------------------------------------------------
mkdir -p "$OUTPUT_DIR" "$FASTQC_RAW_DIR" "$FASTQC_TRIMMED_DIR" "$MULTIQC_RAW_DIR" "$MULTIQC_TRIMMED_DIR"

# -----------------------------------------------------
# STEP 3: Run FastQC on raw FASTQ files
# -----------------------------------------------------
echo "Running FastQC on raw FASTQ files..."
fastqc "$INPUT_DIR"/*.fastq.gz -o "$FASTQC_RAW_DIR"

# -----------------------------------------------------
# STEP 4: Run MultiQC to summarize raw FastQC reports
# -----------------------------------------------------
echo "Summarizing raw FastQC reports with MultiQC..."
multiqc "$FASTQC_RAW_DIR" -o "$MULTIQC_RAW_DIR"

# Check the MultiQC report to decide trimming parameters and sample quality

# -----------------------------------------------------
# STEP 5: Run fastp to trim raw reads
# -----------------------------------------------------
echo "Running fastp for read trimming..."

for f1 in "$INPUT_DIR"/*_1.fastq.gz; do
  f2="${f1/_1.fastq.gz/_2.fastq.gz}"
  name=$(basename "$f1" _1.fastq.gz)

  echo "Processing sample: $name"

  fastp \
    -i "$f1" \
    -I "$f2" \
    -o "$OUTPUT_DIR/${name}_1.trimmed.fastq.gz" \
    -O "$OUTPUT_DIR/${name}_2.trimmed.fastq.gz" \
    -h "$OUTPUT_DIR/${name}_fastp.html" \
    -j "$OUTPUT_DIR/${name}_fastp.json" \
    --thread 4

  echo "Done: $name"
done

# -----------------------------------------------------
# STEP 6: Run FastQC on trimmed reads
# -----------------------------------------------------
echo "Running FastQC on trimmed FASTQ files..."
fastqc "$OUTPUT_DIR"/*_trimmed.fastq.gz -o "$FASTQC_TRIMMED_DIR"

# -----------------------------------------------------
# STEP 7: Summarize trimmed FastQC reports with MultiQC
# -----------------------------------------------------
echo "Performing MultiQC..."
multiqc "$FASTQC_TRIMMED_DIR" "$OUTPUT_DIR" -o "$MULTIQC_TRIMMED_DIR"

echo "Done!"
#echo "Check raw data QC: $MULTIQC_RAW_DIR/multiqc_report.html"
#echo "Check trimmed data QC: $MULTIQC_TRIMMED_DIR/multiqc_report.html"

