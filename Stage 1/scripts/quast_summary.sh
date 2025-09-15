#!/bin/bash

# Define input and output directories
INPUT_DIR="quast_output"
QUAST_SUMMARY="quast_summary.tsv"

# Write header for the summary table
echo -e "Sample\t# contigs\tContigs ≥1kb\tTotal length\tLargest contig\tN50\tGC (%)" > "$QUAST_SUMMARY"

# Loop through each QUAST report.tsv file
for REPORT in "$INPUT_DIR"/*/report.tsv; do
  if [[ -f "$REPORT" ]]; then

    # Extract sample name from parent directory
    SAMPLE=$(basename "$(dirname "$REPORT")")

    # Extract metrics from the report
    CONTIGS=$(grep -m 1 "^# contigs" "$REPORT" | cut -f2)
    CONTIGS_GT_KB=$(grep -m 1 "^# contigs (>= 1000 bp)" "$REPORT" | cut -f2)
    TOTAL_LENGTH=$(grep -m 1 "^Total length" "$REPORT" | cut -f2)
    LARGEST_CONTIG=$(grep -m 1 "^Largest contig" "$REPORT" | cut -f2)
    N50=$(grep -m 1 "^N50" "$REPORT" | cut -f2)
    GC=$(grep -m 1 "^GC (%)" "$REPORT" | cut -f2)

    # Append results to summary file
    echo -e "${SAMPLE}\t${CONTIGS}\t${CONTIGS_GT_KB}\t${TOTAL_LENGTH}\t${LARGEST_CONTIG}\t${N50}\t${GC}" >> "$QUAST_SUMMARY"
  else
    echo "Report.tsv not found in $REPORT"
  fi
done

echo "Summary of QUAST reports written to $QUAST_SUMMARY"
