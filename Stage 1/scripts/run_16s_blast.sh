#!/bin/bash

#==================================
# Extract 16S rRNA sequences using barrnap from each SPAdes contigs.fasta file,
# then BLAST the extracted 16S sequences against a 16S database for species confirmation.
#==================================

CONTIGS_DIR="spades_contigs"         # Input folder containing contigs.fasta files
EXTRACTED_16S_DIR="extracted_16s"    # Where extracted 16S sequences will be saved
BLAST_RESULTS_DIR="blast_16s_results" # Output folder for BLAST results

DB_16S="listeria_db"             # Your 16S BLAST database (makeblastdb output prefix)

mkdir -p "$EXTRACTED_16S_DIR"
mkdir -p "$BLAST_RESULTS_DIR"

# Loop over each contigs fasta
for contigs_fasta in "$CONTIGS_DIR"/*.fasta; do
    sample_name=$(basename "$contigs_fasta" .fasta)
    echo "Processing sample: $sample_name"

    barrnap_out_fasta="$EXTRACTED_16S_DIR/${sample_name}_16s.fasta"

    # Extract 16S rRNA sequences using barrnap's built-in output option
    barrnap --kingdom bac --outseq "$barrnap_out_fasta" "$contigs_fasta" 2> /dev/null

    if [[ ! -s "$barrnap_out_fasta" ]]; then
        echo "No 16S rRNA gene found for $sample_name — skipping 16S BLAST."
        continue
    fi

    # Run BLASTn on extracted 16S sequences
    blast_out="$BLAST_RESULTS_DIR/${sample_name}_16s_blast.tsv"
    blastn -query "$barrnap_out_fasta" \
           -db "$DB_16S" \
           -out "$blast_out" \
           -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
           -max_target_seqs 5

    echo "✅ BLAST complete for 16S rRNA of $sample_name. Results saved to $blast_out"

done

echo "All samples processed."

