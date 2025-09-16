#!/bin/bash

# Directories for input and output
CONTIGS_DIR="spades_contigs"
EXTRACTED_DIR="extracted_contigs"     # extract top contigs from each contigs.fasta into this folder
BLAST_RESULTS_DIR="blast_results"     # output
DB_NAME="listeria_db"

mkdir -p "$EXTRACTED_DIR"
mkdir -p "$BLAST_RESULTS_DIR"

# Number of contigs to extract
NUM_CONTIGS=5

# Loop through input directory, extract sample name
for sample_fasta in "$CONTIGS_DIR"/*.fasta; do
    sample_name=$(basename "$sample_fasta" .fasta)
    echo "Processing sample: $sample_name"

    # Extract top 5 longest contigs
    extracted_fasta="$EXTRACTED_DIR/${sample_name}_top${NUM_CONTIGS}_contigs.fasta"

    # Sorts the sequences in the input FASTA file by length -l, -r reverse the order (longest first)
   # then selects the top N sequences, and saves them to a new FASTA file.
    seqkit sort -l -r "$sample_fasta" | seqkit head -n $NUM_CONTIGS > "$extracted_fasta"

    # if the file doen't exist, skip
    if [[ ! -s "$extracted_fasta" ]]; then
        echo "No sequences extracted for $sample_name, skipping BLAST."
        continue
    fi

    # Run blastn
    blast_out="$BLAST_RESULTS_DIR/${sample_name}_blast.txt"
    blastn -query "$extracted_fasta" \
           -db "$DB_NAME" \
           -out "$blast_out" \
           -outfmt "6 qseqid sseqid pident length evalue bitscore stitle" \
           -max_target_seqs 5
    
    echo "BLAST done for $sample_name. Results: $blast_out"
done

echo "All done!"

