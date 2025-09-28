# Bipolar Disorder Type II (BD-II) RNA-seq Analysis Project

## Brief Background

Bipolar Disorder Type II (BD-II) is a chronic psychiatric condition characterized by episodes of depression and hypomania. It significantly impacts mood regulation, cognition, and overall functioning. Despite advances in psychiatric genomics, the molecular underpinnings distinguishing familial BD-II (FBD) from sporadic BD-II (SBD) remain unclear. Understanding these differences can improve diagnosis, prognosis, and treatment strategies.

## Objectives

- Identify **differentially expressed genes (DEGs)** in:
  - SBD vs matched healthy controls (SHC)
  - FBD vs matched healthy controls (FHC)
  - SBD vs FBD
- Compare **molecular pathways** and **functional enrichment** in each group
- Assess **shared vs unique gene signatures**
- Investigate **synaptic**, **immune**, and **developmental** pathways
- Explore **biomarker candidates** and **therapeutic targets**

##  Methods Overview

| Step | Module | Description |
|------|--------|-------------|
| 1    | Data Acquisition | Download raw FASTQ files from SRA |
| 2    | Quality Control | Trim adapters and low-quality reads |
| 3    | Expression Quantification | STAR for indexing and alignment & featureCounts for quantification |
| 4    | Differential Expression | DESeq2 to identify DEGs between groups |
| 5    | Functional Analysis | GO/KEGG enrichment, UpSet plots, PCA |
| 6    | Integration | Compare SBD vs FBD vs Shared DEGs |
| 7    | Interpretation | Annotate biomarkers and pathways of interest |

## Requirements

- R / RStudio (with Bioconductor)
- R packages: `DESeq2`, `org.Hs.eg.db`, `clusterProfiler`, `ggplot2`, `pheatmap`
- STAR (for genome indexing and read alignment)
- featureCounts (for gene-level quantification)
- wget (for downloading raw sequencing data)



## 🔍 Module 1: Data Download, QC and Trimming

- **Tools**: `fastqc`, `fastp`
- Downloaded 8 RNA-Seq samples: 2 each from SBD, FBD, SHC, and FHC.
- Performed adapter removal and read quality trimming.

```
  #!/bin/bash

# ============================
# CONFIGURATION
# ============================

# Define base and output directories
BASE_DIR="/home/maa/himabindu/Hackbio-Aug25/Stage2/bipolar_rnaseq"
OUTPUT_DIR="$BASE_DIR/data/raw"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# ============================
# DOWNLOAD FASTQ FILES
# ============================

# Each wget downloads a paired-end FASTQ file into the output directory
wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/065/SRR33243165/SRR33243165_1.fastq.gz
wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/065/SRR33243165/SRR33243165_2.fastq.gz

wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/066/SRR33243166/SRR33243166_1.fastq.gz
wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/066/SRR33243166/SRR33243166_2.fastq.gz

wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/064/SRR33243164/SRR33243164_1.fastq.gz
wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/064/SRR33243164/SRR33243164_2.fastq.gz

wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/068/SRR33243168/SRR33243168_1.fastq.gz
wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/068/SRR33243168/SRR33243168_2.fastq.gz

wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/067/SRR33243167/SRR33243167_1.fastq.gz
wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/067/SRR33243167/SRR33243167_2.fastq.gz

wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/070/SRR33243170/SRR33243170_1.fastq.gz
wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/070/SRR33243170/SRR33243170_2.fastq.gz

wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/071/SRR33243171/SRR33243171_1.fastq.gz
wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/071/SRR33243171/SRR33243171_2.fastq.gz

wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/069/SRR33243169/SRR33243169_1.fastq.gz
wget -nc -P "$OUTPUT_DIR" ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR332/069/SRR33243169/SRR33243169_2.fastq.gz

# ============================
# DONE
# ============================

echo "Download complete. Files saved to: $OUTPUT_DIR"
```
```
#!/bin/bash

# ============================================================
# Script: run_fastqc.sh
# Purpose: Run FastQC on raw FASTQ files
# ============================================================

# 1. CONFIGURATION
BASE_DIR="/home/maa/himabindu/Hackbio-Aug25/Stage2/bipolar_rnaseq"
INPUT_DIR="$BASE_DIR/data/raw"
OUTPUT_DIR_FASTQC="$BASE_DIR/reports/qc_rawdata"

# 2. CHECK & CREATE OUTPUT DIRECTORY
echo "Creating FastQC output directory (if not present)..."
mkdir -p "$OUTPUT_DIR_FASTQC"

# 3. RUN FASTQC
echo "Running FastQC on raw FASTQ files in: $INPUT_DIR"

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory '$INPUT_DIR' does not exist."
    exit 1
fi

# Run FastQC on all fastq.gz files
fastqc "$INPUT_DIR"/*.fastq.gz -o "$OUTPUT_DIR_FASTQC"

echo "FastQC completed. Reports saved to: $OUTPUT_DIR_FASTQC"

```
```
#!/bin/bash
# ================================
# Script: run_fastp.sh
# Purpose: Run fastp on raw FASTQ files
# ================================

# Set paths to your directories
BASE_DIR="/home/maa/himabindu/Hackbio-Aug25/Stage2/bipolar_rnaseq"
INPUT_DIR="$BASE_DIR/data/raw"
OUTPUT_DIR="$BASE_DIR/data/trimmed"
REPORT_DIR="$BASE_DIR/reports/fastp"

# Create output directories if they don't exist
mkdir -p "$OUTPUT_DIR"
mkdir -p "$REPORT_DIR"

# Number of threads to use
THREADS=4

echo "Starting fastp quality trimming for paired-end FASTQ files..."

# Loop through all forward reads in INPUT_DIR
for R1 in "$INPUT_DIR"/*_1.fastq.gz; do
    # Derive the reverse read filename by replacing _1.fastq.gz with _2.fastq.gz
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"

    # Extract sample name from filename
    SAMPLE=$(basename "$R1" _1.fastq.gz)

    echo "Processing sample: $SAMPLE"

    # Run fastp with recommended parameters
    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$OUTPUT_DIR/${SAMPLE}_1.trimmed.fastq.gz" \
        -O "$OUTPUT_DIR/${SAMPLE}_2.trimmed.fastq.gz" \
        --detect_adapter_for_pe \
        --thread "$THREADS" \
        -h "$REPORT_DIR/${SAMPLE}_fastp.html" \
        -j "$REPORT_DIR/${SAMPLE}_fastp.json"

    echo "Finished trimming: $SAMPLE"
done

echo "All samples processed. Trimmed files in: $OUTPUT_DIR"
echo "Fastp reports saved to: $REPORT_DIR"

```
### Quality Control Summary

- **FastQC** revealed:
  - High sequence quality across all samples.
  - Presence of duplicated and overrepresented sequences in a few samples, particularly **SRR33243171_1**, which showed the highest duplication rate (**38.9%**).

- **Fastp** filtering results:
  - All samples retained **>42 million** reads after filtering, indicating good sequencing depth.
  - **GC content** ranged from **46.9% to 50.9%**, consistent with expected human transcriptome profiles.
  - **Adapter contamination was low**, with **<1.5% adapter content** in all samples.
  - Over **99% of reads passed filtering (PF)** in every sample, demonstrating high sequencing quality.

## ✅ Module 2: Alignment to Reference Genome

Following quality control, clean reads were aligned to the human reference genome to facilitate accurate transcript quantification.

### Objectives

- Download the **GRCh38 reference genome** and corresponding annotation (GTF).
- Build the **STAR index**.
- Prepare data for alignment in the next module.



### Script: Reference Genome Download & STAR Indexing

The following Bash script performs the full setup for STAR alignment:

---

```
#!/bin/bash

# Variables
BASE_DIR="/home/maa/himabindu/Hackbio-Aug25/Stage2/bipolar_rnaseq"
REF_DIR="$BASE_DIR/reference"
STAR_INDEX_DIR="$REF_DIR/STAR_index"
THREADS=8

# Create directories
mkdir -p "$REF_DIR"
mkdir -p "$STAR_INDEX_DIR"

# Move to reference directory
cd "$REF_DIR"

# Download genome fasta and gtf annotation
echo "Downloading reference genome and annotation..."
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

# Unzip files
echo "Unzipping reference files..."
gunzip -f GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip -f GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

# File paths
GENOME_FA="$REF_DIR/GCF_000001405.40_GRCh38.p14_genomic.fna"
GTF_FILE="$REF_DIR/GCF_000001405.40_GRCh38.p14_genomic.gtf"

# Build STAR index
echo "Building STAR index..."
STAR --runThreadN $THREADS \
     --runMode genomeGenerate \
     --genomeDir "$STAR_INDEX_DIR" \
     --genomeFastaFiles "$GENOME_FA" \
     --sjdbGTFfile "$GTF_FILE" \
     --sjdbOverhang 99

echo "STAR indexing completed!"
```

---

### Output

- Downloaded reference genome: `GCF_000001405.40_GRCh38.p14_genomic.fna`
- Annotation file: `GCF_000001405.40_GRCh38.p14_genomic.gtf`
- STAR genome index created in: `reference/STAR_index/`

---

> Next: Proceed to read alignment using STAR and quantify gene expression using **featureCounts**.



