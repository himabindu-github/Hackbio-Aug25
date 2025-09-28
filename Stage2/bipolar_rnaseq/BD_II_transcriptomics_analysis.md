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
---
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
---
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
---
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
---
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
## Module 3: Read Alignment to Reference Genome

After preparing the STAR index, the next step involves aligning the **trimmed paired-end RNA-seq reads** to the human reference genome. STAR (Spliced Transcripts Alignment to a Reference) was used for its speed and splice-aware alignment capabilities.

---

###  Objectives

- Align each sample to the GRCh38 reference genome using **STAR**.
- Generate **sorted BAM files** and alignments annotated with all SAM attributes.
- Verify successful alignment of all reads.

---

### Script: STAR Alignment of Trimmed Reads

The following Bash script automates the alignment of all paired-end FASTQ files from the trimmed dataset:

```bash
#!/bin/bash

# STAR Alignment Script for paired-end RNA-seq reads
# Aligns trimmed fastq.gz files to a pre-built STAR genome index
# Outputs sorted BAM and gene counts

# Variables
BASE_DIR="/home/maa/himabindu/Hackbio-Aug25/Stage2/bipolar_rnaseq"
FASTQ_DIR="$BASE_DIR/data/trimmed"        # Location of trimmed fastq.gz files
STAR_INDEX_DIR="$BASE_DIR/reference/STAR_index"
OUTPUT_DIR="$BASE_DIR/STAR_alignment_output"
THREADS=6

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Check if STAR index directory exists
if [ ! -d "$STAR_INDEX_DIR" ]; then
    echo "Error: STAR index directory not found at $STAR_INDEX_DIR"
    echo "Please build the STAR index first."
    exit 1
fi

# Check if fastq files exist
if [ ! "$(ls $FASTQ_DIR/*_1.trimmed.fastq.gz 2>/dev/null)" ]; then
    echo "Error: No R1 trimmed fastq.gz files found in $FASTQ_DIR"
    exit 1
fi

# Loop through all R1 reads
for r1 in "$FASTQ_DIR"/*_1.trimmed.fastq.gz; do
    sample=$(basename "$r1" | sed 's/_1\.trimmed\.fastq\.gz//')
    r2="$FASTQ_DIR/${sample}_2.trimmed.fastq.gz"

    # Check if paired R2 exists
    if [ ! -f "$r2" ]; then
        echo "Warning: Paired file for $sample not found, skipping sample."
        continue
    fi

    echo "Starting alignment for sample: $sample"

    # Define output prefix and output directory per sample
    SAMPLE_OUT_PREFIX="$OUTPUT_DIR/${sample}_"
    
    # Run STAR alignment
    STAR --runThreadN $THREADS \
         --genomeDir "$STAR_INDEX_DIR" \
         --readFilesIn "$r1" "$r2" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$SAMPLE_OUT_PREFIX" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes All

    if [ $? -eq 0 ]; then
        echo "Alignment completed successfully for $sample"
    else
        echo "Error during alignment of $sample"
    fi

done

echo "All samples processed."
```

---

### Output Files per Sample

For each sample, STAR produces:

- Sorted BAM file
- Alignment summary
- Splice junctions

---

> Next: Proceed to **quantify gene expression** using `featureCounts` to generate count matrices.
## Module 4: Gene Expression Quantification

Following successful alignment with STAR, this module focuses on **quantifying gene-level expression** using `featureCounts`, a high-performance read summarization tool.

---

### Objectives

- Generate a **gene-level count matrix** from aligned BAM files.
- Use the **GTF annotation** for accurate gene models (from GRCh38.p14).
- Prepare expression data for downstream differential expression analysis.

---

### Script: featureCounts (Quantification)

```bash
#!/bin/bash

# Exit immediately on error
set -e

# Paths
BASE_DIR="/home/maa/himabindu/Hackbio-Aug25/Stage2/bipolar_rnaseq"
ANNOTATION="$BASE_DIR/reference/GCF_000001405.40_GRCh38.p14_genomic.gtf"
OUTPUT_DIR="$BASE_DIR/counts"
BAM_DIR="$BASE_DIR/STAR_alignment_output"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Debug info
echo "Using annotation file: $ANNOTATION"
echo "Looking for BAM files in: $BAM_DIR"
ls "$BAM_DIR"/*_Aligned.sortedByCoord.out.bam || { echo "No BAM files found!"; exit 1; }

# Run featureCounts
CMD="featureCounts -T 4 -p -O -t gene -g gene_id -a \"$ANNOTATION\" -o \"$OUTPUT_DIR\"/counts.txt \"$BAM_DIR\"/*.bam"

echo "Running command:"
echo $CMD

eval $CMD  

echo "featureCounts completed successfully."
```

---

###  Output

- `counts.txt`: A matrix of **raw read counts per gene**, per sample.
- `counts.txt.summary`: Summary of reads assigned and unassigned.
- Counts are grouped using the `gene_id` field from the GTF annotation.
- These counts will be **used in DESeq2** for normalization and statistical testing.

---

> Proceed to **Module 5: Differential Expression Analysis**.

## Module 5: DESeq2 Differential Expression Analysis & Visualization Pipeline

This document provides a full RNA-seq analysis workflow using **DESeq2**, including differential gene expression, PCA, volcano plots, and heatmaps.

---

### 1. Load Required Libraries

```r
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(factoextra)
library(ggrepel)
```

---

###  2. Load and Prepare Count Data

```r
# Load count matrix
raw_counts <- read.delim("counts/counts.txt", header = TRUE)

# Set gene IDs as row names
rownames(raw_counts) <- raw_counts$Geneid

# Select only count columns
count_data <- raw_counts[, 7:ncol(raw_counts)]

# Ensure numeric format
count_data <- as.data.frame(lapply(count_data, as.numeric))
rownames(count_data) <- raw_counts$Geneid

# Preview
head(count_data)
```

---

### 3. Define Sample Metadata

```r
metadata <- data.frame(
  SampleID = c("SRR33243164", "SRR33243165", "SRR33243166", "SRR33243167",
               "SRR33243168", "SRR33243169", "SRR33243170", "SRR33243171"),
  Condition = c("SBD", "SBD", "FBD", "FBD", "SHC", "SHC", "FHC", "FHC")
)

# Match order
metadata <- metadata[match(colnames(count_data), metadata$SampleID), ]
stopifnot(all(metadata$SampleID == colnames(count_data)))

# Define colData for DESeq2
colData <- data.frame(
  sample = colnames(count_data),
  condition = factor(metadata$Condition, levels = c("SHC", "SBD", "FBD", "FHC"))
)
```

---

### 4. Differential Expression Analysis (DESeq2)

```r
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = colData, design = ~ condition)
dds <- DESeq(dds)

# Results
res_sbd_shc <- results(dds, contrast = c("condition", "SBD", "SHC"))
res_fbd_fhc <- results(dds, contrast = c("condition", "FBD", "FHC"))
res_sbd_fbd <- results(dds, contrast = c("condition", "SBD", "FBD"))

# Top DEGs
head(res_sbd_shc[order(res_sbd_shc$padj), ])
head(res_fbd_fhc[order(res_fbd_fhc$padj), ])
head(res_sbd_fbd[order(res_sbd_fbd$padj), ])
```

---

### 5. PCA Analysis

```r
vsd <- vst(dds, blind = FALSE)
pca_res <- prcomp(t(assay(vsd)))

# Variance explained
summary(pca_res)

# Test variation by condition
condition <- colData(dds)$condition
summary(aov(pca_res$x[,1] ~ condition))
summary(aov(pca_res$x[,2] ~ condition))

# Plot
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Samples") +
  theme_minimal()
```

---

###  6. PCA Colored by Age

```r
age <- c(33, 28, 36, 38, 24, 25, 45, 67)
metadata <- data.frame(
  Condition = colData(dds)$condition,
  Age = age
)

fviz_pca_ind(pca_res,
             geom.ind = "point",
             col.ind = age,
             addEllipses = FALSE,
             legend.title = "Age",
             repel = TRUE)

# Correlations
cor.test(age, pca_res$x[, 1])
cor.test(age, pca_res$x[, 2])
```

---

###  7. Volcano Plots

#### a. SBD vs SHC

```r
res_df_s <- as.data.frame(res_sbd_shc) %>%
  mutate(gene = rownames(.)) %>%
  filter(!is.na(padj)) %>%
  mutate(sig = case_when(
    padj < 0.1 & log2FoldChange > 1 ~ "Upregulated",
    padj < 0.1 & log2FoldChange < -1 ~ "Downregulated",
    TRUE ~ "Not significant"
  ))

ggplot(res_df_s, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.7, size = 1.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot: SBD vs SHC")
```

_Repeat for FBD vs FHC and SBD vs FBD._

---

### 8. Heatmaps of DEGs

#### a. SBD vs SHC

```r
deg_sbd_shc <- res_sbd_shc %>% as.data.frame() %>%
  filter(padj < 0.1 & abs(log2FoldChange) > 1)

deg_genes <- rownames(deg_sbd_shc)
norm_counts <- counts(dds, normalized = TRUE)
log_norm_counts_deg <- log2(norm_counts[deg_genes, ] + 1)

sample_annot <- as.data.frame(colData(dds)[, "condition", drop = FALSE])
colnames(sample_annot) <- "Condition"

pheatmap(log_norm_counts_deg,
         annotation_col = sample_annot,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Heatmap: SBD vs SHC DEGs")
```

#### b. FBD vs FHC

_Repeat using `res_fbd_fhc`._

---

###  9. Combined Heatmap (SBD vs SHC & FBD vs FHC)

```r
# Extract DEGs
deg_sbd_shc <- res_sbd_shc %>% as.data.frame() %>% filter(padj < 0.1 & abs(log2FoldChange) > 1)
deg_fbd_fhc <- res_fbd_fhc %>% as.data.frame() %>% filter(padj < 0.1 & abs(log2FoldChange) > 1)

combined_genes <- union(rownames(deg_sbd_shc), rownames(deg_fbd_fhc))

# Log-normalized matrix
norm_counts <- counts(dds, normalized = TRUE)
log_deg_counts <- log2(norm_counts[combined_genes, ] + 1)

# Rename columns
colnames(log_deg_counts) <- c("SBD-1", "SBD-2", "FBD-1", "FBD-2", "SHC-1", "SHC-2", "FHC-1", "FHC-2")
rownames(sample_annot) <- colnames(log_deg_counts)

# Heatmap
pheatmap(log_deg_counts,
         annotation_col = sample_annot,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Combined Heatmap: SBD vs SHC & FBD vs FHC")
```

---

## End of Pipeline








