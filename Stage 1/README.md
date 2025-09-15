# Listeriosis Outbreak WGS Analysis 
Between 2017 and 2018, South Africa experienced a severe Listeriosis outbreak. To understand and control it, scientists leveraged whole-genome sequencing (WGS) to identify the pathogen, track its source, and guide treatment strategies. As bioinformatics students, we analyzed WGS data from over 100 bacterial samples collected during the outbreak. Our work aims to: 
-	Confirm the identity of the pathogen
-	Detect antimicrobial resistance (AMR) genes
-	Identify virulence factors such as toxins
-	Suggest effective antibiotics or treatment strategies

##  Objectives
1.	✅ Confirm the identity of the organism
2.	🧬 Determine antimicrobial resistance (AMR) profiles
3.	☠️ Detect potential toxins or virulence factors
4.	💊 Recommend potential treatments based on genomic evidence

## The Workflow
### Step 1: Data Download & Setup
This step involves collecting the raw sequencing data files (.fastq.gz) from a public repository (hackbio platform, in this case) using automated download commands. For this project, 50 samples are used. Organizing these files in a specific folder with clear names makes it easier to manage the data and avoid confusion. This is especially helpful when working with lots of samples, as it ensures everything is easy to find and use for later steps.
-	Input: List of curl commands to download .fastq.gz files from ENA
-	Tool: Custom bash script to read commands, assign filenames, and save to raw_data_wgs/ directory
-	Outcome: Organized 50 sequencing files ready for analysis

### Step 2: Quality Control
Quality control (QC) is vital to assess the raw sequencing data's reliability. Tools like FastQC and MultiQC provide visual and statistical summaries of read quality, adapter contamination, GC content, and duplication rates. This helps identify technical issues or biases that can affect downstream analyses. If problems are detected, trimming and filtering tools such as fastp or Trimmomatic are used to improve data quality by removing adapters, low-quality bases, and duplicates. High-quality input data ensures accurate genome assembly and reliable biological interpretations.
- Tools: FastQC, MultiQC
-	Purpose: Evaluate raw read quality — adapter content, read length, base quality, GC content
-	Actions: Trim adapters and low-quality bases using tools like fastp or Trimmomatic
-	Outcome: High-quality reads passed for downstream processing
#### Quality Control Summary
#### Initial QC
We analyzed raw reads from 50 samples using MultiQC, focusing on base quality, GC content, adapter contamination, and duplication levels. Most samples showed good base quality and expected GC content. However, some exhibited elevated duplication rates and minor adapter contamination—typical in outbreak datasets with closely related isolates. These results indicated the need for trimming and quality filtering before proceeding.

#### Post-Trimming QC
After applying fastp for adapter removal and quality trimming, we reassessed the data. The trimming process effectively reduced adapter content and duplication levels, while overall read quality improved across samples. This confirmed that the cleaned reads were of high quality and ready for downstream analysis.

#### Note on QC Warnings
The sequencing in this Listeria monocytogenes outbreak study was performed using the Nextera XT DNA Library Prep Kit, which uses tagmentation to simultaneously fragment DNA and add adapters. While efficient, this method can introduce sequence bias in the first ~12 bases of each read. As a result, QC tools like FastQC often generate warnings in the “Per base sequence content” module. These warnings are expected and not a sign of poor sequencing quality in this context.

### Step 3: Genome Assembly
Genome assembly reconstructs the full bacterial genomes from short sequencing reads. Using tools like SPAdes, the cleaned reads are pieced together into contiguous sequences (contigs) that represent the organism’s genome. High-quality assemblies are critical for identifying genetic features such as resistance genes, virulence factors, and for accurate species typing. This step transforms fragmented raw data into meaningful genomic structures ready for analysis.
-	Tool: SPAdes
-	Input: Cleaned paired-end reads
-	Output: Assembled contigs for each isolate
-	Purpose: Reconstruct genomes for typing and gene detection

### Step 4: Assembly Quality Assessment (QUAST Summary)
Assessing assembly quality using QUAST provides metrics such as genome size, N50 (contiguity), GC content, and number of contigs. These indicators help determine how complete and reliable each genome assembly is. Identifying outliers that may indicate contamination or misassemblies is important to ensure only high-quality genomes are used for further analysis.
-	Genome lengths: 3.02–3.13 Mb, consistent with L. monocytogenes
-	N50 values: ~60 kb to >250 kb
-	GC content: 37.8%–37.9%
-	Contig counts: mostly 30–100 per sample
-	Notable outlier: SRR27013337 likely contaminated/misclassified












