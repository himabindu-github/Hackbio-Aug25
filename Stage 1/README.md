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

### Step 5 : BLAST Analysis Workflow for Confirming the Organism Identity
To confirm the identity of the organism, I selected contigs longer than 500 base pairs from each assembled genome. These contigs were then compared against a custom Listeria reference database from NCBI using BLASTn. This process verifies that the assembled sequences match Listeria monocytogenes by identifying high similarity to reference sequences.
- All samples, including the outlier SRR27013337, showed strong matches to Listeria monocytogenes with high identity (~94.7–95%), confirming species identity.
- Although SRR27013337 had atypical assembly metrics in QUAST—such as fragmented contigs or irregular assembly statistics—BLAST was still able to detect Listeria sequences confidently. This is because BLAST can identify even small amounts of Listeria DNA within the sample by aligning contigs to reference genomes.


### Step 6 : Kraken2 Classification for Organism Identification
To further investigate sample composition and clarify any contamination or mixed species presence (especially in outlier samples like SRR27013337), the next step is to run Kraken for taxonomic classification. Kraken offers fast, accurate assignment of reads to taxonomic labels, helping to confirm Listeria dominance or detect other organisms.
#### 🦠 Key Findings:
-  Most samples showed over 70% of reads classified as Listeria monocytogenes, confirming the expected outbreak pathogen.
-  Sample SRR27013337 had approximately 85% Micrococcus and only 12% Listeria, indicating likely contamination or mislabeling; this sample will be excluded from downstream analyses.
-  A few samples showed minor contamination (<20%) with genera like Ralstonia and Streptococcus, which is acceptable for further analysis.
-  Around 10–20% of reads in most samples were classified as “Other” or “Unclassified”, typical in metagenomic data.
-  These unclassified reads are unlikely to significantly impact the analysis results.
Samples with moderate contamination will be included in downstream analyses such as antibiotic resistance profiling and sequence typing, but with caution. While the majority of sequences correspond to Listeria monocytogenes, the presence of contaminant reads may influence the interpretation of results by introducing misleading resistance gene profiles if contaminant species carry different resistance genes or they might affect sequence typing accuracy. Therefore, results from these samples will be interpreted carefully, considering possible contamination effects and cross-validating findings when possible.

#### Conclusion and Next Steps:
-	Exclude SRR27013337 from further analysis due to high contamination.
-	Proceed with downstream analyses for samples dominated by Listeria monocytogenes.
-	Samples with moderate contamination will be monitored closely during interpretation.
#### Additional Visualizations:
Sankey diagrams were generated using Pavian to visualize the taxonomic breakdown of each sample. These interactive visualizations clearly show how sequence reads are distributed among bacterial taxa, making it easier to:
-	Confirm the dominance of Listeria monocytogenes, and
-	Detect any contaminants or unusual patterns across samples.

### MLST Analysis Confirms Single Source of Listeria Outbreak

Multilocus Sequence Typing (MLST) was performed to understand the genetic relationship of Listeria monocytogenes from the samples in this study. This method analyzes specific housekeeping genes to identify the sequence type (ST) of bacterial strains, helping to confirm if cases are linked to the same source.
All samples in this study were identified as Sequence Type 6 (ST6), with a consistent allelic profile across seven key genes: abcZ(3), bglA(9), cat(9), dapE(3), dat(3), ldh(1), lhkA(5). This ST6 belongs to Clonal Complex 6 (CC6), which has been previously found in meat samples and is known to be associated with foodborne outbreaks.
The uniform presence of ST6/CC6 with this allelic profile in all samples indicates that the outbreak was caused by a single, clonal strain. This information is crucial for epidemiologists to track the spread, identify contamination sources, and develop targeted interventions to prevent future outbreaks. Using MLST improves outbreak investigation and enhances public health response by providing precise genetic information about the bacteria involved.

### Step 7 : Identification of Antimicrobial Resistance (AMR) Genes

The next step was to identify antimicrobial resistance (AMR) genes in the Listeria monocytogenes isolates used in this study. This helps us know which antibiotics might not work and how to better control the outbreak.

AMR gene detection was performed using Abricate with the CARD (Comprehensive Antibiotic Resistance Database). This analysis identified four resistance genes present in all 24 samples: fosX, lin, norB, and listeria_monocytogenes_mprF.

fosX: Confers resistance to fosfomycin by inactivating the antibiotic.

lin: Associated with resistance to lincosamide antibiotics through drug modification.

norB: A gene linked to fluoroquinolone resistance via efflux mechanisms.

listeria_monocytogenes_mprF: Confers resistance to peptide antibiotics by modifying membrane charge, reducing drug binding.

#### Prevalance of AMR genes:
All the samples used in this study showed 100% prevalence of the above 4 AMR genes, indicating widespread resistance potential to fluoroquinolone, fosfomycin, lincosamide, and peptide antibiotic classes within the outbreak strains.
The presence of these genes suggests that antibiotics from these classes may be less effective or ineffective for treating infections caused by these isolates. Careful antibiotic selection is essential to ensure effective clinical management.
Importantly, no resistance genes were detected that confer resistance to first-line drugs used for listeriosis treatment, such as ampicillin and gentamicin.
Treatment Recommendations:
Standard therapy with ampicillin combined with gentamicin remains effective and appropriate. Continued monitoring is recommended to detect any future emergence of resistance to these frontline antibiotics.

- We validated four resistance genes identified by Abricate using BLASTn against the CARD database. All genes showed high similarity (≥94% identity, ≥99% coverage), confirming their presence in the genomes. This step reduces the likelihood of false positives and supports the reliability of initial AMR screening

### Virulence Factor Analysis
To better understand the disease-causing potential of the Listeria monocytogenes isolates, we analyzed their virulence genes. This is important because virulence factors determine how well the bacteria can infect host cells, survive inside them, and spread, which impacts the severity of an outbreak.
The analysis was performed using Abricate with the VFDB (Virulence Factor Database), which identifies known virulence genes from genomic data.
- Genes hly, mpl, actA, plcA, plcB, and prfA correspond to Pathogenicity Island 1 (LIPI-1), which is essential for the intracellular life cycle of L. monocytogenes, enabling it to invade and multiply within host cells.
- The InlA-InlB locus, including inlA (Internalin A) and inlB (Internalin B), is involved in bacterial adhesion and entry into host cells, critical early steps in infection.
- Several genes from Listeria Pathogenicity Island 3 (LIPI-3), which codes for Listeriolysin S (LLS)—a bacteriocin and hemolytic cytotoxic factor—were identified, including llsA, llsH, llsX, llsB, llsY, llsD, and llsP. The presence of these genes indicates that the strains are highly virulent with enhanced capacity for cytotoxicity and survival in competitive environments.
Most virulence genes, including those critical for infection and intracellular survival, were present in nearly all isolates, underlining the high pathogenic potential of this outbreak strain.

### Discussion

- Our analysis shows that all Listeria monocytogenes samples come from the same strain (ST6), which is linked to meat products and can cause serious illness. This helps confirm the outbreak source and track how it spreads.

- We found genes that make the bacteria resistant to some antibiotics like fluoroquinolone and fosfomycin, but the usual treatments for listeriosis, such as ampicillin and gentamicin, should still work. It’s important to keep monitoring for any new resistance.

- The bacteria also carry several genes that help them invade and survive inside human cells, making this strain highly virulent.

- Genomic data helps turn complex information into clear insights that public health teams can use.

- Sharing this information in real-time supports outbreak response teams to track how the strain changes and spreads.

- This improves decision-making for controlling the outbreak and protecting the public.

- Overall, ongoing monitoring and clear communication are key to managing outbreaks, ensuring treatments remain effective, and keeping food products safe.













