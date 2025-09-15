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
This step involves collecting the raw sequencing data files (.fastq.gz) from a public repository (hackbio platform, in this case) using automated download commands. Organizing these files into a dedicated directory is essential for efficient management and downstream processing. Properly naming and storing the files ensures traceability and prevents data loss, which is crucial when handling large datasets with many samples.
-	Input: List of curl commands to download .fastq.gz files from ENA
-	Tool: Custom bash script to read commands, assign filenames, and save to raw_data_wgs/ directory
-	Outcome: Organized 100+ sequencing files ready for analysis





