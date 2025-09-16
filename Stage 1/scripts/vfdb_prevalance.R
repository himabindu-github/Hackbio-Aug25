# Load necessary library
library(dplyr)

# Step 1: Load the Abricate output file (VFDB results)
virulence_factor_data <- read.delim("abricate_vfdb_results.tsv", header = TRUE, sep = "\t")

# Step 2: Rename the first column to 'Sample' 
colnames(virulence_factor_data)[1] <- "Sample"

# Step 3: Clean and select relevant columns from the dataset
# Remove header rows accidentally included in the data (if present), and keep only necessary columns
virulence_factor_data_clean <- virulence_factor_data %>%
  select(Sample, GENE, X.COVERAGE, X.IDENTITY, ACCESSION, PRODUCT) %>%
  filter(
    GENE != "GENE",                          # Remove repeated header rows
    X.COVERAGE != "%COVERAGE",
    X.IDENTITY != "%IDENTITY",
    PRODUCT != "PRODUCT",
    ACCESSION != "ACCESSION"
  ) %>%
  rename(
    "%COVERAGE" = X.COVERAGE,
    "%IDENTITY" = X.IDENTITY
  )

# Step 4: Calculate the total number of unique samples
total_samples <- length(unique(virulence_factor_data_clean$Sample))

# Step 5: Summarise virulence factor prevalence across all samples
virulence_factor_summary <- virulence_factor_data_clean %>%
  group_by(GENE) %>%
  summarise(
    sample_count = n_distinct(Sample),                     # Number of unique samples carrying the gene
    prevalence = (sample_count / total_samples) * 100      # Prevalence as a percentage
  ) %>%
  arrange(desc(prevalence))                                # Sort by prevalence in descending order

# Step 6: Save the cleaned detailed data and summary to CSV files
write.csv(virulence_factor_data_clean, "abricate_vfdb_clean.csv", row.names = FALSE)
write.csv(virulence_factor_summary, "virulence_factor_prevalence.csv", row.names = FALSE)
