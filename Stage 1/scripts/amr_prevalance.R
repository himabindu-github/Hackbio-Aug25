# This script processes Abricate CARD output to clean and summarize AMR gene data.
# It calculates resistance gene prevalence across samples and exports summary files.

# Load required library
library(dplyr)

# Step 1: Read in the CARD database results from Abricate
amr_card_data <- read.delim("abricate_card_results.tsv", header = TRUE, sep = "\t")

# Step 2: Rename the first column to 'Sample' 
colnames(amr_card_data)[1] <- "Sample"

# Step 3: Select relevant columns for analysis
amr_card_resgene_data <- amr_card_data %>%
  select(Sample, GENE, X.COVERAGE, X.IDENTITY, PRODUCT, RESISTANCE)

# Step 4: Clean the data
# Remove repeated header-like rows accidentally included during file generation
amr_card_resgene_data_clean <- amr_card_resgene_data %>%
  filter(
    GENE != "GENE",
    X.COVERAGE != "%COVERAGE",
    X.IDENTITY != "%IDENTITY",
    PRODUCT != "PRODUCT",
    RESISTANCE != "RESISTANCE"
  ) %>%
  rename(
    `%COVERAGE` = X.COVERAGE,
    `%IDENTITY` = X.IDENTITY
  )

# Step 5: Count total number of unique samples
total_samples <- n_distinct(amr_card_resgene_data_clean$Sample)

# Step 6: Group data by RESISTANCE class and compute summary statistics
amr_summary <- amr_card_resgene_data_clean %>%
  group_by(RESISTANCE) %>%
  summarise(
    sample_count = n_distinct(Sample),
    prevalence = (sample_count / total_samples) * 100
  ) %>%
  arrange(desc(prevalence))  # Sort by most prevalent resistance class

# Step 7: Export the results to CSV files
write.csv(amr_card_resgene_data_clean, "abricate_card_clean.csv", row.names = FALSE)
write.csv(amr_summary, "resistance_gene_prevalence.csv", row.names = FALSE)


