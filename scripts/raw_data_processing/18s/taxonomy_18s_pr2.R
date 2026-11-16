library(dada2)
library(phyloseq)

# Load your finished ASV table
output_dir <- "/cluster/projects/nn9744k/02_results/51_Stensrud_students/tnt/18s/Figs" # Adjust if needed
seqtab.nochim <- readRDS(file.path(output_dir, "seqtab_nochim.rds"))

# Define database (Make sure this path is correct for your user!)
pr2_file <- "/cluster/projects/nn9745k/03_databases/pr2/5.0.0/pr2_version_5.0.0_SSU_dada2.fasta"

PR2_tax_levels <- c("Kingdom", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species")

# Run taxonomy (No mcparallel, just standard assignment)
print("Starting PR2 Taxonomy Assignment...")
taxa_pr2 <- assignTaxonomy(seqtab.nochim, refFasta = pr2_file, minBoot=80, taxLevels = PR2_tax_levels, multithread=TRUE, tryRC = TRUE)

# Save the result
saveRDS(taxa_pr2, file.path(output_dir, "pr2_tax_final.rds"))
print("Taxonomy successfully saved!"

print("Creating phyloseq object...")

# 1. Convert DADA2 tables into phyloseq components
# Note: DADA2 ASV tables have samples as rows, so we MUST specify taxa_are_rows = FALSE
ps_otu <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)
ps_tax <- tax_table(taxa_pr2)

# 2. Combine them into a single phyloseq object
ps <- phyloseq(ps_otu, ps_tax)

# 3. Save the final phyloseq object
saveRDS(ps, file.path(output_dir, "phyloseq_pr2_raw.rds"))
print("Phyloseq object successfully saved!"))
