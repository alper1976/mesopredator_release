# Raw Data

This directory contains the combined raw input data, metadata, and baseline R data objects required to run the statistical analysis, trophic reconstruction, and bioenergetic modeling pipelines.

## Contents

### 1. Phyloseq Objects (`*.rds`)
These are the baseline, filtered, and rarefied R data objects (RDS files) fed into the downstream compositional and trophic analyses:
* `clean_phyloseq_mifish_lca.rds`
* `clean_phyloseq_mifish_lca_hostadded.rds`
* `ps_1.rds`
* `ps_coi_metazoa.rds`
* `ps_decon_mifish.rds`
* `ps_filtered_mifish.rds`
* `ps_pruned_mifish_100reads.rds`
* `ps_pruned_mifish_500reads.rds`
* `ps_rarefied_mifish_100reads.rds`
* `ps_with_host.rds`

### 2. Input Datasets & Metadata
* **Metadata**: Files like `Metadata_TnT_NEW.csv` and `plasma_metadata.csv` detailing sample information.
* **Taxa Lists**: Taxonomic assignment tables (e.g., `taxa_list_mifish.csv`, `taxa_list_coi_full.csv`).
* **Modeling Inputs**: External data arrays, CPUE datasets, and species tables (e.g., `wrasse_cpue_skagerrak.csv`, `species_table_n.csv`) utilized by the main `Mesopredator_MASTER` script.
* **Isotope/Diet Data**: Raw CN and plasma datasets (e.g., `CN data Eiler_April2026.csv`, `plasma_data.csv`).