# Raw Data Processing

This directory contains the bioinformatics scripts and pipelines used to process the raw sequencing data. The workflows are separated into subdirectories for three distinct genetic markers: **COI**, **18S**, and **12S**.

These scripts handle the upstream processing steps from raw reads to annotated phyloseq objects, which are then utilized by the downstream statistical and modeling pipelines.

## Directory Structure & Contents

### `coi/` (Cytochrome c oxidase I pipeline)
This directory handles the bioinformatics processing for the COI marker.
* **Read Trimming:** Handled by `Cutadapt_260414_COI.sh`.
* **ASV Inference:** Runs DADA2 via `dada2_260414_parallel.R` and its corresponding SLURM submission script `run_dada2_parallel_260414.slurm`.
* **Taxonomy & Annotation:**
  * Uses BLAST against NCBI via `Blast_annotation.R`, `Blast_annotation_NCBI.R`, `Blast_annotation_ncbi_initial.sh`, and `Annotation_script_NCBI_260415.sh`.
  * Lowest Common Ancestor (LCA) assignment is processed using `lca.py` and `blastn_add_taxonomy_lite.py`.
  * The `add_taxonomy_scripts/` subdirectory contains Python scripts for querying various databases, including `gbif.py`, `bold.py`, `genbank.py`, `silva.py`, `unite.py`, and `privatebold.py`.
* **Output Generation:** Compiles the final data into an R object using `Create_phyloseq_object_NCBI.R`.
* **Utilities:** Includes auxiliary scripts like `headerandcoverage.py`.

### `18s/` (18S rRNA pipeline)
This directory processes the 18S marker data, commonly used for eukaryotic phytoplankton and zooplankton.
* **Read Trimming:** Executed via `Cutadapt_260414_18s.sh`.
* **ASV Inference:** Processed using `dada2_260415_parallel.R` and `run_dada2_parallel_260415.slurm`.
* **Taxonomy & Annotation:** Assigned against the PR² database using `taxonomy_18s_pr2.R` and `run_taxonomy_18s_pr2.slurm`.
* **Outputs:** Generates taxonomic tables such as `lca_taxon_table_n.csv`.

### `12s/` (12S rRNA pipeline)
This directory processes the 12S marker data, typically used for fish environmental DNA (eDNA).
* **Read Trimming:** Processed using `Cutadapt_260415.sh`.
* **ASV Inference:** Runs DADA2 using `dada2_260415_parallel.R` and `run_dada2_parallel_260415.slurm`.
* **Taxonomy & Annotation:** Uses `Blast_annotation_fish.R` and `Annotation_script_fish.sh` for alignment. LCA assignment is performed using `lca.py`.
* **Output Generation:** Compiles the data into phyloseq objects via `Create_phyloseq_object.R`. Includes output tables like `lca_taxon_table_n.csv`.

## Usage Notes
* Ensure that sequence data (FASTQ files) are in the expected directories before running the Cutadapt shell scripts.
* Taxonomic assignment scripts may require external database files (e.g., PR2, NCBI nt) to be downloaded and properly linked.
* The generated `.rds` phyloseq objects and CSV tables should be routed to the `data/raw/` and `data/processed/` directories for downstream analysis in the statistical pipeline.