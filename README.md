# Molecular diets of apex predators reveal trophic rewiring in a transitional coastal system
This repository contains the code, configuration files, and analytical pipelines used for the manuscript: "Molecular diets of apex predators reveal dietary partitioning and exploitation of mesopredator release in a transitional coastal system" by Røstum et al. (2026).
## Overview
This project investigates the multi-trophic restructuring of the Skagerrak and Oslofjord coastal ecosystems over a 14-year period (2012–2025). By integrating multi-marker dietary DNA metabarcoding (MiFish 12S and COI Leray-XT), stable isotope analysis ($\delta^{15}N$, $\delta^{13}C$), and stochastic bioenergetic modeling, we track the foraging shifts of two expanding apex predators: the harbour seal (Phoca vitulina) and the great cormorant (Phalacrocorax carbo).
## Key Features of this Repository:
### Bioinformatics Pipeline: 
Amplicon processing via Cutadapt and DADA2, with spatial cross-contamination removal using the SCRUB algorithm.

### Statistical Data Analysis: 
Robust centered log-ratio (rCLR) transformations, PERMANOVA, and temporal modeling using Generalised Linear Latent Variable Models (GLLVMs).

### Bioenergetic Monte Carlo Simulation: 
A unified 5,000-draw stochastic simulation framework predicting absolute biomass consumption and top-down impact ratios, integrating:
ICES Stock Assessment API (icesSAG) for pelagic and demersal standing stocks.
Fishery-independent CPUE integrations (IMR Beach Seine Survey) for data-limited mesopredators (Labridae).

### Repository Structure
├── Bioinformatics/<br>
│   ├── 01_cutadapt_trimming.sh<br>          # Primer removal
│   ├── 02_dada2_pipeline.R<br>              # ASV inference and chimaera removal
│   ├── 03_lca_taxonomy_assignment.py<br>    # Lowest Common Ancestor (LCA) assignment
│   └── 04_scrub_decontamination.R<br>       # SCRUB implementation for plate-well leakage
├── Data/<br>
│   ├── Inputs/<br>                          # Raw inputs for Monte Carlo models
│   │   ├── Metadata_TnT_main.csv<br>        # Sample metadata
│   │   ├── species_table_n.csv<br>          # Filtered ASV abundance matrices
│   │   ├── wrasse_cpue_skagerrak.csv<br>    # IMR Flødevigen beach-seine CPUE index
│   │   └── ...<br> 
├── Scripts/<br>
│   ├── 01_compositional_stats.R<br>         # rCLR, PCoA, PERMANOVA (vegan)
│   ├── 02_gllvm_temporal_models.R<br>       # Taxon-level GLVM modeling (gllvm)
│   ├── 03_trophic_network.R<br>             # Network construction (igraph/ggraph)
│   └── Mesopredator_MASTER_v1_4.R<br>       # Primary Monte Carlo bioenergetic simulation
└── Output/<br>
    ├── Figures/<br>                         # Generated plots (ggplot2)
    └── Tables/<br>                          # Exported QC tables and MC summaries

## Prerequisites & Dependencies
### R Packages
The statistical and bioenergetic modeling scripts were developed using R (version 4.x). Ensure the following core packages are installed:

install.packages(c("tidyverse", "vegan", "gllvm", "MASS", "igraph", "ggraph", "patchwork"))

For dynamic stock extraction, you will need the ICES Stock Assessment Graph API package:

install.packages("icesSAG")

### Bioinformatic Tools
Cutadapt (v4.4+)
DADA2 (v1.28+)
SCRUB (Source-tracking for Contamination Removal in microBiomes)

### Data Availability
Raw Sequence Data: Demultiplexed .fastq files for all 172 samples are available on the NCBI Sequence Read Archive (SRA) under BioProject PRJNA[XXXXXX].

### Reference Databases: 
The ScandiFish database (version 2.0) used for MiFish assignments is available at [Link/DOI].

### Stable Isotope Baselines: 
Regional monitoring baselines extracted from the Norwegian Water Environment database (vannmiljo.no).

## Workflow Instructions
### 1. Bioinformatics & ASV Profiling
Navigate to the Bioinformatics/ directory. Run the shell and R scripts sequentially. The pipeline outputs clean ASV count tables. A strict $0.1\%$ global relative abundance threshold is applied prior to robust Aitchison transformations to filter transient sequencing noise without artificial sample rarefaction.
### 2. Compositional & Trophic Analysis
Run Scripts/01_compositional_stats.R and 02_gllvm_temporal_models.R to reproduce the beta-diversity ordinations and taxon-level temporal trajectories. Secondary prey ingestion networks are reconstructed using 03_trophic_network.R.3. 

### 3. Bioenergetic Modeling
The main simulation is housed in Scripts/Mesopredator_MASTER_v1_4.R.
This script automatically downloads the latest ICES stock assessment parameters via the API (including the multi-substock spatial approximations for recent cod demographics). It integrates the digitized IMR Flødevigen beach-seine CPUE to establish effort-corrected wrasse standing stocks. The script runs a 5,000-iteration Monte Carlo loop to propagate uncertainties across diet model coefficients, daily energetic mass intake ($D$), regional residency time ($T$), and prey biomass approximations.Note: Requires internet access to ping the icesSAG database .

## Citation
If you use this code or data in your work, please cite the corresponding paper:@article{rostum2026molecular,
  title={Molecular diets of apex predators reveal dietary partitioning and exploitation of mesopredator release in a transitional coastal system},
  author={R{\o}stum, Tarjei and Angell, Thorvald and Werner, Even and Lundstr{\"o}m, Karl and Knutsen, Halvor and Stensrud, Eivind and Eiler, Alexander},
  journal={Global Change Biology},
  year={2026},
  doi={10.1111/gcb.XXXXX}
}
