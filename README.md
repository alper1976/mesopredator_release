# Molecular diets of apex predators reveal trophic rewiring in a transitional coastal system
This repository contains the data, scripts, and modeling outputs for the mesopredator release reconstruction and trophic network analysis.

## Repository Structure

* **`data/`**: Contains all data required for the analysis.
  * `raw/`: Combined input data, unprocessed sequencing data, CN data, and species tables.
  * `processed/`: Combined processed data, intermediate outputs, and cleaned phyloseq objects.
* **`docs/`**: 
  * `modeling/`: Modeling documentation and methodology references.
  * `validation/`: Statistical analyses validation, missing figures logs, and manuscript scope.
* **`scripts/`**: 
  * `statistics/`: Scripts for statistical analyses, PERMANOVA, diet composition, and trophic networks. 
  * `modeling/`: Scripts for modeling (e.g., `Mesopredator_MASTER...R`).
  * `raw_data_processing/`: Scripts for raw data processing and formatting.
* **`results/`**: Final model outputs, figures, and compiled tables.

mesopredator_release-main/
├── README.md
├── data/
│   ├── raw/                 # Combined input data
│   └── processed/           # Combined processed data
├── docs/
│   ├── modeling/            # Modelling Documentation
│   └── validation/          # Statistical analyses validation
├── scripts/
│   ├── statistics/          # Scripts for statistical analyses
│   ├── modeling/            # Scripts for modeling
│   └── raw_data_processing/ # Scripts for raw data processing
└── results/                 # final results

## Usage

1. **Reconstruction**: Navigate to `scripts/reconstruction/` and execute `RUN_ALL.R` to process the raw sequencing and CN data.
2. **Modeling**: After reconstruction is complete, navigate to `scripts/modeling/` and run the `Mesopredator_MASTER` script to generate the final projections.
3. **Outputs**: Final tables and manuscript figures will be saved automatically to the `results/` directory.


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
