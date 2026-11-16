# Independent Mesopredator Manuscript Reconstruction

**Version: 1.0** — *Validated standalone reconstruction bundle*

## Purpose
This folder contains the independently reconstructed analyses, figures, tables, processed datasets, and validation records for the Mesopredator manuscript. The reconstruction preserves approved results and records identified differences from the manuscript rather than altering data or results to force agreement.

The bundle has been run successfully from a fresh R session, and the exported figures and tables have been reviewed. Previously documented scientific deviations remain part of the final record.

## Running the Reconstruction
The complete bundle must contain its internal `rawfiles/`, `processed_data/`, and reference artifacts. Keep the folder name `00_final_reconstruction_scripts`.

In a fresh R session, set the working directory to this folder and run:

```R
source("RUN_ALL.R")
```

RUN_ALL.R executes the analysis stages in order and calls COLLECT_MANUSCRIPT_OUTPUTS.R to assemble the final manuscript outputs. No files from earlier project or reconstruction directories are required for the validated workflow.

### Analysis Scripts (Execution Order)
01_Sequencing_QC.R

02_Rarefaction_Threshold_Test.R

03_Phyloseq_Preprocessing.R

04_PERMANOVA_PCoA.R

04B_PERMANOVA_PCoA.R (sensitivity/assumptions audit)

05_Diet_Composition_FOO_RRA.R

06_GLLVM.R

07_Trophic_Network.R

07B_Trophic_Network_Sensitivity_Jackknife.R (sensitivity audit)

08_Stable_Isotope_Analysis.R

09_Functional_Group_Composition.R

10_Manuscript_Tables.R

### Output Structure
rawfiles/
Original pipeline input files, including phyloseq RDS objects, preserved metadata, and supporting source datasets.

figures/main_figures/
Reconstructed figures designated for the main manuscript.

figures/supplementary_figures/
Reconstructed supplementary figures and QC figures.

tables/main_tables/
Reconstructed main-manuscript tables.

tables/supplementary_tables/
Reconstructed supplementary tables and relevant status notes.

processed_data/
Derived phyloseq objects, intermediate outputs, and analysis-specific tables/figures. Includes internal, previously validated reference model objects and stable-isotope datasets.

validation/
Independent checks, comparison results, sensitivity analyses, figure-collection manifest, and documents explaining deviations. PERMANOVA (04B) and trophic-network jackknife (07B) have separate validation directories. The GLLVM methodological limitations and rationale for revised models are documented separately.

### Scientific Provenance and Limitations
#### Starting Inputs:
The phyloseq RDS files are the starting inputs for this analytical reconstruction. They are not raw FASTQ reads, and upstream sequence processing is outside this bundle's scope.

#### GLLVM Models:
The accepted GLLVM models are frozen, independently validated analytical artifacts in processed_data/reference_models/. Script 06 reads and verifies these models to reproduce accepted outputs; it does not refit them from the original data on every run.

#### Stable Isotopes:
The stable-isotope figure uses previously validated isotope tables in processed_data/reference_sia/. It does not repeat the upstream processing of original instrument measurements.

#### SCRuB Dependencies:
SCRuB is refitted during sequencing QC and requires its R-package and Torch/Lantern dependencies. Its results can be stochastic; the manuscript-reference sample-removal comparison retains a documented qualification.

#### Sensitivity Audits:
The 04B PERMANOVA and 07B network/jackknife scripts are independent diagnostic and sensitivity analyses. Their outputs are held in validation/ rather than presented as reconstructed manuscript figures.

#### Functional Groups & Tables:
The functional-group composition analysis retains documented differences in sample counts and taxonomic assignment. The manuscript-table comparison also retains documented mismatches. Successful execution does not imply exact agreement for these items.

#### Out of Scope Items:
The top-down modeling figure, certain population-dynamics figures, and the bioenergetics workflow illustration are outside the scope of the approved reconstruction scripts in this bundle; they are not silently substituted or fabricated.

#### Dependencies:
The bundle is independent of historical project folders, but it still requires appropriate installed R packages and their system dependencies. Package versions and the local R environment can affect reproducibility.

### Manifests and Interpretation
Consult validation/FIGURE_COLLECTION_MANIFEST.csv for the actual source-to-export mappings. Export names do not necessarily match original LaTeX figure filenames exactly. Consult the relevant files under validation/ for details of all identified deviations and methodological decisions.

### Release Policy
This version is the approved, tested reconstruction snapshot. Preserve its scripts, internal input/reference files, outputs, and validation records together. Any substantive change should be documented and released as a new version rather than silently modifying this one.