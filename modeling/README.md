# Mesopredator analysis — collaborator run package

## Purpose
This minimal package contains the revised 2012–2025 mesopredator master R script, its three local input datasets, archived ICES cod assessment tables used as an API fallback, and essential methodological documentation. It intentionally excludes earlier scripts, intermediate results, plots, and large raw DATRAS datasets.

## Important: script version in this package
The distributed master is a **copy of v1.4 with one requested parameter update**:
`survival <- exp(-0.55)` (about 57.7% annual age-1-to-age-2 survival).
The uploaded v1.4 file had `survival <- 0.553`. The copied script was updated to match the agreed literature-informed mortality assumption. The original uploaded file was not overwritten. Results must be regenerated; existing plots/tables do not reflect this update.

## Prerequisites
- R and internet access to the ICES Stock Assessment Graphs API.
- Install these R packages (including packages needed by later sections):

```r
install.packages(c(
  "dplyr", "tidyr", "readr", "ggplot2", "tibble",
  "MASS", "icesSAG", "stringr", "patchwork"
))
```

## Run instructions
1. Extract the ZIP to a **local writable folder**. Do not run the script directly from inside the ZIP.
2. Keep all folders and filenames exactly as supplied. The script expects:
   - `Data/Inputs/Metadata_TnT - metadata_main.csv`
   - `Data/Inputs/species_table_n.csv`
   - `Data/Inputs/område og art_data.csv` (UTF-16 tab-delimited, decimal comma)
   - `final result/Inputs/summary_17652.csv`, `summary_19661.csv`, `summary_19662.csv`, `summary_19663.csv` (fallback only for the four cod assessments)
3. In RStudio, open the master script and click **Source** to run the complete file. There is **no need to edit paths or set a working directory**. Alternatively, from an R session started in the extracted folder:

```r
source("Mesopredator_MASTER_v1_4_2012_2025_all_cormorants.R")
```

Or run from a terminal in the extracted folder:

```sh
Rscript Mesopredator_MASTER_v1_4_2012_2025_all_cormorants.R
```

The script automatically detects its own location and keeps all local inputs and generated outputs inside the extracted package. This also works if you source the script from a different working directory. You may also run individual lines using **Cmd+Enter** in RStudio, provided the master script is saved and active in the editor (or the extracted package is your working directory).

4. Allow the pipeline to finish. The first stage reads metadata and species data, fits annual diet models, then downloads ICES assessments and runs Monte Carlo simulations and plotting.
5. Look for newly created results in `final result v1.4 2012-2025/`:
   - `Tables/Model_comparison/` (model summaries, fitted diet models, QC)
   - `Tables/ICES/` (cod components and assessment data)
   - `Tables/Monte_Carlo/` (simulated outputs and QC)
   - `Figures/` (model plots)
   - `Documentation/model_run_manifest.csv` (run settings)

## External dependencies and reproducibility
- The script downloads *other fish-stock assessments* dynamically through `icesSAG`. The four bundled cod CSVs are only fallback files for the hard-coded cod assessment keys when their API calls fail. **An entirely offline run is not supported.**
- ICES data and API responses can change; compare the generated ICES provenance files with the manuscript's reported assessment versions.
- The current script is not locked to a package-version environment (no renv.lock). R/package versions may affect numerical output. Record `sessionInfo()` after a successful run.
- The script creates new outputs; it does not require previously generated model `.rds` files because it creates them earlier in the same run.
- Cormorant diet estimates for years before observed samples are model-based assumptions (the script holds pre-observation predictions at the earliest observed fitted value), not direct observations.

## Scientific caveats
- Juvenile cod age-specific mass is based on pooled DATRAS observations (see `Documentation/Evidence/DATRAS_age_weight_summary.csv`).
- Cohort survival uses a historical literature-informed total mortality proxy for young coastal cod, `S = exp(-0.55)`; it is not measured annually for the full 2012–2025 study region. Check exact age interpretation in the primary paper before submission.
- RRA-to-biomass and spatial stock scaling are model assumptions. Interpret outputs as predation-pressure indices, not measured mortality rates.
- The accompanying Word document explains the original-versus-master methodological changes and provides references.

## Checks and troubleshooting
- `Missing metadata/species/wrasse file`: verify the exact path, accents and file names.
- `No ICES assessment found`: confirm internet and `icesSAG` availability. Cod assessment fallback CSVs do not cover the other species.
- `Cannot retrieve ICES assessment`: check the four fallback CSVs and assessment keys.
- `Install patchwork`: install the plotting package and rerun.
- After a successful run, use `sessionInfo()` and keep the run manifest and ICES provenance with shared results.

## Verification status
The package's file paths and script references were checked statically against the source code. **The entire R pipeline has not been executed or numerically validated in this environment**, so this is a run-ready *handoff structure*, not a guarantee of a successful reproducible run on every machine.
