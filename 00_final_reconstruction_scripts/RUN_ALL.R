# RUN_ALL.R — run ONLY from the standalone bundle, no external project files.
# Run in a fresh R session: setwd('/your/path/00_final_reconstruction_scripts'); source('RUN_ALL.R')
if(basename(normalizePath(getwd())) != '00_final_reconstruction_scripts') stop('Set working directory to the standalone bundle before calling RUN_ALL.R')
BASE <- normalizePath(getwd())
if(!all(file.exists(file.path(BASE,'rawfiles',c('ps_with_host.rds','ps_decon_mifish.rds','ps_filtered_mifish.rds','ps_coi_metazoa.rds'))))) stop('Raw files missing; run the external prepare_rawfiles.R first.')
steps <- c(
  '01_Sequencing_QC.R',
  '02_Rarefaction_Threshold_Test.R',
  '03_Phyloseq_Preprocessing.R',
  '04_PERMANOVA_PCoA.R',
  '04B_PERMANOVA_PCoA.R',
  '05_Diet_Composition_FOO_RRA.R',
  '06_GLLVM.R',
  '07_Trophic_Network.R',
  '07B_Trophic_Network_Sensitivity_Jackknife.R',
  '08_Stable_Isotope_Analysis.R',
  '09_Functional_Group_Composition.R',
  '10_Manuscript_Tables.R'
)
# -------------------------------------------------------------
# RUN ALL RECONSTRUCTION SCRIPTS
# -------------------------------------------------------------

for (script in steps) {
  
  cat(
    "\n======= RUNNING: ",
    script,
    " =======\n",
    sep = ""
  )
  
  # Execute scripts in the global environment,
  # matching successful standalone execution.
  sys.source(
    file.path(BASE, script),
    envir = .GlobalEnv
  )
  
  cat(
    "======= COMPLETED: ",
    script,
    " =======\n",
    sep = ""
  )
}
source(file.path(BASE,'COLLECT_MANUSCRIPT_OUTPUTS.R'),local=new.env(parent=globalenv()))
cat('\nBUNDLE COMPLETE. Review figures/, tables/, processed_data/, validation/.\n')
