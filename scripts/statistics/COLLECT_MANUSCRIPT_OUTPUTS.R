# COLLECT_MANUSCRIPT_OUTPUTS.R — gathers only mapped, reconstructed exports.
if(basename(normalizePath(getwd())) != '00_final_reconstruction_scripts') stop('Run from bundle root')
B <- normalizePath(getwd())
dir.create(file.path(B,'figures','main_figures'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(B,'figures','supplementary_figures'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(B,'tables','main_tables'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(B,'tables','supplementary_tables'),recursive=TRUE,showWarnings=FALSE)
figmap <- list(
 c('main_figures','05_Diet_Composition_FOO_RRA/figures/diet_FOO_four_panels.png','Figure1_diet_FOO_four_panels.png'),
 c('main_figures','04_PERMANOVA_PCoA/figures/PCoA_diet_composition_four_panels.png','Figure2_PCoA_diet_composition_four_panels.png'),
 c('main_figures','06_GLLVM/figures/gllvm_coefficients_combined.png','Figure3_GLLVM_coefficients.png'),
 c('main_figures','07_Trophic_Network/figures/trophic_network.png','Figure4_trophic_network.png'),
 c('main_figures','08_Stable_Isotope_Analysis/figures/sia_biplot_final.png','Figure7_sia_biplot.png'),
 c('supplementary_figures','01_Sequencing_QC/figures/host_dna.png','host_dna.png'),
 c('supplementary_figures','01_Sequencing_QC/figures/host_dna_by_processor.png','host_dna_by_processor.png'),
 c('supplementary_figures','01_Sequencing_QC/figures/read_depth.png','read_depth.png'),
 c('supplementary_figures','01_Sequencing_QC/figures/scrub_plate_heatmap.png','scrub_plate_heatmap.png'),
 c('supplementary_figures','01_Sequencing_QC/figures/scrub_contamination.png','scrub_contamination.png'),
 c('supplementary_figures','02_Rarefaction_Threshold_Test/figures/rarefaction_curves.png','rarefaction_curves.png'),
 c('supplementary_figures','02_Rarefaction_Threshold_Test/figures/read_threshold_diagnostics.png','read_threshold_diagnostics.png'),
 c('supplementary_figures','09_Functional_Group_Composition/figures/FG_composition_PHOCA_VITULINA_COI.png','FG_composition_PHOCA_VITULINA_COI.png')
)
# betadisper output filenames may be different; include all distinct betadisper figures.
dispdir <- file.path(B,'processed_data','04_PERMANOVA_PCoA','figures')
if(dir.exists(dispdir)) for(f in list.files(dispdir,pattern='^betadisper_.*\\.png$',full.names=TRUE)) {
   figmap[[length(figmap)+1L]] <- c('supplementary_figures',paste0('04_PERMANOVA_PCoA/figures/',basename(f)),basename(f))
}
fig_status <- lapply(figmap,function(m){
  src <- file.path(B,'processed_data',m[[2]]);dest <- file.path(B,'figures',m[[1]],m[[3]])
  if(file.exists(src)) file.copy(src,dest,overwrite=TRUE)
  data.frame(class=m[[1]],source=m[[2]],figure=m[[3]],found=file.exists(src))
})
write.csv(do.call(rbind,fig_status),file.path(B,'validation','FIGURE_COLLECTION_MANIFEST.csv'),row.names=FALSE)
if(any(!vapply(fig_status,function(z) z$found,logical(1)))) warning('Some figures were not produced. Review ../../docs/validation/FIGURE_COLLECTION_MANIFEST.csv')
# Exact manuscript tables, with filenames retained.
tabdir <- file.path(B,'processed_data','manuscript_tables')
for(tp in c('main','supplementary')) {
 src <- file.path(tabdir,tp)
 if(dir.exists(src)) for(f in list.files(src,full.names=TRUE)) {
  if(file.info(f)$isdir) next
  dest <- file.path(B,'tables',paste0(tp,'_tables'),basename(f))
  file.copy(f,dest,overwrite=TRUE)
 }
}
# Ancillary approved tables are recorded under processed_data, not inserted into manuscript set.
cat('Figure and table collection completed. Check ../../docs/validation/FIGURE_COLLECTION_MANIFEST.csv.\n')

# FINAL_CLEANUP_EXPORT_DEDUP_V1
# Keep original downstream model inputs and reference tables intact.
if (exists("figmap")) for (m in figmap) {
  source_png <- file.path(B,"processed_data",m[[2]])
  final_png <- file.path(B,"figures",m[[1]],m[[3]])
  # 07B checks the 07 network figure: preserve it as a dependency.
  if (grepl("^07_Trophic_Network/",m[[2]])) next
  if (file.exists(source_png) && file.exists(final_png) &&
      identical(unname(tools::md5sum(source_png)),
                unname(tools::md5sum(final_png)))) {
    unlink(source_png)
  }
}
# Final manuscript tables are not read by any later stage (collector runs last).
# Keep all non-manuscript intermediate and jackknife-reference tables.
for (tp in c("main","supplementary")) {
  original_dir <- file.path(B,"processed_data","manuscript_tables",tp)
  final_dir <- file.path(B,"tables",paste0(tp,"_tables"))
  if (!dir.exists(original_dir)) next
  for (p in list.files(original_dir,full.names=TRUE)) {
    if (dir.exists(p)) next
    q <- file.path(final_dir,basename(p))
    if (file.exists(q) && identical(unname(tools::md5sum(p)),
                                     unname(tools::md5sum(q)))) unlink(p)
  }
}
# Do not delete intermediate directories: they contain essential RDS/CSV.
