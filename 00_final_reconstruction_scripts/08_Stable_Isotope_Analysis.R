# PORTABLE STANDALONE RECONSTRUCTION. Original calculations retained.
# ============================================================
# 08_Stable_Isotope_Analysis.R
# Mesopredator_check_2: final reconstructed SIA biplot
# Reuses audited SIA reconstruction 31 (no unverified reprocessing)
# Original template: SIA.BIPLOT.ENVFIT.Rmd / sia_biplot_final.png
# ============================================================

required <- c('ggplot2','dplyr')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Missing packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(ggplot2);library(dplyr)})

# Portable bundle root. No file lookups outside this bundle.
wd <- normalizePath(getwd(), mustWork=TRUE)
if (basename(wd) == '00_final_reconstruction_scripts') {
  BASE <- wd
} else if (dir.exists(file.path(wd, '00_final_reconstruction_scripts'))) {
  BASE <- normalizePath(file.path(wd, '00_final_reconstruction_scripts'))
} else {
  stop('Set working directory to the bundle itself or its immediate parent')
}
ROOT <- BASE
FOLDER <- basename(BASE)
SCRIPT_FOLDER <- FOLDER
SOURCE <- file.path(BASE,'processed_data','reference_sia','31_full_reconstruction')
OUT <- file.path(BASE,'processed_data','08_Stable_Isotope_Analysis')
FIG <- file.path(OUT,'figures'); TAB <- file.path(OUT,'tables')
VAL <- file.path(BASE,'validation','08_Stable_Isotope_Analysis')
for(d in c(FIG,TAB,VAL)) dir.create(d,recursive=TRUE,showWarnings=FALSE)

inputs <- c(
 predator=file.path(SOURCE,'31_SIA_predator_samples_FINAL.csv'),
 reference=file.path(SOURCE,'31_SIA_reference_individuals.csv'),
 reference_summary=file.path(SOURCE,'31_SIA_reference_summary_FINAL.csv'),
 predator_summary=file.path(SOURCE,'31_SIA_predator_summary_FINAL.csv')
)
if(any(!file.exists(inputs))) stop('Missing audited SIA input(s): ',paste(inputs[!file.exists(inputs)],collapse='; '),
 '\nRe-run the previously validated 31_SIA_full_reconstruction.R only if source is missing.')

read_checked <- function(path,columns) {
 d <- read.csv(path,stringsAsFactors=FALSE,check.names=FALSE)
 if(!all(columns %in% names(d))) stop('Unexpected columns in ',path,': ',paste(setdiff(columns,names(d)),collapse=', '))
 d
}

pred <- read_checked(inputs[['predator']],c('sample_name','species_label','ocean','X13C_corr','X15N'))
refs <- read_checked(inputs[['reference']],c('VitenskapligNavn','d13C','d15N'))
ref_summary <- read_checked(inputs[['reference_summary']],c('VitenskapligNavn','mean_d13C','sd_d13C','mean_d15N','sd_d15N','n'))
pred_summary <- read_checked(inputs[['predator_summary']],c('species_label','mean_d13C','sd_d13C','mean_d15N','sd_d15N','n'))

# Scientific and file-integrity checks: do not silently drop observations.
if(nrow(pred)!=40L) stop('Expected 40 predator SIA measurements; found ',nrow(pred))
counts <- table(factor(pred$species_label,levels=c('Great cormorant','Harbour seal')))
if(!identical(as.integer(counts),c(18L,22L))) stop('Predator group sizes differ from validated analysis: ',paste(counts,collapse='/'))
if(anyDuplicated(pred$sample_name)) stop('Duplicate SIA sample names')
if(!setequal(pred_summary$species_label,names(counts))) stop('Predator summary labels differ')
if(!identical(sum(as.integer(pred_summary$n)),40L)) stop('Predator summary total differs')
if(!all(c('Clupea harengus','Gadus morhua','Mytilus edulis') %in% ref_summary$VitenskapligNavn)) stop('Expected reference species missing')
if(any(ref_summary$n<10)) stop('Reference display selection (n >= 10) differs')
if(any(ref_summary$VitenskapligNavn=='Pandalus borealis')) stop('Pandalus should not appear in final reference summary')
for(x in list(pred$X13C_corr,pred$X15N,refs$d13C,refs$d15N,
              ref_summary$mean_d13C,ref_summary$mean_d15N,
              pred_summary$mean_d13C,pred_summary$mean_d15N)) {
 if(any(!is.finite(x))) stop('Non-finite isotope measurements detected')
}
# Tissue identity follows Tarjei's final original SIA biplot:
# Oslofjord/oslofjorden = muscle; Skagerrak/SLU = outer stomach wall.
# Use metadata, never isotope positions, to classify tissue.
pred$ocean_clean <- tolower(trimws(as.character(pred$ocean)))
if (anyNA(pred$ocean_clean) || any(!nzchar(pred$ocean_clean))) {
 stop('Missing ocean metadata: cannot classify tissue type safely')
}
pred$tissue_type <- ifelse(
 pred$ocean_clean %in% c('oslofjord','oslofjorden'),
 'Muscle tissue (Oslofjord)',
 'Outer stomach wall (Skagerrak/SLU)'
)
# Follow original source classification: any non-Oslofjord origin is outer stomach wall.
# Preserve origin values in an audit for inspection.
pred$tissue_type <- factor(pred$tissue_type, levels=c(
 'Muscle tissue (Oslofjord)', 'Outer stomach wall (Skagerrak/SLU)'))
cat('Predator tissue type by species:\n')
print(table(pred$species_label,pred$tissue_type))
write.csv(data.frame(sample_name=pred$sample_name, species=pred$species_label,
 ocean=pred$ocean, tissue_type=pred$tissue_type),
 file.path(VAL,'08_TISSUE_TYPE_AUDIT.csv'),row.names=FALSE)

# Original plot shows reference individuals excluding Pandalus.
refs <- refs[refs$VitenskapligNavn!='Pandalus borealis',,drop=FALSE]
# Verify the archived summaries are consistent with individuals.
check_means <- function(raw,group,x,y,summary_df,threshold=1e-5) {
 actual <- aggregate(raw[c(x,y)],by=list(group=raw[[group]]),FUN=mean)
 for(i in seq_len(nrow(summary_df))) {
  name <- if('species_label' %in% names(summary_df)) summary_df$species_label[i] else summary_df$VitenskapligNavn[i]
  row <- actual[actual$group==name,,drop=FALSE]
  if(nrow(row)!=1L) stop('Summary species missing from individuals: ',name)
  if(abs(row[[x]]-summary_df$mean_d13C[i])>threshold || abs(row[[y]]-summary_df$mean_d15N[i])>threshold)
   stop('Mean mismatch between summary and individuals for ',name)
 }
}
check_means(pred,'species_label','X13C_corr','X15N',pred_summary)
check_means(refs,'VitenskapligNavn','d13C','d15N',ref_summary)
cat('Validated SIA: 40 predators (18 cormorant, 22 seal); ',nrow(refs),' reference observations\n',sep='')

pred_colours <- c('Great cormorant'='#d6604d','Harbour seal'='#2166ac')
ref_colours <- c('Clupea harengus'='#4dac26','Gadus morhua'='#8c510a','Mytilus edulis'='#762a83')
colours <- c(pred_colours,ref_colours)

p <- ggplot() +
 geom_point(data=refs,aes(x=d13C,y=d15N,colour=VitenskapligNavn),alpha=.2,size=1.5,shape=16) +
 # Horizontal and vertical SD segments, equivalent to original error bars.
 geom_segment(data=ref_summary,aes(x=mean_d13C-sd_d13C,xend=mean_d13C+sd_d13C,
                                 y=mean_d15N,yend=mean_d15N,colour=VitenskapligNavn),linewidth=.6) +
 geom_segment(data=ref_summary,aes(x=mean_d13C,xend=mean_d13C,
                                 y=mean_d15N-sd_d15N,yend=mean_d15N+sd_d15N,
                                 colour=VitenskapligNavn),linewidth=.6) +
 geom_point(data=ref_summary,aes(x=mean_d13C,y=mean_d15N,colour=VitenskapligNavn),shape=18,size=5) +
 geom_text(data=ref_summary,aes(x=mean_d13C,y=mean_d15N-sd_d15N-.5,label=VitenskapligNavn),
           size=2.8,fontface='italic',colour='grey40',hjust=.5) +
 geom_point(data=pred,aes(x=X13C_corr,y=X15N,colour=species_label,shape=tissue_type),alpha=.6,size=2.5) +
 geom_segment(data=pred_summary,aes(x=mean_d13C-sd_d13C,xend=mean_d13C+sd_d13C,
                                  y=mean_d15N,yend=mean_d15N,colour=species_label),linewidth=.8) +
 geom_segment(data=pred_summary,aes(x=mean_d13C,xend=mean_d13C,
                                  y=mean_d15N-sd_d15N,yend=mean_d15N+sd_d15N,
                                  colour=species_label),linewidth=.8) +
 geom_point(data=pred_summary,aes(x=mean_d13C,y=mean_d15N,colour=species_label),shape=16,size=5) +
 geom_text(data=pred_summary,aes(x=mean_d13C,y=mean_d15N+sd_d15N+.5,label=paste0('n=',n)),
           size=2.8,colour='grey30') +
 scale_colour_manual(values=colours,name=NULL) +
 scale_shape_manual(values=c('Muscle tissue (Oslofjord)'=16,
                             'Outer stomach wall (Skagerrak/SLU)'=17),name='Tissue type') +
 guides(colour=guide_legend(override.aes=list(size=3,alpha=1,shape=16),order=1),
        shape=guide_legend(override.aes=list(size=3,colour='grey30',alpha=1),order=2)) +
 labs(x=expression(paste(delta^{13},'C (‰, lipid-corrected)')),
      y=expression(paste(delta^{15},'N (‰)'))) +
 theme_classic(base_size=11) +
 theme(legend.position='right',legend.text=element_text(face='italic',size=9),
       panel.grid.major=element_line(colour='grey90',linewidth=.3))

figure <- file.path(FIG,'sia_biplot_final.png')
ggsave(figure,p,width=10,height=7,dpi=300,bg='white')
write.csv(pred_summary,file.path(TAB,'08_SIA_predator_summary.csv'),row.names=FALSE)
write.csv(ref_summary,file.path(TAB,'08_SIA_reference_summary.csv'),row.names=FALSE)
write.csv(data.frame(file=unname(inputs),md5=unname(tools::md5sum(inputs))),
          file.path(VAL,'08_INPUT_PROVENANCE.csv'),row.names=FALSE)
writeLines(c(
 'Original source: SIA.BIPLOT.ENVFIT.Rmd, Tarjei Rostum; reconstructed in 31_SIA_full_reconstruction.R.',
 'This script reuses audited individual SIA observations and reference summaries; no correction is recomputed.',
 'Original source applies Kiljunen lipid correction to Oslofjord/ocean-labelled predator samples; check underlying script 31.',
 'Reference records: Oslofjord/Skagerrak selection from WaterRegistrationExport, muscle plus Mytilus BB.',
 'Original final plot uses reference taxa n >= 10 and excludes Pandalus borealis.',
 'Points show individual observations; predator tissue shape from archived ocean metadata (circle Oslofjord muscle; triangle Skagerrak/SLU outer stomach wall).',
 'Mean points remain pooled by predator as in original; means +/- 1 SD; n shown for predators.',
 'A separate envfit analysis exists in earlier reconstructed script 32; not silently integrated into main biplot.',
 'Tissue shape and legend reproduced from the final tissue-type biplot in Tarjei original Rmd; visual comparison still required.'
),file.path(VAL,'08_METHOD_NOTES.txt'))
capture.output(sessionInfo(),file=file.path(VAL,'08_SESSION_INFO.txt'))
cat('\n07 COMPLETE\nFigure: ',figure,'\nTables: ',TAB,'\nValidation: ',VAL,'\n',sep='')
