# PORTABLE STANDALONE RECONSTRUCTION. Original calculations retained.
# =============================================================================
# 01_Sequencing_QC.R — Tarjei MiFish quality-control figure reconstruction
# Five outputs: scrub_plate_heatmap.png, scrub_contamination.png,
# read_depth.png, host_dna.png, host_dna_by_processor.png
# =============================================================================
# Original sources:
#  setup.n.scrub.2026.04.20.Rmd (SCRuB model, contamination, plate heatmap)
#  updated.masterchunk.sample.quality.Rmd (host DNA, read depth, processor)
# IMPORTANT: Original SCRuB model fit objects were not archived. This script
# refits SCRuB with the original seed and documented plate/control metadata;
# fitted p values MAY differ from the original. It never estimates p from PNGs.
# No synthetic observations, no forced model agreement, no edits to scripts 01-08.
# =============================================================================

required <- c('phyloseq','ggplot2','dplyr','tidyr','scales','gridExtra')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install required packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({
  library(phyloseq); library(ggplot2); library(dplyr); library(tidyr)
})

# Locate project and existing original RDS objects
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
OUT <- file.path(BASE,'processed_data','01_Sequencing_QC')
FIG <- file.path(OUT,'figures'); TAB <- file.path(OUT,'tables')
VAL <- file.path(BASE,'validation','01_Sequencing_QC')
for(d in c(FIG,TAB,VAL)) dir.create(d,recursive=TRUE,showWarnings=FALSE)

# Configure an explicit raw-data directory if needed, before sourcing:
# QC_INPUT_DIR <- '/path/to/original/MiFish/files'
input_dirs <- file.path(BASE, 'rawfiles')
locate <- function(fname, required=TRUE){
  p <- file.path(input_dirs,fname); p <- p[file.exists(p)]
  if(length(p)) return(normalizePath(p[[1]]))
  if(required) stop('Required original MiFish input not found: ',fname,'\nSearched:\n',paste(input_dirs,collapse='\n'),
                    '\nRun the external prepare_rawfiles.R to populate ../../data/raw/.')
  NULL
}
path_with_host <- locate('ps_with_host.rds')
path_decon <- locate('ps_decon_mifish.rds')
ps_wh_raw <- readRDS(path_with_host)
ps_de_raw <- readRDS(path_decon)
stopifnot(inherits(ps_wh_raw,'phyloseq'),inherits(ps_de_raw,'phyloseq'))
# The pre-SCRuB filtered object is not archived as an RDS. Reconstruct it from
# the saved raw ps_with_host by removing the four host genera, as in Tarjei.

# -------------------------------------------------------------
# Remove host genera from the original phyloseq object
# -------------------------------------------------------------

host_genera <- c(
  "Phoca",
  "Homo",
  "Phalacrocorax",
  "Halichoerus"
)

tax_df <- as.data.frame(
  phyloseq::tax_table(ps_wh_raw),
  check.names = FALSE
)

if (!"Genus" %in% colnames(tax_df)) {
  stop("Original ps_with_host lacks taxonomy Genus required for host removal.")
}

# Apply the same host-genus filtering as the original analysis.
# Explicit indexing avoids subset_taxa() evaluation issues.

keep_taxa <- is.na(tax_df$Genus) |
  !(as.character(tax_df$Genus) %in% host_genera)

ps_filtered_raw <- phyloseq::prune_taxa(
  rownames(tax_df)[keep_taxa],
  ps_wh_raw
)

# Remove taxa with zero reads after filtering
ps_filtered_raw <- phyloseq::prune_taxa(
  phyloseq::taxa_sums(ps_filtered_raw) > 0,
  ps_filtered_raw
)


# Work with the same OS19 omission documented in the original master chunk.
# OS19 is excluded for final host/read figures. It is NOT silently removed from
# the SCRuB fitting cohort, which is treated separately below.
os19 <- c('d-4-29','d-4-30','d-4-31','d-4-32')
remove_os19 <- function(ps){ prune_samples(!(sample_names(ps) %in% os19),ps) }
true_only <- function(ps){
  md <- as.data.frame(sample_data(ps))
  if(!'Sample_or_Control' %in% names(md)) stop('Missing Sample_or_Control metadata')
  prune_samples(rownames(md)[as.character(md$Sample_or_Control)=='True_Sample'],ps)
}
as_samples <- function(ps){m <- as(otu_table(ps),'matrix');if(taxa_are_rows(ps))m<-t(m);m}
metadata <- function(ps){ x<-as.data.frame(sample_data(ps),check.names=FALSE);x$library_name<-rownames(x);x }
get_species <- function(x){
  x<-as.character(x); x[x=='top_cormorant']<-'cormorant'
  ifelse(x=='cormorant','Great cormorant',ifelse(x=='harbor_seal','Harbour seal',NA_character_))
}
cols_species <- c('Great cormorant'='#d6604d','Harbour seal'='#2166ac')

# ---- A. Host-DNA figure -------------------------------------------------------
ps_host <- true_only(remove_os19(ps_wh_raw))
hostY <- as_samples(ps_host)
tax <- as.data.frame(tax_table(ps_host),check.names=FALSE)
tax <- tax[colnames(hostY),,drop=FALSE]
if(!identical(rownames(tax),colnames(hostY))) stop('Host taxonomy does not align with OTU table')
read_total <- rowSums(hostY)
if(any(read_total<=0)) stop('Zero-read true samples in ps_with_host: inspect data.')
H <- vapply(host_genera,function(genus) rowSums(hostY[,!is.na(tax$Genus) & as.character(tax$Genus)==genus,drop=FALSE]),numeric(nrow(hostY)))
rownames(H) <- rownames(hostY)
H <- sweep(H,1,read_total,'/')*100
hdf <- data.frame(library_name=rownames(H),H,check.names=FALSE) |>
  left_join(metadata(ps_host),by='library_name') |>
  mutate(species_label=get_species(species),
         Non_host=100-rowSums(across(all_of(host_genera))),
         total_host=rowSums(across(all_of(host_genera))))
if(anyNA(hdf$species_label)) stop('Unrecognized species in original host DNA metadata')
if(any(hdf$Non_host < -1e-6)) stop('Host reads exceed total reads')
write.csv(hdf,file.path(TAB,'01_HOST_DNA_BY_SAMPLE.csv'),row.names=FALSE)

# Exactly one bar per library (not sample_name, which may be duplicated).
ord <- hdf |> arrange(species_label,desc(total_host),library_name)
ord$plot_id <- paste(ord$species_label,ord$library_name,sep='__')
long_host <- ord |>
  select(plot_id,species_label,all_of(host_genera),Non_host) |>
  pivot_longer(cols=c(all_of(host_genera),'Non_host'),names_to='genus',values_to='pct')
long_host$genus <- factor(long_host$genus,levels=c('Non_host','Homo','Halichoerus','Phoca','Phalacrocorax'))
# Keep highest total host-DNA percentage on the LEFT, as in Tarjei's figure.
# ord is already sorted with desc(total_host) within each predator.
# Reversing its factor levels would invert the original visual order.
long_host$plot_id <- factor(long_host$plot_id,levels=unique(ord$plot_id))
col_host <- c('Non_host'='grey85','Homo'='#fdae61','Halichoerus'='#4dac26',
              'Phoca'='#2166ac','Phalacrocorax'='#d6604d')
host_labels <- c('Non-host','Homo sapiens','Halichoerus grypus',
                 'Phoca vitulina','Phalacrocorax carbo')
P_host_top <- ggplot(long_host,aes(x=plot_id,y=pct,fill=genus)) +
  geom_col(width=1) + facet_wrap(~species_label,scales='free_x') +
  scale_fill_manual(values=col_host,breaks=names(col_host),labels=host_labels,name='Taxon') +
  scale_y_continuous(limits=c(0,100),expand=expansion(mult=c(0,0))) +
  labs(x=NULL,y='Proportion of reads (%)') + theme_classic() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        strip.background=element_blank(),strip.text=element_text(face='italic',size=11),
        legend.position='right')
P_host_bottom <- ggplot(hdf,aes(x=species_label,y=total_host+0.01,fill=species_label,colour=species_label)) +
  geom_boxplot(width=.3,outlier.shape=NA,alpha=.4) +
  geom_jitter(width=.1,size=1.5,alpha=.6) +
  scale_y_log10(breaks=c(.01,.1,1,10,100),labels=c('0','0.1','1','10','100')) +
  scale_fill_manual(values=cols_species) + scale_colour_manual(values=cols_species) +
  labs(x=NULL,y='Total host DNA (%, log scale)')+theme_classic()+
  theme(legend.position='none',axis.text.x=element_text(face='italic'))
G_host <- gridExtra::arrangeGrob(P_host_top,P_host_bottom,nrow=2,heights=c(2,1))
ggsave(file.path(FIG,'host_dna.png'),plot=G_host,width=10,height=8,dpi=300,bg='white')

# ---- B. Host DNA by processor ------------------------------------------------
if(!'person_emptied' %in% names(hdf)) stop('person_emptied field missing; cannot reproduce processor comparison')
hdf <- hdf |> mutate(processor=ifelse(as.character(person_emptied)=='SLU','SLU','Self-processed'))
write.csv(hdf,file.path(TAB,'01_HOST_DNA_BY_PROCESSOR.csv'),row.names=FALSE)
P_processor <- ggplot(hdf,aes(x=processor,y=total_host+.01,fill=species_label,colour=species_label))+
  geom_boxplot(width=.4,outlier.shape=NA,alpha=.5)+
  geom_jitter(width=.1,size=1.5,alpha=.5)+
  facet_wrap(~species_label,scales='free_x')+
  scale_y_log10(breaks=c(.01,.1,1,10,100),labels=c('0','0.1','1','10','100'))+
  scale_fill_manual(values=cols_species)+scale_colour_manual(values=cols_species)+
  labs(x='Processor',y='Total host DNA (%, log scale)')+theme_classic()+
  theme(legend.position='none',strip.background=element_blank(),
        strip.text=element_text(face='italic',size=11))
ggsave(file.path(FIG,'host_dna_by_processor.png'),plot=P_processor,width=8,height=5,dpi=300,bg='white')

# ---- C. Read-depth comparison ------------------------------------------------
# Same original stage definitions: before host removal (ps_with_host), after
# host removal (ps_filtered), after SCRuB (original saved ps_decon).
ps_fi <- true_only(remove_os19(ps_filtered_raw))
ps_de <- true_only(remove_os19(ps_de_raw))
stages <- list('Before host removal'=ps_host,
               'After host removal'=ps_fi,'After SCRuB'=ps_de)
reads_df <- bind_rows(lapply(names(stages),function(stage){
  obj<-stages[[stage]]; m<-as_samples(obj); mt<-metadata(obj)
  data.frame(library_name=rownames(m),reads=rowSums(m),stage=stage) |>
    left_join(mt[,c('library_name','species'),drop=FALSE],by='library_name') |>
    mutate(species_label=get_species(species))
}))
if(anyNA(reads_df$species_label)) stop('Unrecognized predator species in read-depth dataset')
reads_df$stage <- factor(reads_df$stage,levels=names(stages))
ref_order <- names(sort(rowSums(as_samples(ps_host)),decreasing=TRUE))
reads_df$library_name <- factor(reads_df$library_name,levels=unique(c(ref_order,as.character(reads_df$library_name))))
write.csv(reads_df,file.path(TAB,'01_READ_DEPTH_BY_SAMPLE_AND_STAGE.csv'),row.names=FALSE)
P_reads <- ggplot(reads_df,aes(x=library_name,y=reads,fill=species_label)) +
  geom_col(width=1)+facet_grid(stage~species_label,scales='free',space='free_x')+
  scale_fill_manual(values=cols_species)+scale_y_continuous(labels=scales::comma)+
  labs(x=NULL,y='Number of reads')+theme_classic()+
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position='none',
        strip.background=element_blank(),strip.text.x=element_text(size=11),
        strip.text.y=element_text(angle=0,size=10))
ggsave(file.path(FIG,'read_depth.png'),plot=P_reads,width=10,height=8,dpi=300,bg='white')

# ---- D. Refit SCRuB from original plate & control labels --------------------
# SCRuB is NOT run from the pruned analysis phyloseq: positive and negative
# controls are needed. Refitting can produce different values from archived PNG.
if(!requireNamespace('SCRuB',quietly=TRUE)) {
  warning('SCRuB package is not installed: host/read figures generated, SCRuB figures skipped. ',
          'Install SCRuB from Shenhav-and-Korem-labs/SCRuB to complete all five.')
} else {
  set.seed(1)
  ps_sc <- ps_filtered_raw
  sc_md <- metadata(ps_sc)
  needed <- c('plate_pcr2','plate_well_pcr2','Sample_or_Control','source','species')
  if(!all(needed %in% names(sc_md))) stop('SCRuB metadata lacks: ',paste(setdiff(needed,names(sc_md)),collapse=', '))
  # Record every sample actually used in fit; no arbitrary inferred plate wells.
  contam_list <- list(); plate_audit <- list()
  for(k in 1:2){
    plate_id <- paste0('plate_',k)
    plate_ps <- prune_samples(sc_md$library_name[as.character(sc_md$plate_pcr2)==plate_id],ps_sc)
    dat <- as_samples(plate_ps)
    dat <- dat[rowSums(dat)>0,,drop=FALSE]
    m <- metadata(plate_ps)
    m <- m[match(rownames(dat),m$library_name),,drop=FALSE]
    if(!identical(rownames(dat),m$library_name))stop('SCRuB metadata alignment error, ',plate_id)
    flags <- as.character(m$Sample_or_Control)
    if(!any(flags=='is.neg'))stop('No negative controls for SCRuB in ',plate_id)
    model_meta <- data.frame(is_control=flags=='is.neg',
                             sample_type=as.character(m$source),
                             sample_well=as.character(m$plate_well_pcr2),
                             row.names=m$library_name)
    if(anyNA(model_meta$sample_well) || any(!nzchar(model_meta$sample_well)))
       stop('SCRuB sample well missing in ',plate_id)
    fit <- SCRuB::SCRuB(dat,model_meta,control_order='pcr_neg')
    p <- fit$p
    if(is.null(names(p)) || anyDuplicated(names(p)))
       stop('SCRuB p returned without unique sample IDs; inspect fitted object, ',plate_id)
    if(any(!names(p) %in% rownames(model_meta)))stop('SCRuB returned unknown sample IDs')
    contam_list[[k]] <- data.frame(library_name=names(p),
                  contamination=as.numeric(1-p),plate=paste('Plate',k))
    plate_audit[[k]] <- data.frame(library_name=m$library_name,
                        plate=paste('Plate',k),well=model_meta$sample_well,
                        Sample_or_Control=flags,stringsAsFactors=FALSE)
  }
  C <- bind_rows(contam_list) |>
    left_join(sc_md[,c('library_name','species','Sample_or_Control'),drop=FALSE],by='library_name') |>
    mutate(species_label=case_when(
      Sample_or_Control=='is.neg'~'Negative control',
      Sample_or_Control=='pos_control'~'Positive control',
      TRUE~get_species(species)))
  if(anyNA(C$species_label)) warning('Some SCRuB sample labels not assigned to predators or controls.')
  write.csv(C,file.path(TAB,'01_SCRUB_ESTIMATED_CONTAMINATION.csv'),row.names=FALSE)
  write.csv(bind_rows(plate_audit),file.path(VAL,'01_SCRUB_PLATE_FIT_MEMBERSHIP.csv'),row.names=FALSE)
  Cplot <- C |> filter(!is.na(species_label)) |>
    mutate(contamination_log=pmax(contamination,0)+.0001)
  P_scrub <- ggplot(Cplot,aes(x=species_label,y=contamination_log,colour=species_label))+
    geom_boxplot(width=.4,outlier.shape=NA)+
    geom_jitter(width=.1,size=2,alpha=.6)+
    facet_wrap(~plate,scales='free_x')+
    scale_colour_manual(values=c(cols_species,'Positive control'='darkgreen','Negative control'='black'))+
    scale_y_log10(breaks=c(.0001,.001,.01,.1),labels=c('0','0.1%','1%','10%'))+
    labs(x=NULL,y='Estimated contamination level (%, log scale)')+theme_classic()+
    theme(legend.position='none',strip.background=element_blank(),
          strip.text=element_text(size=11,face='bold'),
          axis.text.x=element_text(size=9,angle=20,hjust=1))
  ggsave(file.path(FIG,'scrub_contamination.png'),plot=P_scrub,width=10,height=5,dpi=300,bg='white')

  # Rebuild 8 x 12 well grids from the original plate_well_pcr2 field.
  # Original Rmd references plate_layout but does not show its construction.
  parse_well <- function(s){
    s<-toupper(trimws(as.character(s)))
    s<-gsub('^[^A-H0-9]*','',s)
    row<-sub('^([A-H])([0-9]{1,2}).*$','\\1',s)
    col<-suppressWarnings(as.integer(sub('^([A-H])([0-9]{1,2}).*$','\\2',s)))
    data.frame(row=row,col=col)
  }
  membership <- bind_rows(plate_audit)
  w <- parse_well(membership$well)
  membership$row<-w$row;membership$col<-w$col
  if(any(!membership$row %in% LETTERS[1:8] | is.na(membership$col) |
         membership$col<1 | membership$col>12))
    stop('Unrecognized plate well format. See 01_SCRUB_PLATE_FIT_MEMBERSHIP.csv')
  if(anyDuplicated(paste(membership$plate,membership$row,membership$col)))
    stop('Duplicated PCR2 plate well assignment; inspect 01_SCRUB_PLATE_FIT_MEMBERSHIP.csv')
  membership <- membership |>
    left_join(C[,c('library_name','contamination')],by='library_name') |>
    mutate(sample_type=case_when(Sample_or_Control=='is.neg'~'Negative',
         Sample_or_Control=='pos_control'~'Positive',TRUE~'Sample'))
  full_grid <- tidyr::expand_grid(plate=c('Plate 1','Plate 2'),
                                  row=LETTERS[1:8],col=1:12) |>
    left_join(membership[,c('plate','row','col','sample_type','contamination','library_name')],
              by=c('plate','row','col')) |>
    mutate(sample_type=coalesce(sample_type,'Empty'),
           row=factor(row,levels=rev(LETTERS[1:8])),
           col=factor(col,levels=as.character(1:12)))
  write.csv(full_grid,file.path(TAB,'01_SCRUB_PLATE_LAYOUT.csv'),row.names=FALSE)
  P_plate <- ggplot(full_grid,aes(x=col,y=row))+
    geom_tile(data=subset(full_grid,sample_type=='Sample'),aes(fill=contamination),
              colour='grey60',linewidth=.4)+
    geom_tile(data=subset(full_grid,sample_type=='Negative'),fill='#FF0000',colour='grey60',linewidth=.4)+
    geom_tile(data=subset(full_grid,sample_type=='Positive'),fill='#4dac26',colour='grey60',linewidth=.4)+
    geom_tile(data=subset(full_grid,sample_type=='Empty'),fill='white',colour='grey60',linewidth=.4)+
    geom_text(data=subset(full_grid,sample_type=='Negative'),aes(label='NEG'),
              size=2,colour='white',fontface='bold')+
    geom_text(data=subset(full_grid,sample_type=='Positive'),aes(label='POS'),
              size=2,colour='white',fontface='bold')+
    facet_wrap(~plate,ncol=1)+
    scale_fill_gradientn(colours=c('white','white','#fdae61','#d6604d'),
      values=scales::rescale(c(0,.005,.05,.11)),limits=c(0,.11),
      oob=scales::squish,na.value='white',name='Contamination\nlevel',
      labels=scales::percent_format(accuracy=1))+
    scale_x_discrete(position='top')+coord_fixed()+labs(x=NULL,y=NULL)+theme_classic()+
    theme(strip.background=element_blank(),strip.text=element_text(size=11,face='bold'),
          axis.text.x=element_text(size=9,face='bold'),
          axis.text.y=element_text(size=9,face='bold'),
          axis.line=element_blank(),axis.ticks=element_blank(),
          legend.position='right',panel.border=element_rect(colour='grey60',fill=NA))
  ggsave(file.path(FIG,'scrub_plate_heatmap.png'),plot=P_plate,width=10,height=8,dpi=300,bg='white')
}

# ---- Provenance / known reproduction limitations -----------------------------
infiles <- c(path_with_host,path_decon)
write.csv(data.frame(file=infiles,md5=unname(tools::md5sum(infiles))),
          file.path(VAL,'01_INPUT_PROVENANCE.csv'),row.names=FALSE)
writeLines(c(
 'SCRIPT 09: Original MiFish sequencing / contamination QC figure reconstruction.',
 'Original code sources: setup.n.scrub.2026.04.20.Rmd and updated.masterchunk.sample.quality.Rmd.',
 'Figures: scrub_plate_heatmap.png; scrub_contamination.png; read_depth.png;',
 '         host_dna.png; host_dna_by_processor.png.',
 'Original ps_with_host.rds and ps_decon_mifish.rds are inputs.',
 'The original ps_filtered intermediate was not saved: reconstructed using host genus exclusions.',
 'OS19 d-4-29 through d-4-32 excluded in host/read figures, per original master chunk.',
 'SCRuB fit uses original per-plate controls and sample wells, seed 1; original SCRuB',
 'parameter objects were not archived, so model p values can differ from Tarjei output.',
 'SCRuB plate layout rebuilt from metadata: original Rmd references undefined plate_layout.',
 'SCRuB fitting cohort is separate from host/read cohorts.',
 'All figures remain candidates until visually compared with the originals.',
 'No original methods or results were changed to force a visual match.'
),file.path(VAL,'01_METHOD_AND_DEVIATION_NOTES.txt'))
capture.output(sessionInfo(),file=file.path(VAL,'01_SESSION_INFO.txt'))
cat('\n09 COMPLETED. Figures located at:\n',FIG,'\n',sep='')
