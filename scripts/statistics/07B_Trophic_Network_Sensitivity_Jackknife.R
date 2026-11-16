# PORTABLE STANDALONE RECONSTRUCTION. Original calculations retained.
# ============================================================
# 07B_Trophic_Network_Sensitivity_Jackknife.R
# Mesopredator_check_2
# Supplementary QA for the approved 07_Trophic_Network.R.
# Does not overwrite the approved Figure 06 or its tables.
# Tests fish-invertebrate co-occurrence thresholds: 5, 10, 15%.
# Predator-fish FOO remains strictly >5%; top-eight selection is
# recalculated on each leave-one-out replicate.
# A co-occurrence edge is NOT evidence of direct predation.
# ============================================================

required <- c('phyloseq','dplyr','ggplot2')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(phyloseq); library(dplyr); library(ggplot2)})

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
INPUT <- file.path(BASE,'processed_data','03_Phyloseq_Preprocessing')
FIG06 <- file.path(BASE,'processed_data','07_Trophic_Network','figures','trophic_network.png')
TAB06 <- file.path(BASE,'processed_data','07_Trophic_Network','tables')
VAL <- file.path(BASE,'validation','07B_Trophic_Network_Sensitivity_Jackknife')
dir.create(VAL,recursive=TRUE,showWarnings=FALSE)
mi_file <- file.path(INPUT,'ps_mifish_global_0_1pct_identified.rds')
co_file <- file.path(INPUT,'ps_coi_global_0_1pct_identified.rds')
ref_fish <- file.path(TAB06,'07_TOP8_FISH_FOO.csv')
ref_pf <- file.path(TAB06,'07_PREDATOR_FISH_EDGES.csv')
ref_fi <- file.path(TAB06,'07_FILTERED_FISH_INVERTEBRATE_EDGES.csv')
required_files <- c(mi_file,co_file,FIG06,ref_fish,ref_pf,ref_fi)
if(!all(file.exists(required_files))) stop('Run approved 06 first. Missing:\n',paste(required_files[!file.exists(required_files)],collapse='\n'))

ps_mf <- readRDS(mi_file)
ps_coi <- readRDS(co_file)
stopifnot(inherits(ps_mf,'phyloseq'),inherits(ps_coi,'phyloseq'))
if(!'LCA_taxon' %in% colnames(as(tax_table(ps_mf),'matrix'))) stop('MiFish LCA_taxon missing')
ps_mf <- subset_taxa(ps_mf,LCA_taxon!='no identification')

extract_data <- function(ps) {
 m <- as(otu_table(ps),'matrix')
 if(taxa_are_rows(ps)) m <- t(m)
 tx <- as.data.frame(tax_table(ps),check.names=FALSE)
 md <- as.data.frame(sample_data(ps),check.names=FALSE)
 stopifnot(all(colnames(m)%in%rownames(tx)),all(rownames(m)%in%rownames(md)))
 tx <- tx[colnames(m),,drop=FALSE]
 md <- md[rownames(m),,drop=FALSE]
 stopifnot(all(c('LCA_taxon','Class') %in% names(tx)),'species'%in%names(md))
 labels <- gsub('_',' ',as.character(tx$LCA_taxon),fixed=TRUE)
 if(anyNA(labels)||any(!nzchar(labels))) stop('Empty taxon label')
 mat <- t(rowsum(t(m),group=labels,reorder=FALSE))
 list(otu=mat,tx=tx,md=md,labels=labels)
}
mi <- extract_data(ps_mf); co <- extract_data(ps_coi)
shared <- intersect(rownames(mi$otu),rownames(co$otu))
if(length(shared)<10) stop('Too few shared samples')
class_by_label <- split(as.character(co$tx$Class),co$labels)
if(any(vapply(class_by_label,function(z)length(unique(z))>1L,logical(1)))) stop('Conflicting COI classes')
invertebrates <- names(class_by_label)[vapply(class_by_label,function(z)identical(unique(z),'Actinopteri')==FALSE,logical(1))]
invertebrates <- intersect(colnames(co$otu),invertebrates)
m <- mi$otu[shared,,drop=FALSE]
v <- co$otu[shared,invertebrates,drop=FALSE]
# Mirror metadata extraction in approved script 06. In particular, avoid
# converting the one-column data.frame directly with as.character().
md_mifish <- as.data.frame(mi$md[shared, , drop=FALSE])
md_coi <- as.data.frame(co$md[shared, , drop=FALSE])
pred <- unname(trimws(as.character(md_mifish[['species']])))
co_pred <- unname(trimws(as.character(md_coi[['species']])))
pred[pred=='top_cormorant'] <- 'cormorant'
co_pred[co_pred=='top_cormorant'] <- 'cormorant'
if(length(pred)!=length(shared) || length(co_pred)!=length(shared)) stop('Predator metadata length mismatch')
metadata_audit <- data.frame(sample_id=shared,
                             MiFish_original=as.character(md_mifish[['species']]),
                             COI_original=as.character(md_coi[['species']]),
                             MiFish_normalized=pred,
                             COI_normalized=co_pred,
                             stringsAsFactors=FALSE)
write.csv(metadata_audit,file.path(VAL,'07B_PREDATOR_METADATA_AUDIT.csv'),row.names=FALSE)
if(anyNA(pred) || anyNA(co_pred) || !identical(pred,co_pred) ||
   !all(pred %in% c('cormorant','harbor_seal'))) {
  print(metadata_audit[is.na(pred)|is.na(co_pred)|pred!=co_pred,,drop=FALSE])
  stop('True predator metadata mismatch; see 07B_PREDATOR_METADATA_AUDIT.csv')
}
cat('Validated matched predator labels: ',paste(names(table(pred)),as.integer(table(pred)),collapse=', '),'\n',sep='')
if(anyNA(m)||anyNA(v)||any(!is.finite(m))||any(!is.finite(v))||any(m<0)||any(v<0)) stop('Invalid abundance matrices')

# Exactly the taxon-shortening rules of approved script 06.
short_fish <- function(x) vapply(x,function(z) {
 if(grepl('Myoxocephalus scorpius Myoxocephalus',z)) return('Myoxocephalus scorpius')
 if(grepl('Pleuronectes platessa Platichthys',z)) return('Pleuronectes platessa')
 if(grepl('Platichthys flesus Pleuronectes platessa Hippoglossoides',z)) return('Platichthys flesus')
 if(grepl('Pollachius virens Pollachius|Pollachius pollachius',z)) return('Pollachius spp.')
 if(grepl('Limanda limanda Hippoglossus',z)) return('Limanda limanda')
 if(grepl('Clupea harengus.*Clupea pallasii|Clupea pallasii.*Clupea harengus',z)) return('Clupea harengus')
 if(grepl('Hyperoplus.*Ammodytes|Ammodytes.*Hyperoplus',z)) return('Ammodytes spp.')
 if(grepl('Sprattus sprattus Clupea|Clupea pallasii Sprattus',z)) return('Sprattus sprattus')
 paste(head(strsplit(z,' ',fixed=TRUE)[[1]],2),collapse=' ')
},character(1),USE.NAMES=FALSE)
short_inv <- function(x) vapply(x,function(z) {
 if(grepl('Evadne',z)) return('Evadne sp.')
 if(grepl('Alitta|Hediste',z)) return('Alitta / Hediste')
 if(grepl('Echinogammarus.*Gammarus',z)) return('Echinogammarus sp.')
 if(grepl('Heteromastus|Capitellidae',z)) return('Heteromastus sp.')
 if(grepl('Palaemon',z)) return('Palaemon adspersus')
 if(grepl('Harmothoe',z)) return('Harmothoe imbricata')
 if(grepl('Liocarcinus',z)) return('Liocarcinus holsatus')
 if(grepl('Penilia',z)) return('Penilia avirostris')
 if(grepl('Acartia',z)) return('Acartia spp.')
 if(grepl('Cyclopoida',z)) return('Cyclopoida sp.')
 paste(head(strsplit(z,' ',fixed=TRUE)[[1]],2),collapse=' ')
},character(1),USE.NAMES=FALSE)
fish_short <- short_fish(colnames(m))
inv_short <- short_inv(colnames(v))
thresholds <- c(5,10,15)

# Repeat the same full-sample network logic as script 06 for a subset.
fit_network <- function(keep) {
 mm <- m[keep,,drop=FALSE]; vv <- v[keep,,drop=FALSE]; pp <- pred[keep]
 corm <- pp=='cormorant'; seal <- pp=='harbor_seal'
 if(!any(corm)||!any(seal)) stop('Missing predator group after jackknife')
 fc <- colSums(mm[corm,,drop=FALSE]>0)/sum(corm)*100
 fs <- colSums(mm[seal,,drop=FALSE]>0)/sum(seal)*100
 fd <- data.frame(taxon_short=fish_short,foo_corm=as.numeric(fc),foo_seal=as.numeric(fs),mean_foo=(as.numeric(fc)+as.numeric(fs))/2)
 top <- fd %>% group_by(taxon_short) %>% summarise(foo_corm=max(foo_corm),foo_seal=max(foo_seal),mean_foo=max(mean_foo),.groups='drop') %>% arrange(desc(mean_foo)) %>% slice_head(n=8)
 if(nrow(top)!=8L) stop('Insufficient top-eight fish')
 pf <- dplyr::bind_rows(data.frame(from='Great cormorant',to=top$taxon_short,weight=top$foo_corm,predator='Great cormorant'),data.frame(from='Harbour seal',to=top$taxon_short,weight=top$foo_seal,predator='Harbour seal')) %>% filter(weight>5)
 raw <- vector('list',0)
 for(fish in top$taxon_short) {
  present <- rowSums(mm[,fish_short==fish,drop=FALSE])>0
  denom <- sum(present)
  if(denom==0L) next
  for(j in seq_len(ncol(vv))) {
   cnt <- sum(present & vv[,j]>0)
   if(cnt>0L) raw[[length(raw)+1L]] <- data.frame(fish=fish,inv=inv_short[j],source_invertebrate=colnames(vv)[j],cooccur=cnt,pct=100*cnt/denom)
  }
 }
 if(length(raw)) {
  raw <- dplyr::bind_rows(raw)
  fi <- raw %>% group_by(fish,inv) %>% summarise(cooccur=max(cooccur),pct=max(pct),.groups='drop')
 } else fi <- data.frame(fish=character(),inv=character(),cooccur=integer(),pct=numeric())
 list(top=top,pf=pf,fi=fi,n_corm=sum(corm),n_seal=sum(seal))
}

key_pf <- function(x) paste(x$from,x$to,sep=' || ')
key_fi <- function(x) paste(x$fish,x$inv,sep=' || ')
jaccard <- function(a,b) {
 un <- union(a,b)
 if(!length(un)) return(1)
 length(intersect(a,b))/length(un)
}

base <- fit_network(rep(TRUE,length(shared)))
ref_top <- read.csv(ref_fish,stringsAsFactors=FALSE)
ref_pf_df <- read.csv(ref_pf,stringsAsFactors=FALSE)
ref_fi_df <- read.csv(ref_fi,stringsAsFactors=FALSE)
# Exact 10%-baseline agreement with already approved module 06.
base10 <- base$fi[base$fi$pct>=10,,drop=FALSE]
if(!setequal(base$top$taxon_short,ref_top$taxon_short) ||
   !setequal(key_pf(base$pf),key_pf(ref_pf_df)) ||
   !setequal(key_fi(base10),key_fi(ref_fi_df))) {
 stop('10% baseline does not match approved script 06. Stop before jackknife.')
}
check_join <- function(a,b,keyfun,cols,label) {
 ka<-keyfun(a); kb<-keyfun(b)
 ia<-match(kb,ka)
 for(col in cols) {
  if(!all(col%in%names(a),col%in%names(b))) stop('Missing comparison field ',col)
  delta <- abs(as.numeric(a[[col]][ia])-as.numeric(b[[col]]))
  if(anyNA(delta)||any(delta>1e-7)) stop('Baseline weights differ: ',label,' ',col)
 }
}
check_join(base$pf,ref_pf_df,key_pf,'weight','predator-fish')
check_join(base10,ref_fi_df,key_fi,c('pct','cooccur'),'fish-invertebrate')
cat('Baseline 10% network MATCHES approved script 06 (edges and weights).\n')

# Sensitivity to 5,10,15% fish-invertebrate co-occurrence threshold.
all_baselines <- lapply(thresholds,function(t) {
 fi <- base$fi[base$fi$pct>=t,,drop=FALSE]
 data.frame(threshold_pct=t,n_shared=length(shared),n_predator_fish=nrow(base$pf),n_fish_invertebrate=nrow(fi),n_invertebrate_nodes=length(unique(fi$inv)),n_top8=nrow(base$top),n_fish_invertebrate_links_from_10=length(intersect(key_fi(fi),key_fi(base10))),jaccard_with_10=jaccard(key_fi(fi),key_fi(base10)))
})
sens <- do.call(rbind,all_baselines)
write.csv(sens,file.path(VAL,'07B_THRESHOLD_SUMMARY.csv'),row.names=FALSE)
for(t in thresholds) {
 fi <- base$fi[base$fi$pct>=t,,drop=FALSE]
 write.csv(fi,file.path(VAL,sprintf('06B_FISH_INVERTEBRATE_EDGES_%02dpct.csv',t)),row.names=FALSE)
}

# Full leave-one-out jackknife: rerank top-eight fish on each subset,
# recalculate predator FOO and fish-invertebrate co-occurrence denominators.
jack_summary <- vector('list',length(shared)*length(thresholds))
jack_pf <- vector('list',length(shared))
jack_fish <- vector('list',length(shared))
jack_fi <- vector('list',length(shared)*length(thresholds))
k <- 0L
for(i in seq_along(shared)) {
 keep <- seq_along(shared)!=i
 one <- fit_network(keep)
 jack_pf[[i]] <- data.frame(left_out=shared[i],key=key_pf(one$pf),weight=one$pf$weight)
 jack_fish[[i]] <- data.frame(left_out=shared[i],fish=one$top$taxon_short)
 for(t in thresholds) {
  k <- k+1L
  fi <- one$fi[one$fi$pct>=t,,drop=FALSE]
  base_fi <- base$fi[base$fi$pct>=t,,drop=FALSE]
  jack_fi[[k]] <- data.frame(left_out=shared[i],threshold_pct=t,key=key_fi(fi),pct=fi$pct)
  jack_summary[[k]] <- data.frame(left_out=shared[i],predator=pred[i],threshold_pct=t,n_samples=sum(keep),n_predator_fish=nrow(one$pf),n_fish_invertebrate=nrow(fi),n_invertebrate_nodes=length(unique(fi$inv)),top8_overlap=length(intersect(one$top$taxon_short,base$top$taxon_short)),pf_jaccard=jaccard(key_pf(one$pf),key_pf(base$pf)),fi_jaccard=jaccard(key_fi(fi),key_fi(base_fi)),links_lost=length(setdiff(key_fi(base_fi),key_fi(fi))),links_gained=length(setdiff(key_fi(fi),key_fi(base_fi))))
 }
 if(i%%25L==0L) cat('Jackknife processed ',i,'/',length(shared),' samples\n',sep='')
}
jack_summary <- do.call(rbind,jack_summary)
jack_pf <- do.call(rbind,jack_pf)
jack_fish <- do.call(rbind,jack_fish)
jack_fi <- do.call(rbind,jack_fi)
write.csv(jack_summary,file.path(VAL,'07B_JACKKNIFE_LEAVE_ONE_OUT.csv'),row.names=FALSE)
write.csv(jack_pf,file.path(VAL,'07B_JACKKNIFE_PREDATOR_FISH_EDGES.csv'),row.names=FALSE)
write.csv(jack_fish,file.path(VAL,'07B_JACKKNIFE_TOP8_FISH.csv'),row.names=FALSE)
write.csv(jack_fi,file.path(VAL,'07B_JACKKNIFE_FISH_INVERTEBRATE_EDGES.csv'),row.names=FALSE)

# Inclusion frequencies across ALL leave-one-out replicates.
# A selected edge has frequency 1 only if present in every replicate.
edge_freq <- function(df,baseline_keys,edge_class,t) {
 observed <- table(df$key)
 all_keys <- union(as.character(names(observed)),baseline_keys)
 count <- setNames(integer(length(all_keys)),all_keys)
 count[names(observed)] <- as.integer(observed)
 data.frame(edge_class=edge_class,threshold_pct=t,edge=all_keys,baseline_selected=all_keys%in%baseline_keys,n_LOO_selected=unname(count),selection_frequency=unname(count)/length(shared))
}
fr <- list(edge_freq(jack_pf,key_pf(base$pf),'predator_fish',NA_real_))
for(t in thresholds) {
 jj <- jack_fi[jack_fi$threshold_pct==t,,drop=FALSE]
 bb <- base$fi[base$fi$pct>=t,,drop=FALSE]
 fr[[length(fr)+1L]] <- edge_freq(jj,key_fi(bb),'fish_invertebrate',t)
}
frequencies <- do.call(rbind,fr)
write.csv(frequencies,file.path(VAL,'07B_JACKKNIFE_EDGE_SELECTION_FREQUENCIES.csv'),row.names=FALSE)
fish_inclusion <- as.data.frame(table(jack_fish$fish),stringsAsFactors=FALSE)
names(fish_inclusion) <- c('fish','n_LOO_top8')
fish_inclusion$selection_frequency <- fish_inclusion$n_LOO_top8/length(shared)
fish_inclusion$baseline_top8 <- fish_inclusion$fish%in%base$top$taxon_short
write.csv(fish_inclusion,file.path(VAL,'07B_JACKKNIFE_TOP8_STABILITY.csv'),row.names=FALSE)

# Aggregate leave-one-out sensitivity by threshold.
jack_agg <- jack_summary %>% group_by(threshold_pct) %>% summarise(n_replicates=n(),mean_fi_edges=mean(n_fish_invertebrate),min_fi_edges=min(n_fish_invertebrate),max_fi_edges=max(n_fish_invertebrate),mean_fi_jaccard=mean(fi_jaccard),min_fi_jaccard=min(fi_jaccard),mean_pf_jaccard=mean(pf_jaccard),min_pf_jaccard=min(pf_jaccard),mean_top8_overlap=mean(top8_overlap),min_top8_overlap=min(top8_overlap),.groups='drop')
write.csv(jack_agg,file.path(VAL,'07B_JACKKNIFE_SUMMARY_BY_THRESHOLD.csv'),row.names=FALSE)

# QA figures in validation only; approved trophic_network.png remains unchanged.
p1 <- ggplot(sens,aes(threshold_pct,n_fish_invertebrate))+geom_line(linewidth=.7)+geom_point(size=2.5)+scale_x_continuous(breaks=thresholds)+labs(x='Fish-invertebrate co-occurrence threshold (%)',y='Selected fish-invertebrate edges',title='Network sensitivity to threshold')+theme_classic()
ggsave(file.path(VAL,'07B_THRESHOLD_SENSITIVITY.png'),p1,width=7,height=4.5,dpi=300)
p2 <- ggplot(jack_summary,aes(factor(threshold_pct),fi_jaccard))+geom_boxplot(outlier.alpha=.4)+labs(x='Co-occurrence threshold (%)',y='Jaccard vs full-data network',title='Leave-one-out stability of fish-invertebrate edges')+coord_cartesian(ylim=c(0,1))+theme_classic()
ggsave(file.path(VAL,'07B_JACKKNIFE_STABILITY.png'),p2,width=7,height=4.5,dpi=300)

write.csv(data.frame(file=required_files,md5=unname(tools::md5sum(required_files))),file.path(VAL,'07B_INPUT_PROVENANCE.csv'),row.names=FALSE)
writeLines(c(
 'APPROVED FIGURE 06 UNCHANGED. This module provides QA only.',
 'Underlying ASV filter is unchanged: global 0.1% identified phyloseq objects from 02.',
 'Predator-fish: strict FOO >5%; does NOT vary with sensitivity threshold.',
 'Fish-invertebrate: inclusive co-occurrence >=5%, >=10%, >=15%.',
 'Original definition: count shared samples with fish present and a given COI taxon present; divide by number of shared samples where fish is present; group abbreviated invertebrate labels with max().',
 'Top eight fish selected using max of original per-LCA FOO and max of per-LCA mean FOO.',
 'Leave-one-out jackknife over MATCHED samples; top eight taxa and every edge are re-estimated in each replicate.',
 'Jaccard=intersection/union of edge sets; empty-empty defined as 1.',
 'Selection frequency=number of leave-one-out networks containing edge / N shared samples.',
 'These are leave-one-out influence diagnostics, not independent replication, Bayesian probabilities, or standard CI.',
 'A co-occurrence edge DOES NOT establish that a fish consumed an invertebrate.',
 'Baseline at 10% must match approved script 06 edge IDs AND numerical weights.',
 'All output files produced under ../../docs/validation/07B_Trophic_Network_Sensitivity_Jackknife.'
 ),file.path(VAL,'07B_METHOD_NOTES.txt'))
capture.output(sessionInfo(),file=file.path(VAL,'07B_SESSION_INFO.txt'))
cat('\n06B COMPLETE: ',VAL,'\n',sep='')
print(sens)
print(jack_agg)
