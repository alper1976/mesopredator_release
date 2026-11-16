# PORTABLE STANDALONE RECONSTRUCTION. Original calculations retained.
# ============================================================================
# 09_Functional_Group_Composition.R
# Mesopredator_check_2 | standalone manuscript figure reconstruction
# TARGET: FG_composition_PHOCA_VITULINA_COI.png (Tarjei thesis Figure 3.15)
#
# INPUT: validated script 02 COI phyloseq object (no re-fitting or upstream edits)
# IMPORTANT SOURCE LIMITATION:
# The available original Rmd does not contain the exact Figure 3.15 export code
# or an exhaustive COI functional-group lookup. The ten explicitly documented
# assignments are retained; other fish taxa have a transparent, provisional
# genus-based ecological assignment and non-fish metazoans use their Class.
# Thus visual/numerical identity with the original MUST be checked before approval.
# ============================================================================

required <- c('phyloseq','ggplot2','dplyr','tidyr')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install required packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({
  library(phyloseq); library(ggplot2); library(dplyr); library(tidyr)
})

# ---- 1. Folder layout ---------------------------------------------------------
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
INPUT <- file.path(BASE,'processed_data','03_Phyloseq_Preprocessing',
                   'ps_coi_global_0_1pct_identified.rds')
OUT <- file.path(BASE,'processed_data','09_Functional_Group_Composition')
FIG <- file.path(OUT,'figures')
TAB <- file.path(OUT,'tables')
VAL <- file.path(BASE,'validation','09_Functional_Group_Composition')
for(p in c(FIG,TAB,VAL)) dir.create(p,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(INPUT)) stop('Script 02 COI phyloseq input missing: ',INPUT)

# ---- 2. Input data, with sample / taxonomy ID alignment ----------------------
ps <- readRDS(INPUT)
stopifnot(inherits(ps,'phyloseq'))
Y <- as(phyloseq::otu_table(ps),'matrix')
if(phyloseq::taxa_are_rows(ps)) Y <- t(Y)
storage.mode(Y) <- 'double'
md <- as.data.frame(phyloseq::sample_data(ps),check.names=FALSE)
tx <- as.data.frame(phyloseq::tax_table(ps),check.names=FALSE)
if(!setequal(rownames(Y),rownames(md)) ||
   !setequal(colnames(Y),rownames(tx))) stop('Phyloseq metadata/ASV IDs not aligned')
md <- md[rownames(Y),,drop=FALSE]
tx <- tx[colnames(Y),,drop=FALSE]
stopifnot(all(c('species','season','sample_year') %in% names(md)),
          all(c('LCA_taxon','Class') %in% names(tx)),
          !anyNA(Y),all(is.finite(Y)),all(Y>=0))

# Harbour seal only; use the COI sample cohort from script 02, not PCoA subset.
sp <- trimws(as.character(md[['species']]))
sp[sp=='top_cormorant'] <- 'cormorant'
keep <- sp=='harbor_seal'
if(!any(keep)) stop('No harbour seal samples found')
Y <- Y[keep,,drop=FALSE]
md <- md[keep,,drop=FALSE]

# Keep samples with positive identified-COI signal. No fictitious 56-sample
# subset is imposed (original figure has 56 labelled valid samples).
pos <- rowSums(Y)>0
zero_ids <- rownames(Y)[!pos]
Y <- Y[pos,,drop=FALSE]
md <- md[pos,,drop=FALSE]
if(!nrow(Y)) stop('No positive-signal seal samples')

# ---- 3. Reproduce the full-LCA sample-level RRA method ----------------------
lca <- gsub('_',' ',trimws(as.character(tx$LCA_taxon)),fixed=TRUE)
if(anyNA(lca)||any(!nzchar(lca))) stop('Missing LCA names')
colnames(Y) <- lca
agg <- t(rowsum(t(Y),group=lca,reorder=TRUE))
stopifnot(identical(rownames(agg),rownames(md)), all(rowSums(agg)>0))
rra <- sweep(agg,1,rowSums(agg),'/')*100
stopifnot(max(abs(rowSums(rra)-100))<1e-8)

# One class per LCA name; flag mixed classes for manual review.
class_by_lca <- split(as.character(tx$Class),lca)
class_txt <- vapply(class_by_lca,function(x)
  paste(sort(unique(x[!is.na(x)])),collapse=' | '),character(1))
class_txt <- class_txt[colnames(rra)]

# ---- 4. Functional assignments: provenance always retained -----------------
# Exact ten-name lookup from original COI diet code (not a complete source map).
original_exact <- c(
  'Clupea harengus'='Pelagic',
  'Sprattus sprattus'='Pelagic',
  'Micromesistius poutassou fish environmental sample'='Pelagic',
  'Pollachius virens Pollachius environmental sample Actinopterygii environmental sample'='Demersal',
  'Merlangius merlangus'='Demersal',
  'Ctenolabrus rupestris'='Demersal',
  'Pleuronectes platessa Pleuronectes sp. SH-2018'='Benthic',
  'Myoxocephalus scorpius Myoxocephalus scorpioides Myoxocephalus aenaeus Myoxocephalus verrucosus'='Benthic',
  'Pomatoschistus minutus Pomatoschistus microps Neogobius melanostomus Gobius niger'='Benthic',
  'Pontophilus norvegicus'='Invertebrates'
)

# Provisional cross-marker ecological assignments of the remaining fish.
# Based on the same functional definitions used in Tarjei's MiFish lookup.
# These are hypotheses, NOT verified as the original Figure 3.15 classification.
fish_genera <- c(
  Ammodytes='Pelagic', Argentina='Pelagic', Belone='Pelagic',
  Clupea='Pelagic', Gasterosteus='Pelagic', Hyperoplus='Pelagic',
  Micromesistius='Pelagic', Scomber='Pelagic', Spinachia='Pelagic',
  Sprattus='Pelagic',
  Aphia='Demersal', Callionymus='Demersal', Centrolabrus='Demersal',
  Ctenolabrus='Demersal', Enchelyopus='Demersal', Gadus='Demersal',
  Gaidropsarus='Demersal', Gobiusculus='Demersal', Labrus='Demersal',
  Melanogrammus='Demersal', Merlangius='Demersal', Merluccius='Demersal',
  Molva='Demersal', Pollachius='Demersal', Raniceps='Demersal',
  Symphodus='Demersal', Trisopterus='Demersal',
  Anguilla='Benthic', Arnoglossus='Benthic', Glyptocephalus='Benthic',
  Gobius='Benthic', Limanda='Benthic', Lumpenus='Benthic',
  Lycodes='Benthic', Microstomus='Benthic', Myoxocephalus='Benthic',
  Platichthys='Benthic', Pleuronectes='Benthic', Pomatoschistus='Benthic',
  Scophthalmus='Benthic', Solea='Benthic', Taurulus='Benthic',
  Zoarces='Benthic'
)

assign_group <- function(label,cls) {
  if(label %in% names(original_exact))
    return(c(group=unname(original_exact[[label]]),source='ORIGINAL_EXACT'))
  if(!identical(cls,'Actinopteri') && !grepl('\\|',cls) && nzchar(cls))
    return(c(group='Invertebrates',source='CLASS_NON_ACTINOPTERI'))
  if(identical(cls,'Actinopteri')) {
    genus <- strsplit(label,' ',fixed=TRUE)[[1]][1]
    if(genus %in% names(fish_genera))
      return(c(group=unname(fish_genera[[genus]]),source='PROVISIONAL_FISH_GENUS'))
  }
  c(group='Unassigned',source='SOURCE_GAP')
}
assignments <- t(vapply(seq_along(class_txt),function(i)
  assign_group(names(class_txt)[i],unname(class_txt[[i]])),character(2)))
lookup <- data.frame(
  LCA_taxon=names(class_txt),
  Class=unname(class_txt),
  functional_group=unname(assignments[,'group']),
  assignment_source=unname(assignments[,'source']),
  stringsAsFactors=FALSE
)
write.csv(lookup,file.path(VAL,'09_TAXON_FUNCTIONAL_MAPPING_AUDIT.csv'),row.names=FALSE)

# ---- 5. Season and time-period metadata -------------------------------------
seasons <- c('Spring','Summer','Autumn','Winter')
periods <- c('Early (2012-2016)','Mid (2019-2021)','Recent (2022-2025)')
year <- suppressWarnings(as.integer(as.character(md[['sample_year']])))
season <- tools::toTitleCase(tolower(trimws(as.character(md[['season']]))))
period <- ifelse(year %in% 2012:2016,periods[1],
          ifelse(year %in% 2019:2021,periods[2],
          ifelse(year %in% 2022:2025,periods[3],NA_character_)))
metadata <- data.frame(sample_id=rownames(md),season=season,time_period=period,
                       sample_year=year,stringsAsFactors=FALSE)
write.csv(metadata,file.path(VAL,'09_ALL_SEAL_SAMPLE_METADATA.csv'),row.names=FALSE)
valid <- !is.na(period) & season %in% seasons
if(!all(valid)) warning(sum(!valid),' seal samples excluded because period/season is unknown; see audit')
metadata <- metadata[valid,,drop=FALSE]
rra <- rra[metadata$sample_id,,drop=FALSE]
stopifnot(identical(rownames(rra),metadata$sample_id))

# ---- 6. Sum group RRA within each sample; THEN average among samples --------
fg_names <- c('Pelagic','Demersal','Benthic','Invertebrates','Unassigned')
M <- sapply(fg_names,function(g) {
  idx <- lookup$functional_group==g
  if(!any(idx)) rep(0,nrow(rra)) else rowSums(rra[,idx,drop=FALSE])
})
rownames(M) <- rownames(rra)
if(any(abs(rowSums(M)-100)>1e-6)) stop('Functional-group RRA does not sum to 100% per sample')
per_sample <- data.frame(sample_id=rownames(M),M,check.names=FALSE)
per_sample <- merge(per_sample,metadata,by='sample_id',sort=FALSE)
write.csv(per_sample,file.path(VAL,'09_SAMPLE_FUNCTIONAL_RRA.csv'),row.names=FALSE)

long <- tidyr::pivot_longer(per_sample,cols=dplyr::all_of(fg_names),
                            names_to='functional_group',values_to='RRA_percent')
summary <- long |>
  dplyr::group_by(season,time_period,functional_group) |>
  dplyr::summarise(mean_RRA=mean(RRA_percent),n_samples=dplyr::n(),.groups='drop')

# Make all 12 season-by-period combinations explicit; 'No data' means n<3.
combos <- tidyr::expand_grid(season=seasons,time_period=periods)
counts <- metadata |>
  dplyr::count(season,time_period,name='n_samples')
combos <- dplyr::left_join(combos,counts,by=c('season','time_period'))
combos$n_samples[is.na(combos$n_samples)] <- 0L
summary <- tidyr::expand_grid(season=seasons,time_period=periods,
                              functional_group=fg_names) |>
  dplyr::left_join(summary[,c('season','time_period','functional_group','mean_RRA')],
                   by=c('season','time_period','functional_group')) |>
  dplyr::left_join(combos,by=c('season','time_period')) |>
  dplyr::mutate(is_plottable=n_samples>=3,
                mean_RRA=ifelse(is_plottable,coalesce(mean_RRA,0),NA_real_))

sums <- summary |>
  dplyr::filter(is_plottable) |>
  dplyr::group_by(season,time_period) |>
  dplyr::summarise(total=sum(mean_RRA),.groups='drop')
if(any(abs(sums$total-100)>1e-5)) stop('Non-100% mean composition in one or more cells')

summary$season <- factor(summary$season,levels=seasons)
summary$time_period <- factor(summary$time_period,levels=periods)
summary$functional_group <- factor(summary$functional_group,levels=fg_names)
combos$season <- factor(combos$season,levels=seasons)
combos$time_period <- factor(combos$time_period,levels=periods)
write.csv(summary,file.path(TAB,'FG_composition_PHOCA_VITULINA_COI_data.csv'),row.names=FALSE)
write.csv(combos,file.path(VAL,'09_SEASON_PERIOD_SAMPLE_COUNTS.csv'),row.names=FALSE)

# Compare against the labelled counts in Tarjei Figure 3.15: 56 samples.
# This is an audit, not a hard-coded subset or a pass condition.
reference <- data.frame(
  season=c('Spring','Spring','Summer','Autumn','Autumn','Winter','Winter'),
  time_period=c(periods[1],periods[3],periods[1],periods[1],periods[2],periods[1],periods[3]),
  original_n=c(8L,5L,4L,11L,11L,8L,9L),stringsAsFactors=FALSE)
check <- dplyr::left_join(reference,counts,by=c('season','time_period'))
check$n_samples[is.na(check$n_samples)] <- 0L
check$matches_original <- check$original_n==check$n_samples
write.csv(check,file.path(VAL,'09_ORIGINAL_FIGURE_N_COMPARISON.csv'),row.names=FALSE)
if(!all(check$matches_original)) warning(
  'Sample counts do not exactly match Tarjei Figure 3.15. ',
  'See ../../docs/validation/09_ORIGINAL_FIGURE_N_COMPARISON.csv. Do not silently remove samples.')

# ---- 7. Thesis-style four-season faceted stacked bars ------------------------
colours <- c(Pelagic='#1f77c9',Demersal='#8b4513',Benthic='#8a008e',
             Invertebrates='#2e8b57',Unassigned='#d3d3d3')
bar_data <- summary[summary$is_plottable,,drop=FALSE]
no_data <- combos[combos$n_samples<3,,drop=FALSE]
labels <- combos[combos$n_samples>=3,,drop=FALSE]

plot <- ggplot() +
  geom_col(data=bar_data,
    aes(x=time_period,y=mean_RRA,fill=functional_group),width=0.9,
    colour='white',linewidth=0.2) +
  geom_col(data=no_data,aes(x=time_period,y=100),
    fill='#eeeeee',colour='#cccccc',width=0.9) +
  geom_text(data=no_data,aes(x=time_period,y=50,label='No data'),
    angle=90,colour='#999999',size=3.5) +
  geom_text(data=labels,aes(x=time_period,y=105,label=paste0('n=',n_samples)),size=3.6) +
  facet_grid(.~season,drop=FALSE) +
  scale_x_discrete(labels=c(
    'Early (2012-2016)'='Early\n(2012–2016)',
    'Mid (2019-2021)'='Mid\n(2019–2021)',
    'Recent (2022-2025)'='Recent\n(2022–2025)'
  )) +
  scale_fill_manual(values=colours,breaks=fg_names,drop=FALSE,name='Functional Group') +
  scale_y_continuous(limits=c(0,115),breaks=c(0,25,50,75,100),
                     expand=expansion(mult=c(0,0))) +
  labs(title=expression('Functional group composition across season and time period — '*italic('Phoca vitulina')*' (COI Leray-XT)'),
       x='Time Period',y='Mean relative abundance (%)') +
  theme_classic(base_size=13) +
  theme(plot.title=element_text(hjust=0.5,size=14,face='bold'),
        strip.background=element_blank(),strip.text=element_text(face='bold',size=13),
        axis.text.x=element_text(size=10,lineheight=0.95,hjust=0.5,margin=margin(t=5)),
        axis.text.y=element_text(size=11),
        legend.position='right',
        panel.spacing.x=grid::unit(0.5,'lines'),
        panel.grid.major.y=element_line(colour='#eeeeee',linewidth=0.3))

fig_path <- file.path(FIG,'FG_composition_PHOCA_VITULINA_COI.png')
ggsave(fig_path,plot,width=19,height=7.5,dpi=300,bg='white')

# ---- 8. Reconstruction differences / transparent QA -------------------------
# Do not remove any samples or modify any assignments to match reference bars.
# This is a presentation/documentation update; numerical analysis is unchanged.
ref_comparison <- check
ref_comparison$extra_in_reconstruction <- (ref_comparison$n_samples -
                                            ref_comparison$original_n)
write.csv(ref_comparison,
          file.path(VAL,'09_SAMPLE_COUNT_DEVIATIONS.csv'),row.names=FALSE)

# Save the identities of *all* samples in discrepant cells as candidates for
# reconciliation. Without Tarjei's original per-sample inclusion table, it is
# not possible to identify which specific seven were omitted.
discrepant_cells <- ref_comparison[ref_comparison$extra_in_reconstruction!=0,
                                  c('season','time_period'),drop=FALSE]
subject_audit <- dplyr::inner_join(metadata,discrepant_cells,
                                   by=c('season','time_period'))
write.csv(subject_audit,
          file.path(VAL,'09_DISCREPANT_GROUP_SAMPLE_IDS.csv'),row.names=FALSE)

# Quantify the abundance assigned provisionally to demersal taxa.
dem_ref <- lookup[lookup$functional_group=='Demersal',,drop=FALSE]
dem_orig <- dem_ref$LCA_taxon[dem_ref$assignment_source=='ORIGINAL_EXACT']
dem_prov <- dem_ref$LCA_taxon[dem_ref$assignment_source=='PROVISIONAL_FISH_GENUS']
if(!all(c(dem_orig,dem_prov)%in%colnames(rra))) stop('Demersal audit taxa not aligned')
sum_taxa <- function(cols) {
  if(length(cols)==0L) return(rep(0,nrow(rra)))
  rowSums(rra[,cols,drop=FALSE])
}
dem_audit <- data.frame(sample_id=rownames(rra),
  ORIGINAL_EXACT=sum_taxa(dem_orig),
  PROVISIONAL_FISH_GENUS=sum_taxa(dem_prov),stringsAsFactors=FALSE)
dem_audit <- dplyr::left_join(dem_audit,metadata,by='sample_id')
dem_period <- dem_audit |>
  dplyr::group_by(season,time_period) |>
  dplyr::summarise(n_samples=dplyr::n(),
    original_exact_mean_RRA=mean(ORIGINAL_EXACT),
    provisional_mean_RRA=mean(PROVISIONAL_FISH_GENUS),
    total_demersal_mean_RRA=mean(ORIGINAL_EXACT+PROVISIONAL_FISH_GENUS),
    .groups='drop')
write.csv(dem_period,file.path(VAL,'09_DEMERSAL_SOURCE_BY_SEASON_PERIOD.csv'),
          row.names=FALSE)

n_original <- sum(ref_comparison$original_n)
n_our <- sum(ref_comparison$n_samples)
n_provisional <- sum(lookup$assignment_source=='PROVISIONAL_FISH_GENUS')
n_origexact <- sum(lookup$assignment_source=='ORIGINAL_EXACT')
notes <- c(
 '08 — FUNCTIONAL GROUP COMPOSITION: METHODS AND RECONSTRUCTION DIFFERENCES',
 '',
 'Status: reconstructed Figure 3.15 candidate; sample cohort and taxon grouping',
 'are not demonstrated to be numerically identical to the original.',
 '',
 'Calculation (unchanged): script 02 COI globally filtered identified ASVs;',
 'positive-signal harbour-seal samples, sum ASVs by full LCA label, scale each',
 'sample to 100%, sum functional groups within samples, then average within',
 'season × time-period cells. Grey No data bars for n < 3.',
 'Only the x-axis labels and export width changed in this revision.',
 '',
 paste0('Original Figure 3.15 labeled sample count = ',n_original,
        '; reconstructed cells = ',n_our,
        '; difference = +',n_our-n_original,'.'),
 'Breakdown of extra samples (source figure counts vs script 02 cohort):',
 '  Summer / Early: original 4; reconstruction 5 (+1)',
 '  Autumn / Early: original 11; reconstruction 12 (+1)',
 '  Autumn / Mid: original 11; reconstruction 13 (+2)',
 '  Winter / Early: original 8; reconstruction 11 (+3)',
 '  Spring / Early, Spring / Recent, Winter / Recent: counts match.',
 'The seven excess observations cannot be assigned to specific sample IDs',
 'without the original per-sample inclusion/exclusion list.',
 'No samples are removed or chosen solely to force a match to the figure.',
 'See 09_SAMPLE_COUNT_DEVIATIONS.csv and 08_DISCREPANT_GROUP_SAMPLE_IDS.csv.',
 '',
 'Functional-group mapping differs from a confirmed source-derived full map.',
 paste0('Direct original COI diet lookup assignments: ',n_origexact,'.'),
 paste0('Provisional fish-genus assignments: ',n_provisional,'.'),
 'Non-fish metazoans are classified using taxonomic Class.',
 'These provisional assignments are explicitly flagged, not claimed to be',
 'Tarjei’s exact Figure 3.15 functional-group dictionary.',
 'Among identified demersal taxa, Molva molva, Labrus bergylta and',
 'Trisopterus minutus are important provisionally assigned contributors.',
 'The extended mapping can affect the fraction labelled Demersal, especially',
 'Autumn Early and Autumn Mid; sample-cohort differences may also contribute.',
 'The difference should not be attributed exclusively to taxonomy.',
 'See 09_TAXON_FUNCTIONAL_MAPPING_AUDIT.csv and',
 '09_DEMERSAL_SOURCE_BY_SEASON_PERIOD.csv for reproducible numerical details.',
 '',
 'No COI model re-fitting, PERMANOVA, manual abundance adjustment, sample',
 'deletion, or changes to functional-group assignments were made.',
 'The source original figure export code/full lookup has not been verified.',
 'Therefore describe this as reconstructed with documented deviations,',
 'not an exact numerical reproduction.'
)
writeLines(notes,file.path(VAL,'09_METHOD_AND_DEVIATION_NOTE.txt'))

# ---- 8. Provenance and explicit reconstruction limitations ------------------
write.csv(data.frame(source=INPUT,md5=unname(tools::md5sum(INPUT))),
          file.path(VAL,'09_INPUT_PROVENANCE.csv'),row.names=FALSE)
writeLines(c(
  '08 Functional Group Composition: Phoca vitulina, COI marker.',
  'Source: Tarjei thesis Figure 3.15 and original COI functional group lookup.',
  'Input: approved script 02 COI global 0.1%-filter, identified ASVs.',
  'Full LCA abundance: aggregate ASVs, normalize per positive-signal sample to 100%,',
  'sum by functional group within sample, average groups within season and time-period.',
  'Periods: 2012-2016 Early, 2019-2021 Mid, 2022-2025 Recent.',
  'Missing/low-n panels: grey No data when n<3.',
  'Ten functional assignments are explicitly supported by original COI diet code.',
  'Other metazoans are grouped using Class; additional fish use a provisional',
  'genus lookup based on prior MiFish ecology assignments, not verified Figure 3.15 code.',
  'All uncertain mapping decisions are recorded in validation mapping audit.',
  'Original figure shows 56 labelled positive/valid samples in populated cells.',
  'Script DOES NOT force 56 samples: matching requires sample ID reconciliation.',
  'No PERMANOVA is rerun and no statistical claim is inferred from this plot.',
  'Figure is reconstructed with documented cohort and taxonomy deviations.',
  '09_METHOD_AND_DEVIATION_NOTE.txt records detailed causes, limitations and audits.'
),file.path(VAL,'09_METHOD_NOTES.txt'))
capture.output(sessionInfo(),file=file.path(VAL,'09_SESSION_INFO.txt'))
cat('\n08 COMPLETE (UNVALIDATED RECONSTRUCTION CANDIDATE)\nFigure: ',
    fig_path,'\nReference sample counts matched: ',all(check$matches_original),
    '\nMapping source tally:\n',sep='')
print(table(lookup$assignment_source))
print(check,row.names=FALSE)
