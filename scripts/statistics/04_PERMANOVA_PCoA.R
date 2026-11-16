# PORTABLE STANDALONE RECONSTRUCTION. Original calculations retained.
# 04_PERMANOVA_PCoA.R
# Final reconstruction candidate — Tarjei Røstum MiFish/COI PERMANOVA,
# seasonal betadisper and predator-specific PCoA plots.
# ORIGINAL SOURCES:
#   0.1.PERMANOVA.PCOA.Rmd (final ps_01 run_analysis block, n>=2)
#   COI.PERMANOVA.Rmd
# INPUT: ONLY derived phyloseq objects saved by 03_Phyloseq_Preprocessing.R.
# Manuscript outputs reproduce original filenames and ggplot geometry.
# Validation tables/diagnostics go to ../../docs/validation/, not outputs/.

required <- c('phyloseq','vegan','ggplot2')
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Install required packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({library(phyloseq);library(vegan);library(ggplot2)})
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
OUT <- file.path(BASE,'processed_data','04_PERMANOVA_PCoA')
FIG <- file.path(OUT,'figures')
VAL <- file.path(BASE,'validation','04_PERMANOVA_PCoA')
for(d in c(FIG,VAL)) dir.create(d,recursive=TRUE,showWarnings=FALSE)
inputs <- c(MiFish='ps_mifish_global_0_1pct.rds',COI='ps_coi_global_0_1pct_identified.rds')
paths <- file.path(INPUT,unname(inputs)); names(paths) <- names(inputs)
if(!all(file.exists(paths))) stop('Missing script-02 outputs: ',paste(paths[!file.exists(paths)],collapse=' | '))
ps_mifish <- readRDS(paths[['MiFish']])
ps_coi <- readRDS(paths[['COI']])
stopifnot(inherits(ps_mifish,'phyloseq'),inherits(ps_coi,'phyloseq'),
          nsamples(ps_mifish)==161L,ntaxa(ps_mifish)==160L,
          nsamples(ps_coi)==145L,ntaxa(ps_coi)==172L)
# Original MiFish ps_01 analysis excludes taxa with LCA_taxon == 'no identification'.
# This object retains the 12 zero-signal samples until analysis-specific filtering.
tax <- as.data.frame(tax_table(ps_mifish),check.names=FALSE)
if(!'LCA_taxon' %in% names(tax)) stop('MiFish lacks LCA_taxon')
no_id <- !is.na(tax$LCA_taxon) & tolower(trimws(as.character(tax$LCA_taxon)))=='no identification'
if(sum(no_id)!=91L) stop('Expected 91 unidentified MiFish ASVs')
ps_mifish <- prune_taxa(rownames(tax)[!no_id],ps_mifish)
stopifnot(ntaxa(ps_mifish)==69L,sum(sample_sums(ps_mifish)==0)==12L)

season_colors <- c(Autumn='#CD853F',Spring='#90EE90',Summer='#2E8B57',Winter='#ADD8E6')
period_colors <- c('Early (2012-2016)'='#08306B','Mid (2019-2021)'='#33A02C',
                   'Recent (2022-2025)'='#E31A1C')
season_shapes <- c(Autumn=17,Spring=16,Summer=15,Winter=18)
# Original GG2 figure: individual YEAR colours (not period colours for points).
# 2012 is included as a defensively defined level; original 11-colour vector
# starts at 2013. The exact year values in each panel are checked at runtime.
year_colors <- c(
  '2012'='#08306B','2013'='#08306B','2014'='#2171B5',
  '2015'='#E31A1C','2016'='#FF7F00','2019'='#33A02C',
  '2020'='#6A3D9A','2021'='#FFD700','2022'='#B15928',
  '2023'='#FB9A99','2024'='#A6CEE3','2025'='#B2DF8A'
)
# Original MiFish source: GG2_PCOA_12s_MIFISH_CLEAN_13052026.Rmd
mifish_functional_group_map <- c(
  "Clupea_harengus_Clupea_pallasii"                                                    = "Pelagic",
  "Myoxocephalus_scorpius_Myoxocephalus_quadricornis"                                  = "Benthic",
  "Gadus_morhua"                                                                       = "Demersal",
  "Pleuronectes_platessa_Platichthys_flesus"                                           = "Benthic",
  "Symphodus_melops"                                                                   = "Demersal",
  "no identification"                                                                  = NA,
  "Pollachius_virens_Pollachius_pollachius"                                            = "Demersal",
  "Gasterosteus_aculeatus"                                                             = "Pelagic",
  "Clupea_pallasii_Sprattus_sprattus_Clupea_harengus"                                 = "Pelagic",
  "Taurulus_bubalis"                                                                   = "Benthic",
  "Zoarces_viviparus"                                                                  = "Benthic",
  "Sprattus_sprattus_Clupea_pallasii"                                                  = "Pelagic",
  "Sprattus_sprattus_Clupea_harengus"                                                  = "Pelagic",
  "Belone_belone"                                                                      = "Pelagic",
  "Pomatoschistus_minutus"                                                             = "Benthic",
  "Platichthys_flesus_Pleuronectes_platessa_Hippoglossoides_platessoides"              = "Benthic",
  "Clupea_harengus"                                                                    = "Pelagic",
  "Limanda_limanda_Hippoglossus_hippoglossus"                                          = "Benthic",
  "Labrus_bergylta"                                                                    = "Demersal",
  "Anguilla_anguilla"                                                                  = "Benthic",
  "Merluccius_merluccius"                                                              = "Demersal",
  "Pholis_gunnellus"                                                                   = "Benthic",
  "Gobius_niger"                                                                       = "Benthic",
  "Scomber_scombrus"                                                                   = "Pelagic",
  "Hyperoplus_lanceolatus_Ammodytes_tobianus_Ammodytes_marinus_Hyperoplus_immaculatus" = "Pelagic",
  "Merlangius_merlangus"                                                               = "Demersal",
  "Gobiusculus_flavescens"                                                             = "Demersal",
  "Scophthalmus_rhombus"                                                               = "Benthic",
  "Salmo_trutta"                                                                       = "Pelagic",
  "Lumpenus_lampretaeformis_Lumpenus_fabricii"                                         = "Benthic",
  "Ammodytes_marinus_Hyperoplus_lanceolatus_Ammodytes_tobianus"                        = "Pelagic",
  "Lumpenus_lampretaeformis"                                                           = "Benthic",
  "Solea_solea"                                                                        = "Benthic",
  "Spinachia_spinachia"                                                                = "Pelagic",
  "Pleuronectes_platessa"                                                              = "Benthic",
  "Pollachius_pollachius_Pollachius_virens"                                            = "Demersal",
  "Raniceps_raninus"                                                                   = "Benthic",
  "Enchelyopus_cimbrius"                                                               = "Benthic",
  "Molva_molva"                                                                        = "Demersal",
  "Leptoclinus_maculatus"                                                              = "Benthic",
  "Trisopterus_minutus"                                                                = "Demersal",
  "Argentina"                                                                          = "Demersal",
  "Argentina_silus"                                                                    = "Demersal",
  "Sprattus_sprattus"                                                                  = "Pelagic"
)
# COI source: previously audited top-taxon original functional-group lookup
# (25_COI_diet_FOO_RRA_functional_reconstruction_FIXED.R). This is a
# restricted explicit assignment, not a blanket fish/invertebrate guess.
coi_functional_group_map <- c(
  "Clupea harengus" = "Pelagic",
  "Sprattus sprattus" = "Pelagic",
  "Micromesistius poutassou fish environmental sample" = "Pelagic",
  "Pollachius virens Pollachius environmental sample Actinopterygii environmental sample" = "Demersal",
  "Merlangius merlangus" = "Demersal",
  "Ctenolabrus rupestris" = "Demersal",
  "Pleuronectes platessa Pleuronectes sp. SH-2018" = "Benthic",
  "Myoxocephalus scorpius Myoxocephalus scorpioides Myoxocephalus aenaeus Myoxocephalus verrucosus" = "Benthic",
  "Pomatoschistus minutus Pomatoschistus microps Neogobius melanostomus Gobius niger" = "Benthic"
)

# The original PCoA figure is a FOUR-PANEL figure. Individual panels are
# diagnostic files; only the combined plot is promoted to outputs/figures.
if (!requireNamespace('patchwork',quietly=TRUE)) {
  stop('Install the patchwork package for the original four-panel layout')
}
# Run envfit on scores and EXACT sample x ASV matrix used for PCoA.
# Only mapped taxa contribute to a group. Do not rescale relative abundances.
get_functional_arrows <- function(marker, ps, mm, pcoa) {
  tx <- as.data.frame(tax_table(ps),check.names=FALSE)
  stopifnot('LCA_taxon' %in% names(tx))
  stopifnot(all(colnames(mm) %in% rownames(tx)))
  lookup <- if(marker=='MiFish') mifish_functional_group_map else coi_functional_group_map
  groups <- unname(lookup[as.character(tx[colnames(mm),'LCA_taxon'])])
  keep <- !is.na(groups) & groups %in% c('Pelagic','Demersal','Benthic')
  covered <- sum(keep)
  if(!covered) {
    warning('No taxonomically mapped functional groups for ',marker)
    return(data.frame(Dim1=numeric(),Dim2=numeric(),group=character(),p=numeric()))
  }
  dat <- data.frame(Pelagic=rowSums(mm[,keep & groups=='Pelagic',drop=FALSE]),
                    Demersal=rowSums(mm[,keep & groups=='Demersal',drop=FALSE]),
                    Benthic=rowSums(mm[,keep & groups=='Benthic',drop=FALSE]))
  # Exclude unvarying groups from envfit. This is needed for some predator subsets.
  varying <- vapply(dat,function(z) is.finite(stats::sd(z)) && stats::sd(z)>0,logical(1))
  dat <- dat[,varying,drop=FALSE]
  if(!ncol(dat)) return(data.frame(Dim1=numeric(),Dim2=numeric(),group=character(),p=numeric()))
  set.seed(42)
  ef <- vegan::envfit(pcoa,dat,permutations=999)
  vec <- as.data.frame(vegan::scores(ef,display='vectors'))
  if(!nrow(vec)) return(data.frame(Dim1=numeric(),Dim2=numeric(),group=character(),p=numeric()))
  vec <- data.frame(Dim1=vec[,1],Dim2=vec[,2],
                    group=rownames(vec),p=as.numeric(ef$vectors$pvals),
                    stringsAsFactors=FALSE)
  # Tarjei: only p <= 0.05 and 10x arrow scale.
  sig <- vec[is.finite(vec$p) & vec$p<=0.05,,drop=FALSE]
  sig$Dim1 <- 10 * sig$Dim1
  sig$Dim2 <- 10 * sig$Dim2
  attr(sig,'all_results') <- vec
  attr(sig,'mapped_ASVs') <- covered
  sig
}


fits <- list(); audit <- list(); perms <- list(); dispersions <- list(); panels <- list()
plain_md <- function(ps) {
  z <- methods::as(sample_data(ps),'data.frame')
  z <- base::data.frame(z,check.names=FALSE)
  rownames(z) <- sample_names(ps)
  z
}
mat_of <- function(ps) {
  z <- as(otu_table(ps),'matrix')
  if(taxa_are_rows(ps)) z <- t(z)
  storage.mode(z) <- 'double'
  z
}
perm_export <- function(fit,marker,predator,model,n) {
  z <- as.data.frame(fit); z$term <- rownames(z)
  data.frame(marker=marker,predator=predator,model=model,n=n,
             term=z$term,df=z$Df,SumOfSqs=z$SumOfSqs,R2=z$R2,
             F=z$F,p=z[['Pr(>F)']],row.names=NULL,check.names=FALSE)
}
run_one <- function(source,marker,predator,label,expected_n=NA_integer_) {
  md_all <- plain_md(source)
  if(!all(c('species','season','sample_year') %in% names(md_all)))
    stop('Missing script-02 metadata in ',marker)
  species <- as.character(md_all$species)
  species[species=='top_cormorant'] <- 'cormorant'
  ids <- rownames(md_all)[!is.na(species) & species==predator]
  if(!length(ids)) stop('No samples for ',marker,' ',predator)
  ps <- prune_samples(ids,source)
  md <- plain_md(ps)
  md$season <- as.character(md$season)
  md$season[is.na(md$season) | md$season=='Unknown' | !nzchar(md$season)] <- NA_character_
  md$ocean <- if('ocean' %in% names(md)) as.character(md$ocean) else NA_character_
  md$ocean[md$ocean=='Unknown'] <- NA_character_
  md$sample_year <- suppressWarnings(as.numeric(as.character(md$sample_year)))
  # Exactly the final original Rmd blocks: known season, season_year n>=2,
  # prune absent taxa and then zero-abundance samples.
  keep <- !is.na(md$season)
  ps <- prune_samples(rownames(md)[keep],ps); md <- md[keep,,drop=FALSE]
  md$season_year <- paste(md$season,md$sample_year,sep='_')
  groups <- table(md$season_year)
  keep <- md$season_year %in% names(groups)[groups>=2]
  ps <- prune_samples(rownames(md)[keep],ps); md <- md[keep,,drop=FALSE]
  ps <- prune_taxa(taxa_sums(ps)>0,ps)
  mm <- mat_of(ps)
  nonempty <- rowSums(mm)>0
  ps <- prune_samples(rownames(mm)[nonempty],ps)
  mm <- mat_of(ps); md <- md[rownames(mm),,drop=FALSE]
  if(!identical(rownames(mm),rownames(md))) stop('Sample alignment failure')
  if(anyNA(mm) || any(mm<0) || nrow(mm)<6 || ncol(mm)<2) stop('Bad OTU matrix')
  n <- nrow(mm)
  cat('\n==== ',marker,' | ',label,' | n=',n,' ASVs=',ncol(mm),' ====\n',sep='')
  status <- if(is.na(expected_n)) 'NOT_BENCHMARKED' else if(n==expected_n) 'MATCH' else 'REVIEW'
  if(marker=='MiFish' && status!='MATCH') stop('MiFish Table 3.4 sample mismatch for ',predator)
  if(marker=='COI' && status=='REVIEW') warning('COI thesis PCoA n=',expected_n,
     ' but source-code reconstruction n=',n,' for ',predator,
     '. DO NOT force-exclude biological samples.')
  d <- vegan::vegdist(mm,method='robust.aitchison')
  if(any(!is.finite(as.vector(d)))) stop('Nonfinite robust Aitchison distances')
  # TWO SOURCE VARIANTS ARE KEPT DISTINCT:
  # (1) 0.1.PERMANOVA.PCOA.Rmd / COI.PERMANOVA.Rmd test dispersion by
  #     SEASON ONLY. Retain this existing inferential result for continuity.
  # (2) GG2_PCOA_12s_MIFISH_CLEAN_13052026.Rmd plots the dispersion from
  #     SEASON-YEAR GROUPS, with year and sample n on the x-axis.
  #     Use that source's visual geometry for all four requested figures.
  bd <- vegan::betadisper(d, md$season)
  set.seed(42)
  bt <- vegan::permutest(bd, permutations=999)
  cat('\nOriginal season-only betadisper (inferential comparison):\n')
  print(bt)

  season_year_groups <- factor(md$season_year)
  if (anyNA(season_year_groups)) stop('Missing season-year group in dispersion')
  bd_sy <- vegan::betadisper(d, season_year_groups)
  set.seed(42)
  bt_sy <- vegan::permutest(bd_sy, permutations=999)
  cat('\nSeason-year betadisper (source of displayed boxplots):\n')
  print(bt_sy)

  dd <- data.frame(distance=as.numeric(bd_sy$distances),
                   group=as.character(bd_sy$group), stringsAsFactors=FALSE)
  # Reproduce original labels: "Autumn_2023\n(n=5)" and order by group.
  group_n <- table(dd$group)
  levels_group <- sort(names(group_n))
  dd$group <- factor(dd$group, levels=levels_group)
  dd$label <- factor(paste0(as.character(dd$group), '\n(n=',
                            as.integer(group_n[as.character(dd$group)]), ')'),
                     levels=paste0(levels_group, '\n(n=',
                                   as.integer(group_n[levels_group]), ')'))
  dd$season <- sub('_.*$', '', as.character(dd$group))
  if (!all(dd$season %in% names(season_colors)))
    stop('Unexpected season label in season-year betadisper')

  p_disp <- ggplot(dd, aes(x=label, y=distance, fill=season)) +
    geom_boxplot(color='black', outlier.shape=1) +
    scale_fill_manual(values=season_colors, name='Season',
                      breaks=c('Autumn','Spring','Summer','Winter')) +
    theme_classic(base_size=12) +
    labs(title=paste0('Multivariate Dispersion by Season-Year — ',label,
                     if(marker=='COI') ' (COI Leray-XT)' else ' (MiFish)'),
         x=NULL, y='Distance to Centroid', fill='Season') +
    theme(axis.text.x=element_text(angle=45,hjust=1,size=9),
          axis.text.y=element_text(size=10),
          axis.title.y=element_text(size=11),
          plot.title=element_text(hjust=.5,face='bold'),
          legend.position='right')
  suffix <- tolower(gsub(' ','_',label,fixed=TRUE))
  disp_name <- if(marker=='COI') paste0('betadisper_coi_',suffix,'.png') else paste0('betadisper_',suffix,'.png')
  ggsave(file.path(FIG,disp_name),p_disp,width=12,height=6,dpi=300,bg='white')
  # Keep season-year diagnostic separate from the season-only inferential test.
  sy_tab <- as.data.frame(bt_sy$tab)
  write.csv(data.frame(marker=marker,predator=predator,
                       season_year=levels_group,
                       n=as.integer(group_n[levels_group]),
                       label=paste0(levels_group,' (n=',as.integer(group_n[levels_group]),')')),
            file.path(VAL,paste0('03_BETADISPER_GROUPS_',marker,'_',predator,'.csv')),
            row.names=FALSE)
  write.csv(data.frame(marker=marker,predator=predator,
                       grouping='season_year',F=sy_tab[['F']][1],
                       p=sy_tab[['Pr(>F)']][1],n=n),
            file.path(VAL,paste0('03_BETADISPER_SEASON_YEAR_TEST_',marker,'_',predator,'.csv')),
            row.names=FALSE)
  # Original adonis2: year numeric, by='term', 999 permutations, seed 42.
  complete <- stats::complete.cases(md[,c('season','sample_year'),drop=FALSE])
  mperm <- base::data.frame(season=factor(md$season[complete]),
    sample_year=md$sample_year[complete],row.names=rownames(md)[complete])
  dm <- as.matrix(d)
  dperm <- stats::as.dist(dm[complete,complete,drop=FALSE])
  stopifnot(identical(rownames(mperm),attr(dperm,'Labels')))
  set.seed(42)
  pa <- vegan::adonis2(dperm~season+sample_year,data=mperm,by='term',permutations=999)
  set.seed(42)
  pi <- vegan::adonis2(dperm~season*sample_year,data=mperm,by='term',permutations=999)
  cat('\nAdditive PERMANOVA:\n'); print(pa)
  cat('\nInteraction PERMANOVA:\n'); print(pi)
  perms[[length(perms)+1L]] <<- perm_export(pa,marker,predator,'additive',nrow(mperm))
  perms[[length(perms)+1L]] <<- perm_export(pi,marker,predator,'interaction',nrow(mperm))
  bt_df <- as.data.frame(bt$tab)
  dispersions[[length(dispersions)+1L]] <<- data.frame(marker=marker,predator=predator,
    F=bt_df[['F']][1],p=bt_df[['Pr(>F)']][1],n=n)
  # Original base cmdscale and positive-eigenvalue variance definition.
  pc <- stats::cmdscale(d,k=2,eig=TRUE)
  eigen <- pc$eig
  var_exp <- round(eigen/sum(eigen[eigen>0])*100,1)
  xy <- pc$points
  pcoa_df <- data.frame(sample_id=rownames(xy),PC1=xy[,1],PC2=xy[,2],
    season=md$season,sample_year=md$sample_year,ocean=md$ocean,
    stringsAsFactors=FALSE)
  stopifnot(identical(pcoa_df$sample_id,rownames(md)))
  pcoa_df$year_period <- ifelse(pcoa_df$sample_year %in% 2012:2016,'Early (2012-2016)',
    ifelse(pcoa_df$sample_year %in% 2019:2021,'Mid (2019-2021)',
    ifelse(pcoa_df$sample_year %in% 2022:2025,'Recent (2022-2025)',NA_character_)))
  # Plot styling follows original GG2_PCOA source: year-coloured points,
  # period-coloured ellipses, season shapes and significant envfit vectors.
  pcoa_df$sample_year <- factor(as.character(pcoa_df$sample_year))
  year_levels <- levels(pcoa_df$sample_year)
  if(!all(year_levels %in% names(year_colors))) {
    stop('Missing original year colours: ',paste(setdiff(year_levels,names(year_colors)),collapse=', '))
  }
  pcoa_df$year_period <- factor(pcoa_df$year_period,
    levels=c('Early (2012-2016)','Mid (2019-2021)','Recent (2022-2025)'))
  if(anyNA(pcoa_df$year_period)) warning('Some PCoA samples have no time-period ellipse assignment')
  # Match taxonomy to the current predator-specific ASV subset.
  fg_ps <- prune_taxa(colnames(mm),ps)
  arrows <- get_functional_arrows(marker,fg_ps,mm,pc)
  arrow_trace <- attr(arrows,'all_results')
  if(is.null(arrow_trace)) arrow_trace <- data.frame(Dim1=numeric(),Dim2=numeric(),group=character(),p=numeric())
  write.csv(arrow_trace,file.path(VAL,paste0('03_ENVFIT_FUNCTIONAL_',marker,'_',predator,'.csv')),
    row.names=FALSE)
  cat('Envfit significant functional arrows for ',marker,' ',predator,': ',
    if(nrow(arrows)) paste(arrows$group,collapse=', ') else '(none)', '\n',sep='')
  p_pcoa <- ggplot(pcoa_df,aes(x=PC1,y=PC2))+
    stat_ellipse(aes(group=year_period,color=year_period),geom='polygon',fill=NA,
                 linewidth=1.2,linetype=2,show.legend=FALSE)+
    stat_ellipse(aes(group=year_period,fill=year_period),geom='polygon',alpha=.06,
                 linewidth=0,linetype=0,color=NA)+
    geom_point(aes(color=sample_year,shape=season),size=4,alpha=.95)+
    geom_segment(data=arrows,aes(x=0,y=0,xend=Dim1,yend=Dim2),
                 inherit.aes=FALSE,
                 arrow=grid::arrow(length=grid::unit(.4,'cm'),type='closed'),
                 color='black',linewidth=.7)+
    geom_label(data=arrows,aes(x=Dim1*1.1,y=Dim2*1.1,label=group),
               inherit.aes=FALSE,size=3.5,fontface='bold',fill='white',
               color='black',alpha=.9)+
    scale_color_manual(values=c(year_colors,period_colors),breaks=year_levels,
      name='Year')+
    scale_fill_manual(values=period_colors,name='Time Period (Ellipses)',drop=TRUE)+
    scale_shape_manual(values=season_shapes,name='Season')+
    guides(color=guide_legend(title='Year',order=1,
                 override.aes=list(shape=16,size=3)),
           fill=guide_legend(title='Time Period (Ellipses)',order=2,
                 override.aes=list(alpha=.4,size=5)),
           shape=guide_legend(title='Season',order=3))+
    theme_classic(base_size=10)+
    labs(title=paste0(if(marker=='MiFish') {
      if(predator=='cormorant') 'A' else 'B'
    } else if(predator=='cormorant') 'C' else 'D',
      '  ',label,', ',if(marker=='MiFish') 'MiFish' else 'COI Leray-XT'),
      x=paste0('PC1 (',var_exp[1],'% variance explained)'),
      y=paste0('PC2 (',var_exp[2],'% variance explained)'))+
    theme(plot.title=element_text(hjust=0,face='bold',size=10),
          legend.position='right',legend.key.size=grid::unit(.8,'lines'))
  key <- paste(marker,predator,sep='_')
  panels[[key]] <<- p_pcoa
  fig_name <- paste0('PCoA_',tolower(gsub(' ','_',label,fixed=TRUE)),
    if(marker=='MiFish') '_mifish.png' else '_coi.png')
  # Keep standalone panels strictly in validation; manuscript output is composite.
  ggsave(file.path(VAL,fig_name),p_pcoa,width=12,height=6,dpi=300,bg='white')
  audit[[length(audit)+1L]] <<- data.frame(marker=marker,predator=predator,
    n_initial=length(ids),n_final=n,n_taxa=ncol(mm),thesis_PCoA_n=expected_n,
    count_status=status,PC1_percent=var_exp[1],PC2_percent=var_exp[2],
    n_unknown_year_period=sum(is.na(pcoa_df$year_period)),
    figure_PCoA=fig_name,figure_betadisper=disp_name)
  write.csv(data.frame(sample_id=rownames(mm),season=md$season,year=md$sample_year,
            row.names=NULL),file.path(VAL,paste0('03_SAMPLE_TRACE_',marker,'_',predator,'.csv')),
            row.names=FALSE)
  write.csv(pcoa_df,file.path(VAL,paste0('03_PCOA_COORDS_',marker,'_',predator,'.csv')),row.names=FALSE)
  invisible(list(permanova=pa,interaction=pi,betadisper=bt,pcoa=pc))
}
# 1. MiFish expected n verified against thesis Table 3.4.
fits$MiFish_cormorant <- run_one(ps_mifish,'MiFish','cormorant','Phalacrocorax carbo',68L)
fits$MiFish_seal <- run_one(ps_mifish,'MiFish','harbor_seal','Phoca vitulina',59L)
# 2. COI: original cormorant PCoA caption n=64; seal caption n=56.
# Prior reconstruction observed seal n=60. Treat as review, not forced match.
fits$COI_cormorant <- run_one(ps_coi,'COI','cormorant','Phalacrocorax carbo',64L)
fits$COI_seal <- run_one(ps_coi,'COI','harbor_seal','Phoca vitulina',56L)
# Four-panel manuscript layout: A/B MiFish, C/D COI, in original order.
required_panels <- c('MiFish_cormorant','MiFish_harbor_seal',
                     'COI_cormorant','COI_harbor_seal')
stopifnot(all(required_panels %in% names(panels)))
combined <- (panels[['MiFish_cormorant']] + panels[['MiFish_harbor_seal']]) /
            (panels[['COI_cormorant']] + panels[['COI_harbor_seal']])
combined_name <- 'PCoA_diet_composition_four_panels.png'
ggsave(file.path(FIG,combined_name),combined,width=18,height=11,dpi=300,bg='white')
cat('Combined PCoA: ',file.path(FIG,combined_name),'\n',sep='')
# Old, superseded standalone plot exports remain in validation (above) only.
# Explicitly remove only the four filenames owned by prior script-03 runs.
old_pcoa <- c('PCoA_phalacrocorax_carbo_mifish.png',
              'PCoA_phoca_vitulina_mifish.png',
              'PCoA_phalacrocorax_carbo_coi.png',
              'PCoA_phoca_vitulina_coi.png')
unlink(file.path(FIG,old_pcoa))
# Analysis tables are available for direct comparison and potential manuscript
# use, but extra audit exports never enter manuscript output/figures.
permanent <- do.call(rbind,perms)
write.csv(permanent,file.path(OUT,'PERMANOVA_results.csv'),row.names=FALSE,na='')
write.csv(do.call(rbind,dispersions),file.path(VAL,'04_BETADISPER_RESULTS.csv'),row.names=FALSE)
write.csv(do.call(rbind,audit),file.path(VAL,'04_SAMPLE_AND_PCOA_AUDIT.csv'),row.names=FALSE)
write.csv(data.frame(input=names(paths),source=paths,md5=unname(tools::md5sum(paths))),
          file.path(VAL,'04_INPUT_HASHES.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(VAL,'04_SESSION_INFO.txt'))
writeLines(c('Original sources: 0.1.PERMANOVA.PCOA.Rmd (final ps_01 block) and COI.PERMANOVA.Rmd.',
   'Original method: season-year groups n>=2; robust Aitchison; adonis2 by=term, 999 permutations; cmdscale.',
   'Season-only betadisper statistics preserved from PERMANOVA scripts; season-year betadisper is used for plots following GG2_PCOA original.',
   'COI PCoA thesis n=56 seals versus reconstructed n=60 (if confirmed) is an OPEN discrepancy; no samples silently removed.',
   'Betadisper season-year figures accepted by user. PCoA year colours, ellipses, season symbols and original MiFish functional envfit map restored. COI limited functional mapping is reconstructed from prior original top-taxon diet lookup; envfit requires visual approval. Four panels combined in a single output figure; individual panel plots stored in validation.',
   'Only script-02 generated RDS objects were used as input. No old numbered reconstruction scripts are sourced.',
   'Figures require manual visual verification against original thesis/manuscript before approval.'),
   file.path(VAL,'04_METHOD_AND_LIMITATIONS.txt'))
cat('\n03 COMPLETE\nOriginal-figure candidates: ',FIG,'\nPERMANOVA table: ',file.path(OUT,'PERMANOVA_results.csv'),
    '\nValidation: ',VAL,'\n',sep='')
