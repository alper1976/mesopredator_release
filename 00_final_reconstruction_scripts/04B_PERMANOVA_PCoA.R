# PORTABLE STANDALONE RECONSTRUCTION. Original calculations retained.
# 04B_PERMANOVA_PCoA.R
# Independent PERMANOVA sensitivity / assumption audit for approved script 03.
# Never changes script 03, its figures, its input RDS files, or manuscript results.
# Reference: user-supplied 04_PERMANOVA_PCoA.R (03 reconstruction).
# Inputs: approved script 02 RDS files (same cohorts and filtering as script 03).
# Sensitivities: sequential order, additive marginal effects, overall effect,
# interaction order, dispersion, group sizes, season-year/region confounding,
# and optional within-year restricted season permutations.

needed <- c('phyloseq', 'vegan', 'permute')
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop('Install missing packages: ', paste(missing, collapse=', '))
suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
})

# --- 1. Same project structure as 03 ---------------------------------------
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
INPUT <- file.path(BASE, 'processed_data', '03_Phyloseq_Preprocessing')
OUT <- file.path(BASE, 'validation', '04B_PERMANOVA_PCoA')
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
PLOTS <- file.path(OUT, 'figures')
dir.create(PLOTS, recursive=TRUE, showWarnings=FALSE)

inputs <- c(MiFish='ps_mifish_global_0_1pct.rds',
            COI='ps_coi_global_0_1pct_identified.rds')
paths <- file.path(INPUT, inputs)
names(paths) <- names(inputs)
if (any(!file.exists(paths))) stop('Missing script-02 files: ', paste(paths[!file.exists(paths)], collapse='; '))

# --- 2. Utility functions --------------------------------------------------
plain_md <- function(ps) {
  z <- as.data.frame(phyloseq::sample_data(ps), check.names=FALSE)
  z[phyloseq::sample_names(ps), , drop=FALSE]
}
mat_of <- function(ps) {
  z <- as(phyloseq::otu_table(ps), 'matrix')
  if (phyloseq::taxa_are_rows(ps)) z <- t(z)
  storage.mode(z) <- 'double'
  z
}
write_table <- function(x, name) {
  write.csv(x, file.path(OUT, name), row.names=FALSE, na='')
}
extract_fit <- function(fit, marker, predator, test_name, n) {
  z <- as.data.frame(fit)
  data.frame(marker=marker, predator=predator, n=n, test=test_name,
             term=rownames(z), df=as.numeric(z$Df),
             SumOfSqs=as.numeric(z$SumOfSqs),
             R2=as.numeric(z$R2), F=as.numeric(z$F),
             p=as.numeric(z[['Pr(>F)']]), row.names=NULL,
             stringsAsFactors=FALSE)
}
get_term <- function(z, name) {
  x <- z[z$term == name, , drop=FALSE]
  if (nrow(x) != 1L) return(data.frame(R2=NA_real_,F=NA_real_,p=NA_real_))
  x[, c('R2','F','p'), drop=FALSE]
}
safe_disp <- function(d, g, marker, predator, grouping, seed) {
  g <- droplevels(factor(g))
  if (nlevels(g) < 2L || any(table(g) < 2L)) {
    return(data.frame(marker=marker,predator=predator,grouping=grouping,
                      n_groups=nlevels(g),min_group_n=min(table(g)),
                      F=NA_real_,p=NA_real_,status='NOT_TESTED_SINGLETON_OR_ONE_GROUP'))
  }
  set.seed(seed)
  fit <- vegan::betadisper(d, g, type='median', bias.adjust=FALSE)
  test <- vegan::permutest(fit, permutations=999)
  tab <- as.data.frame(test$tab)
  data.frame(marker=marker,predator=predator,grouping=grouping,
             n_groups=nlevels(g),min_group_n=min(table(g)),
             F=as.numeric(tab[['F']][1]),p=as.numeric(tab[['Pr(>F)']][1]),
             status='TESTED',stringsAsFactors=FALSE)
}
all_tests <- list(); all_diffs <- list(); all_disps <- list()
all_groups <- list(); all_samples <- list(); all_region <- list()
all_restrictions <- list()

run_audit <- function(source, marker, predator, expected_n) {
  # Reproduce selection order EXACTLY as script 03, including keeping
  # season-year groups n>=2 before dropping zero-abundance samples.
  md_all <- plain_md(source)
  stopifnot(all(c('species','season','sample_year') %in% names(md_all)))
  sp <- as.character(md_all$species)
  sp[sp=='top_cormorant'] <- 'cormorant'
  ids <- rownames(md_all)[!is.na(sp) & sp == predator]
  if (!length(ids)) stop('No samples found for ',marker,' ',predator)
  ps <- phyloseq::prune_samples(ids, source)
  md <- plain_md(ps)
  md$season <- as.character(md$season)
  md$season[is.na(md$season) | md$season=='Unknown' | !nzchar(md$season)] <- NA_character_
  md$sample_year <- suppressWarnings(as.numeric(as.character(md$sample_year)))
  keep <- !is.na(md$season)
  ps <- phyloseq::prune_samples(rownames(md)[keep], ps)
  md <- md[keep, , drop=FALSE]
  md$season_year <- paste(md$season, md$sample_year, sep='_')
  groups <- table(md$season_year)
  keep <- md$season_year %in% names(groups)[groups >= 2L]
  ps <- phyloseq::prune_samples(rownames(md)[keep], ps)
  md <- md[keep, , drop=FALSE]
  ps <- phyloseq::prune_taxa(phyloseq::taxa_sums(ps)>0, ps)
  mm <- mat_of(ps)
  ps <- phyloseq::prune_samples(rownames(mm)[rowSums(mm)>0], ps)
  mm <- mat_of(ps)
  md <- md[rownames(mm), , drop=FALSE]
  stopifnot(identical(rownames(mm),rownames(md)), all(is.finite(mm)), all(mm>=0))
  n <- nrow(mm)
  if (n != expected_n) {
    msg <- paste('Cohort difference for',marker,predator,':',n,'versus expected',expected_n)
    if (marker=='MiFish') stop(msg) else warning(msg, '; continuing without exclusions')
  }
  if (n < 8 || ncol(mm)<2) stop('Too few data for PERMANOVA')
  dist <- vegan::vegdist(mm, method='robust.aitchison')
  stopifnot(all(is.finite(as.vector(dist))))
  complete <- complete.cases(md[,c('season','sample_year'),drop=FALSE])
  mperm <- data.frame(season=factor(md$season[complete]),
                      sample_year=md$sample_year[complete],
                      row.names=rownames(md)[complete])
  dm <- as.matrix(dist)
  dperm <- as.dist(dm[complete,complete,drop=FALSE])
  stopifnot(identical(rownames(mperm), attr(dperm,'Labels')))
  n_perm <- nrow(mperm)
  mperm$season <- droplevels(mperm$season)
  if (nlevels(mperm$season)<2L || length(unique(mperm$sample_year))<2L) {
    stop('Insufficient season or year diversity for ',marker,' ',predator)
  }
  cat('\n==== ',marker,' / ',predator,' / n=',n_perm,' ====\n',sep='')

  # Both additive term orders with the same starting RNG seed.
  # Reversing the model term order reallocates sequential sums of squares.
  seed <- 42L
  set.seed(seed)
  original <- vegan::adonis2(dperm ~ season + sample_year,
                             data=mperm, permutations=999, by='terms')
  set.seed(seed)
  reversed <- vegan::adonis2(dperm ~ sample_year + season,
                             data=mperm, permutations=999, by='terms')
  set.seed(seed)
  marginal <- vegan::adonis2(dperm ~ season + sample_year,
                             data=mperm, permutations=999, by='margin')
  set.seed(seed)
  overall <- vegan::adonis2(dperm ~ season + sample_year,
                            data=mperm, permutations=999, by=NULL)
  # Interaction: preserve the original reference order and test order sensitivity.
  set.seed(seed)
  interaction <- vegan::adonis2(dperm ~ season * sample_year,
                                data=mperm, permutations=999, by='terms')
  set.seed(seed)
  interaction_reversed <- vegan::adonis2(dperm ~ sample_year * season,
                                         data=mperm, permutations=999, by='terms')
  set.seed(seed)
  interaction_marginal <- vegan::adonis2(dperm ~ season * sample_year,
                                         data=mperm, permutations=999, by='margin')

  tests <- rbind(
    extract_fit(original,marker,predator,'additive_season_first_ORIGINAL',n_perm),
    extract_fit(reversed,marker,predator,'additive_year_first',n_perm),
    extract_fit(marginal,marker,predator,'additive_marginal',n_perm),
    extract_fit(overall,marker,predator,'additive_overall',n_perm),
    extract_fit(interaction,marker,predator,'interaction_season_first_ORIGINAL',n_perm),
    extract_fit(interaction_reversed,marker,predator,'interaction_year_first',n_perm),
    extract_fit(interaction_marginal,marker,predator,'interaction_marginal_HIGHEST_ORDER_ONLY',n_perm))
  all_tests[[length(all_tests)+1L]] <<- tests

  # Side-by-side comparison: same named effect across sequential orders.
  for (term_name in c('season','sample_year')) {
    a <- get_term(tests[tests$test=='additive_season_first_ORIGINAL',],term_name)
    b <- get_term(tests[tests$test=='additive_year_first',],term_name)
    c <- get_term(tests[tests$test=='additive_marginal',],term_name)
    all_diffs[[length(all_diffs)+1L]] <<- data.frame(
      marker=marker,predator=predator,term=term_name,n=n_perm,
      R2_original=a$R2,R2_reversed=b$R2,R2_marginal=c$R2,
      delta_R2=b$R2-a$R2,
      F_original=a$F,F_reversed=b$F,F_marginal=c$F,
      p_original=a$p,p_reversed=b$p,p_marginal=c$p,
      delta_p=b$p-a$p,
      sig_original=ifelse(is.na(a$p),NA,a$p<0.05),
      sig_reversed=ifelse(is.na(b$p),NA,b$p<0.05),
      sig_marginal=ifelse(is.na(c$p),NA,c$p<0.05))
  }

  # Group dispersion: season is tested as in script 03. Season-year is
  # exploratory: the script 03 selection starts with >=2 per season-year,
  # but zero-read removal can subsequently create a singleton.
  all_disps[[length(all_disps)+1L]] <<- safe_disp(
    dperm,mperm$season,marker,predator,'season',seed)
  md_perm <- md[rownames(mperm),,drop=FALSE]
  all_disps[[length(all_disps)+1L]] <<- safe_disp(
    dperm,md_perm$season_year,marker,predator,'season_year',seed)

  group_tbl <- as.data.frame(table(season=mperm$season,
                                   sample_year=mperm$sample_year),stringsAsFactors=FALSE)
  names(group_tbl)[3] <- 'n'
  group_tbl$marker <- marker; group_tbl$predator <- predator
  all_groups[[length(all_groups)+1L]] <<- group_tbl

  sample_tbl <- data.frame(marker=marker,predator=predator,
                           sample_id=rownames(mperm),season=as.character(mperm$season),
                           sample_year=mperm$sample_year,stringsAsFactors=FALSE)
  if ('ocean' %in% names(md_perm)) sample_tbl$ocean <- as.character(md_perm$ocean)
  all_samples[[length(all_samples)+1L]] <<- sample_tbl
  if ('ocean' %in% names(sample_tbl)) {
    region_tbl <- as.data.frame(table(sample_year=sample_tbl$sample_year,
                                     ocean=sample_tbl$ocean),stringsAsFactors=FALSE)
    names(region_tbl)[3] <- 'n'; region_tbl$marker <- marker
    region_tbl$predator <- predator
    all_region[[length(all_region)+1L]] <<- region_tbl
  }

  # With-year restricted permutations are only informative about season if
  # there is within-year seasonal variability. Do not report year p-values
  # from this restricted design: year is constant within blocks.
  year_blocks <- factor(mperm$sample_year)
  tab <- table(year_blocks, mperm$season)
  varying <- sum(rowSums(tab>0) >= 2L)
  eligible <- varying >= 2L && sum(table(year_blocks)>1L) >= 2L
  if (eligible) {
    ctl <- permute::how(nperm=999, blocks=year_blocks)
    restricted <- tryCatch({
      set.seed(seed)
      vegan::adonis2(dperm ~ sample_year + season,
                     data=mperm, by='terms', permutations=ctl)
    },error=function(e)e)
    if (!inherits(restricted,'error')) {
      rz <- extract_fit(restricted,marker,predator,'within_year_restricted_SEASON_ONLY',n_perm)
      rz <- rz[rz$term=='season',,drop=FALSE]
      all_restrictions[[length(all_restrictions)+1L]] <<- data.frame(
        marker=marker,predator=predator,n=n_perm,
        within_years_with_multiple_seasons=varying,
        test='season_conditioned_on_numeric_year',
        R2=rz$R2,F=rz$F,p=rz$p,
        status='EXPLORATORY_WITHIN_YEAR_PERMUTATIONS')
    } else {
      all_restrictions[[length(all_restrictions)+1L]] <<- data.frame(
        marker=marker,predator=predator,n=n_perm,
        within_years_with_multiple_seasons=varying,
        test='season_conditioned_on_numeric_year',R2=NA_real_,F=NA_real_,p=NA_real_,
        status=paste0('NOT_RUN: ',conditionMessage(restricted)))
    }
  } else {
    all_restrictions[[length(all_restrictions)+1L]] <<- data.frame(
      marker=marker,predator=predator,n=n_perm,
      within_years_with_multiple_seasons=varying,
      test='season_conditioned_on_numeric_year',R2=NA_real_,F=NA_real_,p=NA_real_,
      status='NOT_RUN_INSUFFICIENT_WITHIN_YEAR_SEASONAL_VARIATION')
  }
  invisible(tests)
}

# --- 3. Run the four APPROVED script-03 cohorts -----------------------------
ps_mifish <- readRDS(paths[['MiFish']])
ps_coi <- readRDS(paths[['COI']])
stopifnot(phyloseq::nsamples(ps_mifish)==161L,
          phyloseq::ntaxa(ps_mifish)==160L,
          phyloseq::nsamples(ps_coi)==145L,
          phyloseq::ntaxa(ps_coi)==172L)
tx <- as.data.frame(phyloseq::tax_table(ps_mifish),check.names=FALSE)
stopifnot('LCA_taxon' %in% names(tx))
no_id <- !is.na(tx$LCA_taxon) &
  tolower(trimws(as.character(tx$LCA_taxon)))=='no identification'
if (sum(no_id)!=91L) stop('Expected 91 unidentified MiFish ASVs as in script 03')
ps_mifish <- phyloseq::prune_taxa(rownames(tx)[!no_id],ps_mifish)
stopifnot(phyloseq::ntaxa(ps_mifish)==69L)
run_audit(ps_mifish,'MiFish','cormorant',68L)
run_audit(ps_mifish,'MiFish','harbor_seal',59L)
run_audit(ps_coi,'COI','cormorant',64L)
run_audit(ps_coi,'COI','harbor_seal',60L)  # approved reconstruction n=60; reported PCoA n=56

# --- 4. Export dedicated QC outputs (never write into script 03 outputs) ----
tests <- do.call(rbind,all_tests)
diffs <- do.call(rbind,all_diffs)
disps <- do.call(rbind,all_disps)
groups <- do.call(rbind,all_groups)
samples <- do.call(rbind,all_samples)
restrictions <- do.call(rbind,all_restrictions)
write_table(tests,'04B_ALL_PERMANOVA_TESTS.csv')
write_table(diffs,'04B_TERM_ORDER_SENSITIVITY.csv')
write_table(disps,'04B_BETADISPER_TESTS.csv')
write_table(groups,'04B_SEASON_BY_YEAR_COUNTS.csv')
write_table(samples,'04B_SAMPLE_TRACE.csv')
write_table(restrictions,'04B_WITHIN_YEAR_PERMUTATION_SENSITIVITY.csv')
if (length(all_region)) write_table(do.call(rbind,all_region),'04B_YEAR_BY_REGION_COUNTS.csv')

# Numeric equality check against the existing approved 03 export (where present).
original_csv <- file.path(BASE,'processed_data','04_PERMANOVA_PCoA','PERMANOVA_results.csv')
if (file.exists(original_csv)) {
  baseline <- read.csv(original_csv,stringsAsFactors=FALSE,check.names=FALSE)
  # Baseline uses names additive / interaction, unlike the QC test names.
  check <- subset(tests,test %in% c('additive_season_first_ORIGINAL',
                                   'interaction_season_first_ORIGINAL'))
  check$model <- ifelse(check$test=='additive_season_first_ORIGINAL','additive','interaction')
  need <- c('marker','predator','model','term','df','SumOfSqs','R2','F','p')
  if (all(need %in% names(baseline))) {
    baseline$key <- paste(baseline$marker,baseline$predator,baseline$model,baseline$term,sep='|')
    check$key <- paste(check$marker,check$predator,check$model,check$term,sep='|')
    ix <- match(check$key,baseline$key)
    parity <- data.frame(key=check$key,matched_row=!is.na(ix),stringsAsFactors=FALSE)
    for (field in c('df','SumOfSqs','R2','F','p')) {
      parity[[paste0(field,'_qc')]] <- check[[field]]
      parity[[paste0(field,'_script03')]] <- baseline[[field]][ix]
      parity[[paste0(field,'_abs_difference')]] <- abs(check[[field]]-baseline[[field]][ix])
    }
    write_table(parity,'04B_SCRIPT03_PARITY.csv')
    if (any(!parity$matched_row)) warning('Some script 03 rows could not be matched; inspect parity audit')
  } else warning('Script 03 baseline column schema differs; parity check skipped')
} else warning('Script 03 PERMANOVA_results.csv not found; original-model parity check skipped')

write_table(data.frame(input=names(paths),path=unname(paths),
                       md5=unname(tools::md5sum(paths))),
            '04B_INPUT_HASHES.csv')
writeLines(c(
  '03B: SEPARATE PERMANOVA SENSITIVITY / ASSUMPTIONS QC (NOT REPLACEMENT OF 03)',
  'Identical input files, original sample filters, robust Aitchison distance and original model formulas as script 03.',
  'Unrestricted permutations: 999; set.seed(42) before each adonis2 call.',
  'Original sequential model: season + numeric sample_year, by="terms".',
  'Reordered sequential model: numeric sample_year + season, by="terms".',
  'For nonorthogonal designs, sequential R2/F/p can depend on term order.',
  'Additive marginal test: each main term adjusted for the other main term.',
  'Marginal test on interaction model tests highest-order term only, not simple main effects.',
  'Season betadisper: median-centred, default settings as in 03; test differences in spread (not a complete assumption test).',
  'Season-year betadisper is exploratory and is skipped when any final group is singleton.',
  'Restricted year-block permutations: exploratory season test, only when season varies sufficiently within years.',
  'Do NOT interpret year p-values under within-year restricted permutations.',
  'Year is treated as a NUMERIC covariate; this assumes linearity of the temporal effect.',
  'Season-year- and year-region-imbalance diagnostics help assess exchangeability and confounding.',
  'Region and year confounding cannot be resolved by simply reordering formula terms.',
  'Unequal group dispersion can influence PERMANOVA; significant betadisper calls for interpretation caution.',
  'Multiple permutation tests are exploratory; p values have limited resolution with 999 permutations.',
  'COI seal cohort is n=60 from script 03 reconstruction; original PCoA description n=56 remains unresolved.',
  'The approved script 03 tables, plots and outputs remain unchanged.'
),file.path(OUT,'04B_METHOD_AND_INTERPRETATION.txt'))
capture.output(sessionInfo(),file=file.path(OUT,'04B_SESSION_INFO.txt'))
cat('\n03B QC COMPLETE\nOutput: ',OUT,
    '\n\nSEQUENTIAL ORDER DIFFERENCES:\n',sep='')
print(diffs[,c('marker','predator','term','R2_original','R2_reversed','p_original','p_reversed','p_marginal')],
      row.names=FALSE)
cat('\nBETADISPER:\n')
print(disps,row.names=FALSE)
cat('\nWITHIN-YEAR PERMUTATION SENSITIVITY:\n')
print(restrictions,row.names=FALSE)
