# PORTABLE STANDALONE RECONSTRUCTION. Original calculations retained.
# 06_GLLVM.R
# Mesopredator_check_2 — validated GLLVM, Tarjei-inspired plotting style
# Uses already verified 5%-prevalence MiFish Tweedie LV0 half-year models.
# DOES NOT refit Tarjei's failed original four-season GLLVM.
# Original style: two-by-two predator x effect coefficient panels; seasonal
# colour yellow (#e6ab02), significant year colours predator red/blue;
# nonsignificant grey, dashed zero line, Wald CI, italic taxon labels.
# IMPORTANT: The fitted half-year contrast is NOT equivalent to Tarjei's
# individual Summer/Autumn/Winter contrasts. Label the effect honestly.
# Wald intervals remain preliminary conditional on fit and assumptions.

required <- c('gllvm', 'ggplot2', 'patchwork')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Missing R packages: ', paste(missing, collapse=', '))
suppressPackageStartupMessages({
  library(gllvm)
  library(ggplot2)
  library(patchwork)
})

# Locate project root: run either from project root or final scripts directory.
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
SOURCE20 <- file.path(BASE, 'processed_data', 'reference_models', 'MiFish', '20_GLLVM_decision')
SOURCE21B <- file.path(BASE, 'processed_data', 'reference_models', 'MiFish', '21B_prediction_identity')
SOURCE_COI <- file.path(BASE, 'processed_data', 'reference_models', 'GLLVM', '29A_COI_original')
OUT <- file.path(BASE, 'processed_data', '06_GLLVM')
FIG <- file.path(OUT, 'figures')
TAB <- file.path(OUT, 'tables')
VAL <- file.path(BASE, 'validation', '06_GLLVM')
for(p in c(FIG, TAB, VAL)) dir.create(p, recursive=TRUE, showWarnings=FALSE)

ids <- c(cormorant='20_cormorant_p05_tweedie_halfyear_lv0',
         harbor_seal='20_harbor_seal_p05_tweedie_halfyear_lv0')
master_path <- file.path(SOURCE20, '20_GLLVM_MASTER.csv')
decision_path <- file.path(SOURCE20, '20_GLLVM_DECISION_TABLE.csv')
alignment_path <- file.path(SOURCE21B, '21B_PREDICTION_IDENTITY_DECISION.csv')
model_paths <- file.path(SOURCE20, paste0(unname(ids), '_fit.rds'))
names(model_paths) <- names(ids)
all_files <- c(master_path, decision_path, alignment_path, unname(model_paths))
if(any(!file.exists(all_files))) stop('Missing verified source file(s):\n',
  paste(all_files[!file.exists(all_files)], collapse='\n'))

master <- read.csv(master_path, stringsAsFactors=FALSE, check.names=FALSE)
decision <- read.csv(decision_path, stringsAsFactors=FALSE, check.names=FALSE)
alignment <- read.csv(alignment_path, stringsAsFactors=FALSE, check.names=FALSE)
need_master <- c('model_id','record_type','taxon','term','estimate','std_error','ci_low','ci_high')
need_dec <- c('model_id','convergence','CI_screen','n_samples','n_taxa',
              'prevalence','family','lv','season_mode')
if(!all(need_master %in% names(master)) || !all(need_dec %in% names(decision)) ||
   !all(c('model_id','alignment_status') %in% names(alignment)))
  stop('Validation files do not have expected schema; do not infer columns')

is_true <- function(x) length(x)==1L && !is.na(x) && toupper(as.character(x))=='TRUE'
checked <- vector('list', length(ids)); names(checked) <- names(ids)
coef_rows <- vector('list', length(ids)); names(coef_rows) <- names(ids)

for(pred in names(ids)) {
  id <- ids[[pred]]
  d <- decision[decision$model_id == id, , drop=FALSE]
  a <- alignment[alignment$model_id == id, , drop=FALSE]
  if(nrow(d)!=1L || nrow(a)!=1L) stop('Missing or duplicate validation record: ', id)
  if(!is_true(d$convergence[[1]]) || !is_true(d$CI_screen[[1]]))
    stop('Source model did not pass convergence/CI screen: ', id)
  if(!identical(a$alignment_status[[1]], 'VERIFIED_BY_COVARIATE_RECONSTRUCTION'))
    stop('Prediction identity not verified: ', id)
  if(as.numeric(d$prevalence[[1]])!=5 || as.character(d$family[[1]])!='tweedie' ||
     as.numeric(d$lv[[1]])!=0 || as.character(d$season_mode[[1]])!='halfyear')
    stop('Incorrect original validation specification: ', id)

  saved <- readRDS(model_paths[[pred]])
  if(!is.list(saved) || is.null(saved$model) || is.null(saved$spec) ||
     !inherits(saved$model,'gllvm') ||
     !identical(as.character(saved$spec$model_id), id))
    stop('Saved gllvm archive is not the validated model: ', id)
  m <- saved$model
  if(length(m$family)!=ncol(m$y) || !all(as.character(m$family)=='tweedie') ||
     !identical(as.numeric(m$num.lv),0) ||
     !identical(attr(terms(m$formula),'term.labels'),c('season_halfyear','year_scaled')) ||
     !all(c('season_halfyear','year_scaled') %in% names(m$X)))
    stop('Saved model formula/family mismatch: ', id)
  if(nrow(m$y)!=as.integer(d$n_samples[[1]]) ||
     ncol(m$y)!=as.integer(d$n_taxa[[1]]) ||
     !identical(rownames(m$y), rownames(m$X)) || is.null(colnames(m$y)))
    stop('Sample/taxon provenance mismatch: ', id)
  if(!isTRUE(m$convergence) || !is.finite(m$logL))
    stop('Archived model does not report finite converged fit: ', id)

  x <- master[master$model_id==id & master$record_type=='COEFFICIENT', , drop=FALSE]
  needed_terms <- c('(Intercept)', 'season_halfyearsummer_AprSep', 'year_scaled')
  if(nrow(x)!=ncol(m$y)*3L ||
     anyDuplicated(paste(x$taxon,x$term,sep='|')) ||
     !setequal(unique(x$term), needed_terms) ||
     !setequal(unique(x$taxon), colnames(m$y)))
    stop('Coefficient identities do not match validated model: ', id)
  for(v in c('estimate','std_error','ci_low','ci_high')) {
    x[[v]] <- suppressWarnings(as.numeric(x[[v]]))
    if(any(!is.finite(x[[v]]))) stop('Invalid ',v,' in ',id)
  }
  if(any(x$std_error <= 0 | x$ci_high <= x$ci_low))
    stop('Invalid saved Hessian intervals: ', id)

  # Independently compare CSV estimates with saved model coefficients.
  cf <- coef(m)
  if(is.null(cf$Xcoef) || is.null(cf$Intercept))
    stop('Fitted model coefficients inaccessible: ', id)
  B <- cbind('(Intercept)'=as.numeric(cf$Intercept), cf$Xcoef)
  rownames(B) <- rownames(cf$Xcoef)
  ix <- cbind(match(x$taxon,rownames(B)), match(x$term,colnames(B)))
  if(anyNA(ix) || max(abs(x$estimate - B[ix])) > 1e-5)
    stop('CSV coefficient estimate differs from archived fit: ', id)
  # Confirm the stored 95% Wald CI really uses the stated SE.
  if(max(abs(x$ci_low-(x$estimate-1.96*x$std_error)),
         abs(x$ci_high-(x$estimate+1.96*x$std_error))) > 0.005)
    stop('Wald CI does not match recorded coefficient and SE: ', id)

  # Verify covariate ordering used for checked predictions.
  D <- model.matrix(~ season_halfyear + year_scaled, data=m$X)
  if(!identical(rownames(D),rownames(m$y)) ||
     !identical(colnames(D),colnames(B)) ||
     !identical(rownames(B),colnames(m$y)))
    stop('Design matrix / coefficient ordering mismatch: ', id)
  mu <- exp(D %*% t(B))
  if(any(!is.finite(mu))) stop('Nonfinite model predictions: ',id)

  x$predator <- pred
  x$model_id <- id
  coef_rows[[pred]] <- x[,c('predator','model_id','taxon','term','estimate',
                             'std_error','ci_low','ci_high')]
  checked[[pred]] <- data.frame(predator=pred,model_id=id,
    n_samples=nrow(m$y),n_taxa=ncol(m$y),logL=m$logL,
    convergence=as.logical(m$convergence),
    RMSE_training=sqrt(mean((m$y-mu)^2)),
    MAE_training=mean(abs(m$y-mu)),
    stringsAsFactors=FALSE)
  cat('Verified: ',pred,' | n=',nrow(m$y),' taxa=',ncol(m$y),'\n',sep='')
}

coef_out <- do.call(rbind,coef_rows)
row.names(coef_out) <- NULL
write.csv(coef_out,file.path(TAB,'gllvm_validated_coefficients.csv'),row.names=FALSE)
quality <- do.call(rbind,checked)
write.csv(quality,file.path(VAL,'06_MODEL_VALIDATION_SUMMARY.csv'),row.names=FALSE)

# COI evidence is distinct from independently validated MiFish.
# Earlier 29A reconstruction retained coefficient rows for cormorant, but
# its convergence code was not recorded. Do not label COI as converged.
coi_coef_path <- file.path(SOURCE_COI, '29A_COI_COEFFICIENTS_cormorant.csv')
coi_diag_path <- file.path(SOURCE_COI, '29A_COI_MODEL_DIAGNOSTICS.csv')
if(!file.exists(coi_coef_path) || !file.exists(coi_diag_path))
  stop('COI 29A reconstruction outputs missing; cannot create COI panels.\n',
       paste(c(coi_coef_path,coi_diag_path), collapse='\n'))
coi_diag <- read.csv(coi_diag_path, stringsAsFactors=FALSE, check.names=FALSE)
if(!all(c('predator','fit_returned','error') %in% names(coi_diag)))
  stop('Unexpected COI diagnostics schema')
coi_corm_diag <- coi_diag[coi_diag$predator=='cormorant',,drop=FALSE]
if(nrow(coi_corm_diag)!=1L || !isTRUE(coi_corm_diag$fit_returned[[1]]) ||
   nzchar(coi_corm_diag$error[[1]])) stop('COI cormorant fit did not return cleanly')
coi_raw <- read.csv(coi_coef_path, stringsAsFactors=FALSE, check.names=FALSE)
coi_required <- c('predator','predictor','taxon','estimate','se','p','reliable')
if(!all(coi_required %in% names(coi_raw))) stop('Unexpected COI coefficient schema')
coi_raw <- coi_raw[coi_raw$predator=='cormorant',,drop=FALSE]
for(v in c('estimate','se','p')) coi_raw[[v]] <- as.numeric(coi_raw[[v]])
coi_raw$reliable <- as.character(coi_raw$reliable) %in% c('TRUE','T','1')
coi_plot <- coi_raw[coi_raw$reliable & is.finite(coi_raw$estimate) &
                     is.finite(coi_raw$se) & coi_raw$se>0 &
                     (grepl('^season',coi_raw$predictor) |
                      coi_raw$predictor=='year_scaled'),,drop=FALSE]
if(!nrow(coi_plot) || !any(grepl('^season',coi_plot$predictor)) ||
   !any(coi_plot$predictor=='year_scaled'))
  stop('No reliable COI seasonal and year coefficients in 29A')
coi_plot$ci_low <- coi_plot$estimate - 1.96*coi_plot$se
coi_plot$ci_high <- coi_plot$estimate + 1.96*coi_plot$se
coi_plot$significant <- coi_plot$ci_low>0 | coi_plot$ci_high<0
coi_plot$season <- sub('^season','',coi_plot$predictor)
write.csv(coi_plot, file.path(VAL,'06_COI_CORMORANT_SCREENED_COEFFICIENTS.csv'),row.names=FALSE)
write.csv(coi_diag, file.path(VAL,'06_COI_29A_ORIGINAL_DIAGNOSTICS.csv'),row.names=FALSE)

# -------------------------------------------------------------------------
# FIGURE STYLE: two independent legends, exactly one per marker group.
# MiFish half-year: coefficient is summer minus winter. A significant
# positive estimate is summer-associated; a significant negative estimate
# is winter-associated. Grey means the 95% CI spans zero.
# This sign-colouring does not imply that the model fitted two coefficients.
# COI: season contrasts (Summer/Autumn/Winter) all use Spring as reference.
# No model estimates or confidence intervals are changed by these settings.
# -------------------------------------------------------------------------

predator_names <- c(cormorant='Great cormorant', harbor_seal='Harbour seal')
predator_colors <- c(cormorant='#d6604d', harbor_seal='#2166ac')
mifish_half_cols <- c('Summer half-year'='#1b9e77',
                      'Winter half-year'='#377eb8',
                      'Not significant'='grey75')
coi_season_cols <- c('Summer'='#1b9e77', 'Autumn'='#d95f02',
                     'Winter'='#377eb8', 'Not significant'='grey75')

coef_out$effect <- ifelse(coef_out$term=='season_halfyearsummer_AprSep',
                          'Half-year', ifelse(coef_out$term=='year_scaled',
                          'Year','Intercept'))
coef_out$significant <- (coef_out$ci_low > 0 | coef_out$ci_high < 0)
coef_out$z <- coef_out$estimate / coef_out$std_error
coef_out$p_wald <- 2*pnorm(abs(coef_out$z),lower.tail=FALSE)
coef_out$stars <- ifelse(coef_out$p_wald<0.001,'***',
                  ifelse(coef_out$p_wald<0.01,'**',
                  ifelse(coef_out$p_wald<0.05,'*',
                  ifelse(coef_out$p_wald<0.1,'.',''))))
short_label <- function(x) vapply(strsplit(as.character(x),' ',fixed=TRUE),
  function(w) paste(head(w,2L),collapse=' '),character(1))
coef_out$taxon_short <- short_label(coef_out$taxon)

base_gllvm_theme <- function() theme_classic(base_size=10) +
  theme(axis.text.y=element_text(face='italic',size=8),
        plot.title=element_text(face='bold',size=10,hjust=.5),
        axis.line.y=element_blank(),axis.ticks.y=element_blank(),
        panel.grid.major.x=element_line(colour='grey90',linewidth=.3),
        legend.title=element_text(size=9),legend.text=element_text(size=8))

panel_mifish_halfyear <- function(pred,show_legend=FALSE) {
  a <- coef_out[coef_out$predator==pred & coef_out$effect=='Half-year',,drop=FALSE]
  if(!nrow(a)) stop('No MiFish half-year coefficients for ',pred)
  ord <- unique(a$taxon[order(a$estimate)])
  a$taxon_display <- factor(a$taxon,levels=ord,labels=make.unique(short_label(ord)))
  a$display_group <- ifelse(!a$significant,'Not significant',
                             ifelse(a$estimate>0,'Summer half-year','Winter half-year'))
  a$display_group <- factor(a$display_group,levels=names(mifish_half_cols))
  ggplot(a,aes(x=estimate,y=taxon_display)) +
    geom_vline(xintercept=0,linetype='dashed',colour='grey50',linewidth=.5) +
    geom_segment(aes(x=ci_low,xend=ci_high,y=taxon_display,yend=taxon_display,
                     colour=display_group),linewidth=.6,lineend='butt') +
    geom_point(aes(colour=display_group),size=2.5) +
    scale_colour_manual(values=mifish_half_cols,limits=names(mifish_half_cols),
      drop=FALSE,name='MiFish half-year',
      guide=if(show_legend) guide_legend(order=1) else 'none') +
    labs(x='Coefficient: summer (Apr-Sep) vs winter (Oct-Mar) (95% Wald CI)',y=NULL,
         title=paste0('Half-year effects — ',predator_names[[pred]],' (MiFish)')) +
    base_gllvm_theme()
}

panel_mifish_year <- function(pred) {
  a <- coef_out[coef_out$predator==pred & coef_out$effect=='Year',,drop=FALSE]
  if(!nrow(a)) stop('No MiFish year coefficients for ',pred)
  ord <- unique(a$taxon[order(a$estimate)])
  a$taxon_display <- factor(a$taxon,levels=ord,labels=make.unique(short_label(ord)))
  a$display_group <- ifelse(a$significant,'Year','Not significant')
  year_cols <- c('Year'=unname(predator_colors[[pred]]),'Not significant'='grey75')
  ggplot(a,aes(x=estimate,y=taxon_display)) +
    geom_vline(xintercept=0,linetype='dashed',colour='grey50',linewidth=.5) +
    geom_segment(aes(x=ci_low,xend=ci_high,y=taxon_display,yend=taxon_display,
                     colour=display_group),linewidth=.6,lineend='butt') +
    geom_point(aes(colour=display_group),size=2.5) +
    scale_colour_manual(values=year_cols,guide='none') +
    labs(x='Coefficient: sampling year (standardized) (95% Wald CI)',y=NULL,
         title=paste0('Year trend — ',predator_names[[pred]],' (MiFish)')) +
    base_gllvm_theme()
}

pA <- panel_mifish_halfyear('cormorant',show_legend=TRUE)
pB <- panel_mifish_year('cormorant')
pC <- panel_mifish_halfyear('harbor_seal',show_legend=FALSE)
pD <- panel_mifish_year('harbor_seal')

# COI cormorant: original four-season contrasts, Spring reference.
# Earlier archived COI reconstruction did not independently establish
# convergence; only screened coefficients are shown.
coi_short <- short_label
coi_season <- coi_plot[grepl('^season',coi_plot$predictor),,drop=FALSE]
coi_year <- coi_plot[coi_plot$predictor=='year_scaled',,drop=FALSE]
coi_season$season <- sub('^season','',coi_season$predictor)
if(anyNA(coi_season$season) ||
   !all(coi_season$season %in% c('Summer','Autumn','Winter')))
  stop('Unexpected COI seasonal term; inspect 29A coefficient file')
coi_season$display_group <- ifelse(coi_season$significant,
                                    coi_season$season,'Not significant')
coi_season$display_group <- factor(coi_season$display_group,
                                    levels=names(coi_season_cols))

# Explicit distinct numeric y positions and identical y/yend for horizontal CIs.
season_taxa <- sort(unique(coi_season$taxon))
coi_season$y_base <- match(coi_season$taxon,season_taxa)
y_offset <- c(Summer=.24,Autumn=0,Winter=-.24)
coi_season$y_pos <- coi_season$y_base + unname(y_offset[coi_season$season])
if(anyNA(coi_season$y_pos)) stop('COI season y positions contain NA')

pE <- ggplot(coi_season) +
  geom_vline(xintercept=0,linetype='dashed',colour='grey50',linewidth=.5) +
  geom_segment(aes(x=ci_low,xend=ci_high,y=y_pos,yend=y_pos,
                   colour=display_group),linewidth=.6,lineend='butt') +
  geom_point(aes(x=estimate,y=y_pos,colour=display_group),size=2.2) +
  scale_colour_manual(values=coi_season_cols,limits=names(coi_season_cols),
                      drop=FALSE,name='COI season',guide=guide_legend(order=2)) +
  scale_y_continuous(breaks=seq_along(season_taxa),labels=coi_short(season_taxa),
                     expand=expansion(add=.5)) +
  labs(x='Season coefficient vs spring (95% Wald CI)',y=NULL,
       title='E   Season effects — Great cormorant (COI Leray-XT)') +
  base_gllvm_theme()

year_taxa <- sort(unique(coi_year$taxon))
coi_year$y_pos <- match(coi_year$taxon,year_taxa)
coi_year$display_group <- ifelse(coi_year$significant,'Year','Not significant')
coi_year_cols <- c('Year'='#d6604d','Not significant'='grey75')
pF <- ggplot(coi_year) +
  geom_vline(xintercept=0,linetype='dashed',colour='grey50',linewidth=.5) +
  geom_segment(aes(x=ci_low,xend=ci_high,y=y_pos,yend=y_pos,
                   colour=display_group),linewidth=.6,lineend='butt') +
  geom_point(aes(x=estimate,y=y_pos,colour=display_group),size=2.4) +
  scale_colour_manual(values=coi_year_cols,guide='none') +
  scale_y_continuous(breaks=seq_along(year_taxa),labels=coi_short(year_taxa),
                     expand=expansion(add=.5)) +
  labs(x='Year coefficient (95% Wald CI)',y=NULL,
       title='F   Year trend — Great cormorant (COI Leray-XT)') +
  base_gllvm_theme()

write.csv(coi_season,file.path(VAL,'06_COI_SEASON_PLOT_POSITIONS.csv'),row.names=FALSE)
# One manuscript figure only, 3 rows by 2 columns.
combined <- (pA+pB)/(pC+pD)/(pE+pF) +
  patchwork::plot_layout(heights=c(1,1,1), guides='keep')
fig_path <- file.path(FIG,'gllvm_coefficients_combined.png')
ggsave(fig_path,plot=combined,width=17,height=18,dpi=300,bg='white')
# Retire the previous four-panel figure from manuscript outputs.
unlink(file.path(FIG,'gllvm_coef_mifish.png'))

# Plot input is preserved for independent cross-reference and plotting audits.
write.csv(coef_out,file.path(VAL,'06_PLOT_DATA_MIFISH.csv'),row.names=FALSE)
source_files <- c(master_path,decision_path,alignment_path,unname(model_paths),
                  coi_coef_path,coi_diag_path)
write.csv(data.frame(source_file=source_files,md5=unname(tools::md5sum(source_files))),
          file.path(VAL,'06_INPUT_PROVENANCE.csv'),row.names=FALSE)
writeLines(c(
  '05 GLLVM analysis: independently validated MiFish, prevalence >=5%,',
  'two half-years summer Apr-Sep vs winter Oct-Mar; Tweedie, LV0, scaled year.',
  'Models are loaded from script 20 and checked against its saved coefficient',
  'CSV and script 21B prediction-order validation; no model refit occurs.',
  'Tarjei GLLVM plot aesthetics retained (dashed zero, 95% horizontal intervals,',
  'classic theme, italic taxonomy; MiFish green summer/blue winter/grey nonsignificant.',
  'COI green summer/orange autumn/blue winter/grey nonsignificant.',
  'Year coefficients remain red (cormorant), blue (seal) if significant.',
  'Two independent legends: MiFish on panel A, COI on panel E.',
  'Original four-season model is NOT reproduced and must not be claimed as such.',
  'All taxa displayed, grey if Wald CI spans zero; original significance filter',
  'would remove additional taxa, but completeness is prioritized here.',
  '95% Wald CI are numerical Hessian intervals from original live fits,',
  'not independent evidence of CI coverage or causal effects.',
  'COI cormorant: 29A original-method reconstruction, screened reliable coefficients.',
  'COI archive does NOT record a confirmed convergence flag; do not call this model validated.',
  'COI seal failed; not included in manuscript coefficient panels.',
  'COI seasonal panel has Summer/Autumn/Winter vs Spring; MiFish has half-years.',
  'COI season and year CI drawn at fixed numerical y coordinates, never dodged.',
  'Both season legends are mapped directly to their significance display groups.',
  'MiFish validated coefficient table in output; COI screened table in validation.',
  'One combined six-panel manuscript figure; QA and provenance in validation.'),
  file.path(VAL,'06_METHOD_NOTES.txt'))
capture.output(sessionInfo(),file=file.path(VAL,'06_SESSION_INFO.txt'))
cat('\n05 COMPLETE\nFigure: ',fig_path,'\nTable: ',file.path(TAB,'gllvm_validated_coefficients.csv'),
    '\nValidation: ',VAL,'\n',sep='')
