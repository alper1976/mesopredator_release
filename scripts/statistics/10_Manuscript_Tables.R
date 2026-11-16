# PORTABLE STANDALONE RECONSTRUCTION. Original calculations retained.
# ============================================================================
# 10_Manuscript_Tables.R — Mesopredator article manuscript (NOT master's thesis)
# Source: Mesopredator_release(2).zip: main.tex + supplementary_material.tex
# Reuses validated outputs from modules 01, 03, 04 and optionally original
# MiFish RDS for SCRuB removal membership. Does NOT refit models.
# Two main tables and five supplementary tables; no manuscript values are
# substituted for independently computed results.
# ============================================================================

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
# Output folder is created automatically on each run.
OUT <- file.path(BASE,'processed_data','manuscript_tables')
MAIN <- file.path(OUT,'main'); SUPP <- file.path(OUT,'supplementary')
VAL <- file.path(BASE,'validation','10_Manuscript_Tables')
for(d in c(MAIN,SUPP,VAL)) dir.create(d,recursive=TRUE,showWarnings=FALSE)

# Inputs expected from accepted stages. The original source files stay read-only.
source_files <- c(
 threshold=file.path(BASE,'processed_data','02_Rarefaction_Threshold_Test','tables','threshold_summary.csv'),
 permanova=file.path(BASE,'processed_data','04_PERMANOVA_PCoA','PERMANOVA_results.csv'),
 mifish=file.path(BASE,'processed_data','05_Diet_Composition_FOO_RRA','tables','taxa_list_mifish.csv'),
 coi_invertebrates=file.path(BASE,'processed_data','05_Diet_Composition_FOO_RRA','tables','taxa_list_coi_invertebrates.csv'),
 coi_full=file.path(BASE,'processed_data','05_Diet_Composition_FOO_RRA','tables','taxa_list_coi_full.csv')
)
missing <- source_files[!file.exists(source_files)]
if(length(missing)) stop('Missing approved module outputs:\n',paste(missing,collapse='\n'),
                         '\nPlease run the relevant earlier approved modules; no fallback values are invented.')
read_source <- function(name) read.csv(source_files[[name]],check.names=FALSE,stringsAsFactors=FALSE)
mi <- read_source('mifish'); inv <- read_source('coi_invertebrates')
co <- read_source('coi_full'); pm <- read_source('permanova')
th <- read_source('threshold')
assert_columns <- function(dat,cols,file) {
  absent <- setdiff(cols,names(dat))
  if(length(absent)) stop(file,' missing columns: ',paste(absent,collapse=', '))
}
for(nm in c('mifish','coi_invertebrates','coi_full')) {
 dat <- get(c(mifish='mi',coi_invertebrates='inv',coi_full='co')[[nm]])
 assert_columns(dat,c('Taxon','LCA_rank','FOO_cormorant','FOO_seal'),nm)
}
assert_columns(inv,'Class','coi_invertebrates')
assert_columns(co,c('Class','prey_type'),'coi_full')
assert_columns(pm,c('marker','predator','model','n','term','df','R2','F','p'),'PERMANOVA_results')
assert_columns(th,c('Threshold','Samples','ASVs','Unique_LCA_labels','Identified_LCA_taxa','Species_level_ASVs'),'threshold_summary')

# Comparison helper: mismatches are saved and WARNED, never silently fixed.
checks <- list()
record_check <- function(table_id,measure,computed,manuscript,note='') {
 a <- suppressWarnings(as.numeric(computed)); b <- suppressWarnings(as.numeric(manuscript))
 same <- if(is.na(a)||is.na(b)) identical(as.character(computed),as.character(manuscript)) else
  abs(a-b) <= if(grepl('PERMANOVA',table_id)) 0.0011 else 1e-9
 checks[[length(checks)+1L]] <<- data.frame(table_id=table_id,measure=measure,
      computed=as.character(computed),manuscript=as.character(manuscript),
      matches=isTRUE(same),note=note,stringsAsFactors=FALSE)
}

# --------------------- MAIN TABLE: prey taxa summary -------------------------
# Method selection established by audits 10B–10E:
# MiFish => distinct detected ASV IDs (script 02 identified/global 0.1% object)
# COI => distinct detected FULL LCA strings, independently calculated by script
# 04 (counts from the accepted FOO table). Never divide or edit sequences to
# force the manuscript values. The main table mixes these units; document it.
if (!requireNamespace('phyloseq',quietly=TRUE))
 stop('Please install phyloseq: MiFish main summary needs the audited 02 ASV object')
mi_rds <- file.path(BASE,'processed_data','03_Phyloseq_Preprocessing',
                    'ps_mifish_global_0_1pct_identified.rds')
if (!file.exists(mi_rds)) stop('Missing approved MiFish ASV input: ',mi_rds)
ps <- readRDS(mi_rds)
if (!inherits(ps,'phyloseq')) stop('MiFish input is not a phyloseq object')
mat <- as(phyloseq::otu_table(ps),'matrix')
if (phyloseq::taxa_are_rows(ps)) mat <- t(mat)
md <- as.data.frame(phyloseq::sample_data(ps),check.names=FALSE)
tx <- as.data.frame(phyloseq::tax_table(ps),check.names=FALSE)
if (!setequal(rownames(mat),rownames(md)) ||
    !setequal(colnames(mat),rownames(tx))) stop('MiFish IDs are not aligned')
md <- md[rownames(mat),,drop=FALSE]
tx <- tx[colnames(mat),,drop=FALSE]
if (!'species' %in% names(md) || !'LCA_rank' %in% names(tx))
 stop('MiFish input missing species or LCA_rank')
sp <- trimws(as.character(md$species))
sp[sp %in% c('top_cormorant','Phalacrocorax carbo')] <- 'cormorant'
sp[sp %in% c('Phoca vitulina','harbour_seal')] <- 'harbor_seal'
c_ids <- rownames(md)[sp=='cormorant']
s_ids <- rownames(md)[sp=='harbor_seal']
if (!length(c_ids) || !length(s_ids)) stop('MiFish predator groups missing')
if (anyNA(mat) || any(!is.finite(mat)) || any(mat<0)) stop('Invalid ASV abundances')
seen_c <- colSums(mat[c_ids,,drop=FALSE] > 0) > 0
seen_s <- colSums(mat[s_ids,,drop=FALSE] > 0) > 0
rk <- tolower(trimws(as.character(tx$LCA_rank)))
asv_count <- function(mask=rep(TRUE,ncol(mat))) {
 a <- colnames(mat)[which(mask & seen_c)]
 b <- colnames(mat)[which(mask & seen_s)]
 c(cormorant=length(a),seal=length(b),shared=length(intersect(a,b)),
   total=length(union(a,b)))
}
counts_by <- function(d,predicate=rep(TRUE,nrow(d))) {
 d <- d[which(predicate),,drop=FALSE]
 a <- unique(as.character(d$Taxon[as.numeric(d$FOO_cormorant)>0]))
 b <- unique(as.character(d$Taxon[as.numeric(d$FOO_seal)>0]))
 c(cormorant=length(a),seal=length(b),shared=length(intersect(a,b)),
   total=length(union(a,b)))
}
make_row <- function(marker,group,v) {
 data.frame(marker=marker,group=group,cormorant=unname(v['cormorant']),
            harbour_seal=unname(v['seal']),shared=unname(v['shared']),
            total=unname(v['total']),stringsAsFactors=FALSE)
}
main_taxa <- do.call(rbind,list(
 make_row('MiFish','Actinopteri',asv_count()),
 make_row('MiFish','Species-level',asv_count(rk=='species')),
 make_row('MiFish','Genus-level',asv_count(rk=='genus')),
 make_row('COI Leray-XT','Actinopteri',counts_by(co,as.character(co$Class)=='Actinopteri')),
 make_row('COI Leray-XT','Polychaeta',counts_by(co,as.character(co$Class)=='Polychaeta')),
 make_row('COI Leray-XT','Malacostraca',counts_by(co,as.character(co$Class)=='Malacostraca')),
 make_row('COI Leray-XT','Hexanauplia',counts_by(co,as.character(co$Class)=='Hexanauplia')),
 make_row('COI Leray-XT','Branchiopoda',counts_by(co,as.character(co$Class)=='Branchiopoda')),
 make_row('COI Leray-XT','Total invertebrate',
          counts_by(co,!is.na(co$Class) & as.character(co$Class)!='Actinopteri'))
))
main_taxa$counting_unit <- ifelse(main_taxa$marker=='MiFish','ASV_ID','FULL_LCA_LABEL')
write.csv(main_taxa,file.path(MAIN,'Table_prey_taxa_summary.csv'),row.names=FALSE)
ref_taxa <- data.frame(
 marker=c(rep('MiFish',3),rep('COI Leray-XT',6)),
 group=c('Actinopteri','Species-level','Genus-level','Actinopteri',
         'Polychaeta','Malacostraca','Hexanauplia','Branchiopoda','Total invertebrate'),
 cormorant=c(45,43,2,38,6,7,8,1,22),
 harbour_seal=c(42,41,1,34,5,5,3,1,16),
 shared=c(18,18,0,19,4,4,3,1,12),
 total=c(69,67,2,53,6,10,8,2,26))
for(i in seq_len(nrow(main_taxa))) for(v in c('cormorant','harbour_seal','shared','total'))
 record_check('MAIN_prey_taxa',paste(main_taxa$marker[i],main_taxa$group[i],v,sep=' / '),
              main_taxa[[v]][i],ref_taxa[[v]][i])
ref_taxa$expected_union <- with(ref_taxa,cormorant+harbour_seal-shared)
ref_taxa$union_gap <- ref_taxa$total-ref_taxa$expected_union
ref_taxa$internally_consistent <- ref_taxa$union_gap==0
write.csv(ref_taxa,file.path(VAL,'10_MANUSCRIPT_MAIN_REFERENCE_AND_UNION_AUDIT.csv'),row.names=FALSE)
main_audit <- merge(main_taxa,ref_taxa,by=c('marker','group'),suffixes=c('_computed','_manuscript'),sort=FALSE)
write.csv(main_audit,file.path(VAL,'10_MAIN_TAXA_METHOD_AND_DEVIATIONS.csv'),row.names=FALSE)
asv_lineage <- data.frame(ASV_ID=colnames(mat),LCA_rank=rk,
 LCA_taxon=if('LCA_taxon' %in% names(tx)) as.character(tx$LCA_taxon) else NA_character_,
 present_cormorant=unname(seen_c),present_seal=unname(seen_s))
write.csv(asv_lineage,file.path(VAL,'10_MIFISH_ASV_LINEAGE.csv'),row.names=FALSE)

# --------------------- MAIN TABLE: additive PERMANOVA ------------------------
pa <- pm[pm$model=='additive' & pm$term %in% c('season','sample_year','Residual'),,drop=FALSE]
pa <- pa[order(match(pa$marker,c('MiFish','COI')),
               match(pa$predator,c('cormorant','harbor_seal')),
               match(pa$term,c('season','sample_year','Residual'))),,drop=FALSE]
if(nrow(pa)!=12L) stop('Expected 12 rows in additive PERMANOVA, got ',nrow(pa))
pa$Species <- ifelse(pa$predator=='cormorant','Phalacrocorax carbo','Phoca vitulina')
pa$Term <- ifelse(pa$term=='sample_year','Year',ifelse(pa$term=='season','Season','Residual'))
pa$Marker <- ifelse(pa$marker=='MiFish','MiFish','COI Leray-XT')
pa_out <- pa[,c('Marker','Species','n','Term','df','R2','F','p')]
write.csv(pa_out,file.path(MAIN,'Table_PERMANOVA_combined.csv'),row.names=FALSE,na='')
# Source: article main.tex table label tab:permanova_combined. Compare at
# manuscript's printed 3-decimal precision (no numerical values substituted).
perma_ref <- data.frame(
 marker=rep(c('MiFish','COI'),each=6),
 predator=rep(rep(c('cormorant','harbor_seal'),each=3),2),
 term=rep(c('season','sample_year','Residual'),4),
 n=rep(c(68,59,64,60),each=3),
 df=c(3,1,63,3,1,54,3,1,59,3,1,55),
 R2=c(.047,.010,.943,.133,.000,.866,.129,.001,.871,.101,.032,.867),
 F=c(1.057,.659,NA,2.772,.029,NA,2.902,.052,NA,2.138,2.030,NA),
 p=c(.387,.513,NA,.046,.950,NA,.017,.982,NA,.059,.130,NA),
 stringsAsFactors=FALSE)
for(i in seq_len(nrow(pa))) {
 tag <- paste(pa$marker[i],pa$predator[i],pa$term[i],sep=' / ')
 for(v in c('n','df','R2','F','p')) {
  obs <- pa[[v]][i]; ref <- perma_ref[[v]][i]
  if(is.na(obs)&&is.na(ref)) next
  record_check('MAIN_PERMANOVA',paste(tag,v),obs,ref,
    if(v %in% c('R2','F','p')) 'Compared to published rounded precision (3 decimals)' else '')
 }
}

# ------------------- SUPPLEMENT: threshold retention -------------------------
threshold_ref <- data.frame(Threshold=c('None','0.1%','0.5%','1%','2%','5%'),
 Samples=c(161,161,161,161,161,161),ASVs=c(167,160,145,136,130,114),
 Unique_LCA_labels=c(50,50,47,46,43,40),
 Species_level_ASVs=c(67,67,62,58,54,45))
if(nrow(th)!=6L || !identical(as.character(th$Threshold),threshold_ref$Threshold))
 stop('Threshold rows/order differ; inspect 01 threshold_summary.csv before using')
write.csv(th,file.path(SUPP,'Table_threshold_summary.csv'),row.names=FALSE)
for(v in setdiff(names(threshold_ref),'Threshold')) for(i in seq_len(nrow(th)))
 record_check('SUPP_threshold',paste(th$Threshold[i],v),th[[v]][i],threshold_ref[[v]][i],
  if(v=='Unique_LCA_labels') 'Includes label no identification, which is not a biological taxon' else '')

# -------------------- SUPPLEMENT: taxa lists ---------------------------------
# Full, computed taxon names are preserved: no cosmetic species-name splitting
# or manual row elimination to mimic the manuscript's editorial taxon labels.
write.csv(mi,file.path(SUPP,'Table_MiFish_fish_taxa.csv'),row.names=FALSE)
write.csv(inv,file.path(SUPP,'Table_COI_invertebrate_taxa.csv'),row.names=FALSE)
coi_fish <- co[!is.na(co$Class) & as.character(co$Class)=='Actinopteri',
               c('Taxon','LCA_rank','FOO_cormorant','FOO_seal'),drop=FALSE]
coi_fish <- coi_fish[order(-pmax(coi_fish$FOO_cormorant,coi_fish$FOO_seal),coi_fish$Taxon),]
write.csv(coi_fish,file.path(SUPP,'Table_COI_fish_taxa_FULL.csv'),row.names=FALSE)
# These are coverage indicators, NOT equivalence checks against rendered
# supplement: original COI fish table manually merged labels and omitted six.
record_check('SUPP_MiFish_taxa','distinct LCA table rows',nrow(mi),47,
             '47 supplementary full-LCA rows; compare individual taxon names where required')
record_check('SUPP_COI_invertebrates','distinct taxa',nrow(inv),26,
             'Compare row-wise taxon labels and rounded FOO before claiming exact match')
record_check('SUPP_COI_fish','full distinct fish taxa',nrow(coi_fish),53,
             'Manuscript displays 47 COI fish rows; six extra full-LCA rows are documented, exclusion rule not confirmed')

# ------------------ SUPPLEMENT: removed SCRuB samples ------------------------
# A manuscript list is not equivalent to a computationally verified removal.
# Recover removal candidates from archived pre/post-SCRuB raw phyloseq objects
# only if both raw objects and phyloseq are available. Otherwise record PENDING.
scrub_reference <- data.frame(
 library=c('d-2-30','d-2-37','d-3-16','d-3-19','d-3-20','d-4-18','d-4-33'),
 original_sample=c('A2015/05068','KS_2023_13_B','A2020/05619','HL0320_05','HL0320_04','A2021/05114','OS41_D'),
 species=c('Harbour seal','Harbour seal','Harbour seal','Cormorant','Cormorant','Harbour seal','Harbour seal'),
 origin=c('SLU','SLU','SLU','SLU','SLU','SLU','Oslofjord'),
 extraction_session=c('10','11','11','12','12','15','Pilot'),
 stringsAsFactors=FALSE)
raw_dirs <- file.path(BASE, 'rawfiles')
locate_raw <- function(name) {
 p <- file.path(raw_dirs,name); p <- p[file.exists(p)]
 if(length(p)) normalizePath(p[[1]]) else NULL
}
prepath <- locate_raw('ps_filtered_mifish.rds')
postpath <- locate_raw('ps_decon_mifish.rds')
scrub_state <- 'NOT_REPRODUCED: insufficient independent removal-membership evidence'
scrub_candidate_path <- file.path(VAL,'10_SCRUB_ZERO_READ_CANDIDATES.csv')
if(!is.null(prepath) && !is.null(postpath) && requireNamespace('phyloseq',quietly=TRUE)) {
 tryCatch({
  pre <- readRDS(prepath); post <- readRDS(postpath)
  pre_ids <- phyloseq::sample_names(pre); post_ids <- phyloseq::sample_names(post)
  pre_counts <- phyloseq::sample_sums(pre); post_counts <- phyloseq::sample_sums(post)
  post_common <- names(pre_counts)[names(pre_counts) %in% names(post_counts)]
  # Candidate rule is explicitly provisional: zero after SCRuB, positive before.
  candidate <- post_common[pre_counts[post_common]>0 & post_counts[post_common]==0]
  write.csv(data.frame(library=candidate,pre_reads=as.numeric(pre_counts[candidate]),
       post_reads=as.numeric(post_counts[candidate])),
       scrub_candidate_path,row.names=FALSE)
  scrub_audit <- data.frame(library=union(candidate,scrub_reference$library),
     independently_zeroed=union(candidate,scrub_reference$library) %in% candidate,
     listed_in_manuscript=union(candidate,scrub_reference$library) %in% scrub_reference$library)
  write.csv(scrub_audit,file.path(VAL,'10_SCRUB_MEMBERSHIP_AUDIT.csv'),row.names=FALSE)
  for(s in scrub_reference$library) record_check('SUPP_SCRUB_removed',s,
        s %in% candidate,TRUE,'Zero-read rule is a candidate, not independently verified SCRuB removal status')
  scrub_state <- if(setequal(candidate,scrub_reference$library))
     'CANDIDATES_MATCH: zero-read source comparison; verify SCRuB status and session metadata' else
     'CANDIDATES_DIFFER: inspect independent zero-read source comparison'
 },error=function(e) {
  scrub_state <<- paste('NOT_REPRODUCED: raw-object comparison failed:',conditionMessage(e))
 })
}
# Reference is a benchmark ONLY (in validation). Not emitted as a reproduced
# main/supplement table unless independently computed membership is available.
write.csv(scrub_reference,file.path(VAL,'10_SCRUB_MANUSCRIPT_REFERENCE_ONLY.csv'),row.names=FALSE)
# Export the seventh table explicitly. It is a transcription, not a verified
# reconstruction, and the independent status is separately documented.
write.csv(scrub_reference,file.path(SUPP,'Table_SCRUB_removed_MANUSCRIPT_REFERENCE.csv'),row.names=FALSE)

# A positive source-object comparison can generate independent membership
# candidates, but not reproduce the complete manuscript sample/session table.
writeLines(scrub_state,file.path(SUPP,'Table_SCRUB_removed_STATUS.txt'))

# ------------------------------- QA ------------------------------------------
check_df <- if(length(checks)) do.call(rbind,checks) else data.frame()
write.csv(check_df,file.path(VAL,'10_MANUSCRIPT_VALUE_COMPARISON.csv'),row.names=FALSE)
write.csv(data.frame(source=names(source_files),path=unname(source_files),
                     md5=unname(tools::md5sum(source_files))),
          file.path(VAL,'10_SOURCE_PROVENANCE.csv'),row.names=FALSE)
manifest <- data.frame(
 section=c('MAIN','MAIN',rep('SUPPLEMENTARY',5)),
 table=c('Prey taxa summary','Additive PERMANOVA','Removed SCRuB samples',
         'Global relative abundance threshold','MiFish taxa list',
         'COI invertebrate taxa list','COI fish taxa list'),
 output=c(file.path(MAIN,'Table_prey_taxa_summary.csv'),
          file.path(MAIN,'Table_PERMANOVA_combined.csv'),
          file.path(SUPP,'Table_SCRUB_removed_MANUSCRIPT_REFERENCE.csv'),
          file.path(SUPP,'Table_threshold_summary.csv'),
          file.path(SUPP,'Table_MiFish_fish_taxa.csv'),
          file.path(SUPP,'Table_COI_invertebrate_taxa.csv'),
          file.path(SUPP,'Table_COI_fish_taxa_FULL.csv')),
 stringsAsFactors=FALSE)
manifest$state <- ifelse(manifest$table!='Removed SCRuB samples',
                        'RECONSTRUCTED_WITH_AUDIT','MANUSCRIPT_REFERENCE_ONLY')
write.csv(manifest,file.path(VAL,'10_TABLE_MANIFEST.csv'),row.names=FALSE)
writeLines(c(
 'Scope: ARTICLE MANUSCRIPT Mesopredator_release(2).zip (main.tex and supplementary_material.tex).',
 'NOT Tarjei Røstum master thesis tables 2.1-2.13, 3.1-3.5 or A.1-A.3.',
 'There are two tables in article main.tex and five tables in supplementary_material.tex.',
 'All candidate numbers come from prior reconstructed outputs, never from the manuscript benchmark.',
 'Manuscript numbers are used for checking only, saved to ../../docs/validation/10_MANUSCRIPT_VALUE_COMPARISON.csv.',
 'Main prey taxa: MiFish counts distinct ASV IDs from approved 02 0.1%-identified object.',
 'Main prey taxa: COI counts distinct full-LCA labels from approved 04 FOO export.',
 'These are different counting units and should not be described as directly equivalent.',
 'Manuscript MiFish 45/42/18/69 matches ASVs, ranks differ by one skarv ASV.',
 'Manuscript COI Actinopteri 38/34/19/53 and invertebrate total 22/16/12/26 match LCA.',
 'Several manuscript class-level shared / predator counts differ or fail union arithmetic.',
 'Manuscript internal gaps are in 10_MANUSCRIPT_MAIN_REFERENCE_AND_UNION_AUDIT.csv.',
 'Full COI fish taxa are retained (53 LCA rows), manuscript supplement displays 47;',
 ' six differences are documented without inferring an unverified removal rule.',
 'Supplementary FOO row labels and taxon matching still require manual audit.',
 paste('SCRuB:',scrub_state),
 'SCRuB table in supplementary is explicitly MANUSCRIPT_REFERENCE, not independent reproduction;',
 ' zero-read pre/post comparison is a candidate and not full decontamination proof.',
 'PERMANOVA comparisons tolerate 0.0011 for manuscript figures rounded to 3 decimals;',
 ' statistically substantive differences require re-checking analysis and sample cohorts.',
 'No modelling or refitting is performed by script 10.'
 ),file.path(VAL,'10_METHOD_AND_DEVIATION_NOTE.txt'))
capture.output(sessionInfo(),file=file.path(VAL,'10_SESSION_INFO.txt'))
cat('\n10 COMPLETE — ARTICLE MANUSCRIPT TABLE CANDIDATES\n')
cat('Main tables: ',MAIN,'\nSupplement: ',SUPP,'\nValidation: ',VAL,'\n',sep='')
cat('Comparison mismatches: ',sum(!check_df$matches),' of ',nrow(check_df),'\n',sep='')
cat('SCRuB: ',scrub_state,'\n',sep='')
if(any(!check_df$matches)) warning('Manuscript table comparisons have differences; inspect 10_MANUSCRIPT_VALUE_COMPARISON.csv. No values were altered.')

# Additional concise report documenting fixed vs unresolved findings from 10A–10E.
writeLines(c(
 'ARTICLE MANUSCRIPT TABLE RECONSTRUCTION – FINAL METHOD / DEVIATION NOTE',
 '',
 'Scope: Mesopredator manuscript, main.tex + supplementary_material.tex; NOT master thesis.',
 'All derived values come from accepted script 01/02/03/04 exports; no model is refitted.',
 '',
 'MiFish main-table counts: unique detected ASV IDs from script 02, 0.1% identified object.',
 'The manuscript all-MiFish counts 45 cormorant, 42 seal, 18 shared, 69 total match.',
 'The rank split matches totals 67 species-rank ASVs and 2 genus-rank ASVs, but',
 'cormorant split here is 44 species and 1 genus versus manuscript 43 and 2.',
 'MiFish supplementary table: 47 LCA-labelled rows, same FOO-pair distribution',
 'as the manuscript; ASV counts should not be presented as biological species counts.',
 '',
 'COI main-table counts: unique full-LCA strings from accepted script 04.',
 'Actinopteri: 38 cormorant, 34 seal, 19 shared, 53 total – exact main-table match.',
 'All invertebrate: 22 cormorant, 16 seal, 12 shared, 26 total – exact match.',
 'Selected COI class-level predator/shared cells differ from the manuscript.',
 'Within manuscript, some class rows violate total = cormorant + seal - shared.',
 'See 10_MANUSCRIPT_MAIN_REFERENCE_AND_UNION_AUDIT.csv for the exact rows.',
 '',
 'COI fish supplementary: full computed export has 53 rows, published table 47.',
 'Six rows not displayed in manuscript; no deletion is imposed on computed data.',
 'COI invertebrate supplementary: 26 rows and matching distribution of FOO pairs.',
 'All supplement-level FOO matching is pair-distribution evidence, not necessarily',
 'taxon-label equivalence; compound LCA labels can prevent one-to-one name matching.',
 '',
 'SCRuB removed samples: supplied table is marked MANUSCRIPT_REFERENCE_ONLY.',
 'An independent zero-read comparison does NOT establish original SCRuB-removal status.',
 '',
 'Comparison/QA outputs retained separately under ../../docs/validation/10_Manuscript_Tables.',
 'No source data, approved upstream scripts, or original manuscript are altered.'
),file.path(VAL,'10_FINAL_METHOD_AND_DEVIATIONS.txt'))
