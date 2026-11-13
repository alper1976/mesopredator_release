# PORTABLE STANDALONE RECONSTRUCTION. Original calculations retained.
# ============================================================================
# 02_PHYLOSEQ_PREPROCESSING.R
# Mesopredator | final reconstruction | candidate v1.0
#
# METHOD SOURCES (Tarjei Røstum originals):
#   setup.n.scrub.2026.04.20.Rmd
#   pruned100.general.analysis.Rmd
#   percentage.test.Rmd
#   prevalence.distribution.coi.mifish.Rmd
#
# CROSS-CHECKED RECONSTRUCTIONS:
#   04_reconstruct_mifish_objects.R
#   07_MiFish_global_abundance_filter.R
#   23_COI_original_reconstruction_audit.R
#
# PURPOSE: generate the ORIGINAL-derived phyloseq inputs for later modules.
# No figures, exploratory analyses, or manuscript tables are generated here.
# Derived RDS files are replaceable; original input RDS files are READ ONLY.
# Run from project root OR from 00_final_reconstruction_scripts.
# ============================================================================

if (!requireNamespace("phyloseq", quietly = TRUE)) {
  stop("Install the phyloseq package before running script 02.")
}
suppressPackageStartupMessages(library(phyloseq))

# ---- 0. Project paths --------------------------------------------------------
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
OUTPUT_DIR <- file.path(BASE, "processed_data", "03_Phyloseq_Preprocessing")
VALID_DIR <- file.path(BASE, "validation", "03_Phyloseq_Preprocessing")
for (d in c(OUTPUT_DIR, VALID_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Exact-basename search, restricted to known ORIGINAL data folders. No search in
# old reconstruction directories: we must never feed earlier derived RDS objects
# into the new final pipeline without explicitly saying so.
find_original <- function(filename, relative_dirs) {
  paths <- file.path(ROOT, relative_dirs, filename)
  paths <- paths[file.exists(paths)]
  if (!length(paths)) stop("Missing original ", filename, ". Checked: ",
                           paste(file.path(ROOT, relative_dirs, filename), collapse = " | "))
  if (length(paths) > 1L && length(unique(unname(tools::md5sum(paths)))) != 1L) {
    stop("Conflicting original copies of ", filename, ": ", paste(paths, collapse = " | "))
  }
  normalizePath(paths[[1]], mustWork = TRUE)
}

MIFISH_DIRS <- "rawfiles"
COI_DIRS <- "rawfiles"
source_files <- c(
  ps_decon_mifish = find_original("ps_decon_mifish.rds", MIFISH_DIRS),
  ps_filtered_mifish = find_original("ps_filtered_mifish.rds", MIFISH_DIRS),
  ps_coi_metazoa = find_original("ps_coi_metazoa.rds", COI_DIRS)
)

ps_decon_original <- readRDS(source_files[["ps_decon_mifish"]])
ps_filtered_original <- readRDS(source_files[["ps_filtered_mifish"]])
ps_coi_original <- readRDS(source_files[["ps_coi_metazoa"]])
stopifnot(inherits(ps_decon_original, "phyloseq"),
          inherits(ps_filtered_original, "phyloseq"),
          inherits(ps_coi_original, "phyloseq"))

OS19 <- c("d-4-29", "d-4-30", "d-4-31", "d-4-32")
if (!all(OS19 %in% sample_names(ps_decon_original))) stop("Four OS19 IDs absent from original MiFish")

# ---- 1. Explicit metadata repairs from original scripts ---------------------
fix_metadata <- function(ps, fix_ocean = FALSE) {
  md <- as.data.frame(sample_data(ps), check.names = FALSE)
  required <- c("species", "sample_month", "sample_year", "Sample_or_Control",
                "person_emptied")
  if (!all(required %in% names(md))) {
    stop("Missing metadata columns: ", paste(setdiff(required, names(md)), collapse = ", "))
  }
  md$species <- gsub("top_cormorant", "cormorant", as.character(md$species), fixed = TRUE)
  if ("origin" %in% names(md)) {
    md$origin <- gsub("oslofjorden", "oslofjord", as.character(md$origin), fixed = TRUE)
  }
  md$sample_month <- tolower(trimws(as.character(md$sample_month)))
  md$sample_month <- gsub("augsut", "august", md$sample_month, fixed = TRUE)
  md$sample_month <- gsub("oktober", "october", md$sample_month, fixed = TRUE)
  md$processor <- ifelse(as.character(md$person_emptied) == "SLU", "SLU", "Self_processed")
  md$season <- ifelse(md$sample_month %in% c("march", "april", "may"), "Spring",
               ifelse(md$sample_month %in% c("june", "july", "august"), "Summer",
               ifelse(md$sample_month %in% c("september", "october", "november"), "Autumn",
               ifelse(md$sample_month %in% c("december", "january", "february"), "Winter", "Unknown"))))
  if (fix_ocean && "ocean" %in% names(md)) {
    md$ocean <- gsub("oslofjorden", "oslofjord", as.character(md$ocean), fixed = TRUE)
    md$ocean <- gsub("hvaler", "oslofjord", md$ocean, fixed = TRUE)
  }
  rownames(md) <- sample_names(ps)
  sample_data(ps) <- sample_data(md)
  ps
}

keep_biological <- function(ps) {
  md <- as.data.frame(sample_data(ps), check.names = FALSE)
  ids <- rownames(md)[!is.na(md$Sample_or_Control) &
                        as.character(md$Sample_or_Control) == "True_Sample"]
  prune_samples(ids, ps)
}
remove_os19 <- function(ps) prune_samples(!(sample_names(ps) %in% OS19), ps)
prune_empty_taxa <- function(ps) prune_taxa(taxa_sums(ps) > 0, ps)

# Return counts as sample x ASV for independent abundance checks.
count_matrix <- function(ps) {
  z <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) z <- t(z)
  storage.mode(z) <- "double"
  z
}

# Record only provenance/QA, not a new manuscript table.
records <- list()
record <- function(name, object, rule) {
  mat <- count_matrix(object)
  records[[length(records) + 1L]] <<- data.frame(
    object = name, samples = nsamples(object), ASVs = ntaxa(object),
    abundance_sum = sum(mat), zero_samples = sum(rowSums(mat) == 0),
    transformation = rule, stringsAsFactors = FALSE
  )
  cat(sprintf("%-47s n=%3d ASVs=%3d zero=%2d\n", name, nsamples(object),
              ntaxa(object), sum(rowSums(mat) == 0)))
}
save_derived <- function(ps, name, rule) {
  path <- file.path(OUTPUT_DIR, paste0(name, ".rds"))
  saveRDS(ps, path)
  record(name, ps, rule)
  invisible(path)
}

# ---- 2. MiFish: original read-depth objects -------------------------------
# The original setup generated 100-/500-read objects *before* later biological
# sample filtering. Retain those as separate intermediate objects.
ps_pruned_100 <- prune_empty_taxa(prune_samples(sample_sums(ps_decon_original) >= 100,
                                               ps_decon_original))
ps_pruned_500 <- prune_empty_taxa(prune_samples(sample_sums(ps_decon_original) >= 500,
                                               ps_decon_original))
save_derived(ps_pruned_100, "ps_pruned_mifish_100reads", "Original setup: >=100 reads")
save_derived(ps_pruned_500, "ps_pruned_mifish_500reads", "Original setup: >=500 reads")

# In pruned100.general.analysis.Rmd, ps_pruned metadata was replaced with the
# corresponding records from ps_filtered before harmonisation and rarefaction.
md_filtered <- as.data.frame(sample_data(ps_filtered_original), check.names = FALSE)
if (!all(sample_names(ps_pruned_100) %in% rownames(md_filtered))) {
  stop("Original ps_filtered metadata lacks one or more 100-read samples")
}
ps_pruned_analysis <- ps_pruned_100
sample_data(ps_pruned_analysis) <- sample_data(md_filtered[sample_names(ps_pruned_analysis), , drop = FALSE])
ps_pruned_analysis <- fix_metadata(ps_pruned_analysis)
ps_pruned_analysis <- prune_empty_taxa(keep_biological(ps_pruned_analysis))
if (!all(sample_sums(ps_pruned_analysis) >= 100)) stop("100-read analysis object contains low-depth samples")
save_derived(ps_pruned_analysis, "ps_pruned_mifish_100reads_analysis",
             "Metadata restored from ps_filtered; true biological samples")

# Exact rarefaction parameters from Tarjei's pruned100.general.analysis.Rmd.
set.seed(123)
ps_rarefied <- rarefy_even_depth(ps_pruned_analysis, sample.size = 100,
                                  rngseed = 123, replace = FALSE, verbose = FALSE)
if (!all(sample_sums(ps_rarefied) == 100)) stop("Rarefaction did not yield 100 reads/sample")
save_derived(ps_rarefied, "ps_rarefied_mifish_100reads",
             "rarefy_even_depth(sample.size=100,rngseed=123,replace=FALSE)")

# ---- 3. MiFish: original within-sample 0.1% and 1% objects ----------------
# IMPORTANT: percentage.test.Rmd zeroes individual ASV values below a sample
# threshold, then drops entirely empty taxa/samples. This is NOT the same as
# the GLOBAL maximum-ASV filter used for the 161 x 160 thesis threshold table.
# This route includes four OS19 libraries because certain later original Rmd
# blocks explicitly reload ps_1.rds with them present (see GG2 workflow).
ps_all_bio <- fix_metadata(keep_biological(ps_decon_original))
if (any(sample_sums(ps_all_bio) <= 0)) stop("Zero-depth sample in original biological MiFish")
ps_all_pct <- transform_sample_counts(ps_all_bio, function(x) x / sum(x) * 100)
within_sample_threshold <- function(ps, cutoff) {
  z <- transform_sample_counts(ps, function(x) ifelse(x < cutoff, 0, x))
  prune_empty_taxa(prune_samples(sample_sums(z) > 0, z))
}
ps_01_within <- within_sample_threshold(ps_all_pct, 0.1)
ps_1_within <- within_sample_threshold(ps_all_pct, 1)
save_derived(ps_01_within, "ps_mifish_within_sample_0_1pct",
             "Zero <0.1% in each sample, then prune empty samples/taxa; OS19 included")
save_derived(ps_1_within, "ps_mifish_within_sample_1pct",
             "Zero <1% in each sample, then prune empty samples/taxa; OS19 included")

# ---- 4. MiFish: GLOBAL 0.1% taxon retention for threshold/PERMANOVA/FOO ----
# Remove OS19, keep true samples, and normalize BEFORE global selection.
ps_mf_global_base <- prune_empty_taxa(keep_biological(remove_os19(ps_decon_original)))
if (nsamples(ps_mf_global_base) != 161L || ntaxa(ps_mf_global_base) != 167L ||
    sum(sample_sums(ps_mf_global_base)) != 580002) {
  stop("MiFish 161 / 167 / 580002 baseline differs from independently validated script 01")
}
ps_mf_global_base <- fix_metadata(ps_mf_global_base, fix_ocean = TRUE)
ps_mf_percent <- transform_sample_counts(ps_mf_global_base, function(x) 100 * x / sum(x))
# Strict >0.1 cutoff matches six-row global threshold table confirmed in 01.
ps_mf_global <- filter_taxa(ps_mf_percent, function(x) max(x) > 0.1, prune = TRUE)
if (nsamples(ps_mf_global) != 161L || ntaxa(ps_mf_global) != 160L) {
  stop("MiFish global 0.1% count differs from 161 samples / 160 ASVs")
}
# Preserve within-sample percentage values of retained ASVs, without re-scaling.
check <- count_matrix(ps_mf_percent)[sample_names(ps_mf_global), taxa_names(ps_mf_global), drop = FALSE]
if (!isTRUE(all.equal(unname(check), unname(count_matrix(ps_mf_global)), tolerance = 1e-10))) {
  stop("Global MiFish filtering unexpectedly changed retained abundances")
}
save_derived(ps_mf_global, "ps_mifish_global_0_1pct",
             "ASV maximum across samples >0.1%; no within-sample zeroing")

# Some downstream figures explicitly exclude LCA_taxon == 'no identification'.
# Make the resulting object available without silently dropping empty samples.
mi_tax <- as.data.frame(tax_table(ps_mf_global), check.names = FALSE)
if (!"LCA_taxon" %in% names(mi_tax)) stop("MiFish taxonomy missing LCA_taxon")
ident <- !is.na(mi_tax$LCA_taxon) &
  tolower(trimws(as.character(mi_tax$LCA_taxon))) != "no identification"
if (sum(!ident) != 91L) stop("Expected exactly 91 unidentified MiFish ASVs")
ps_mf_identified <- prune_taxa(rownames(mi_tax)[ident], ps_mf_global)
if (ntaxa(ps_mf_identified) != 69L ||
    sum(sample_sums(ps_mf_identified) == 0) != 12L) {
  stop("MiFish identified-prey object differs from 69 taxa / 12 zero-prey samples")
}
save_derived(ps_mf_identified, "ps_mifish_global_0_1pct_identified",
             "Global 0.1%; remove no identification; retain 12 zero-prey samples")

# ---- 5. COI: original metazoa-to-prey preprocessing -----------------------
ps_coi <- fix_metadata(remove_os19(ps_coi_original), fix_ocean = TRUE)
keep_classes <- c("Actinopteri", "Malacostraca", "Mollusca", "Polychaeta",
                  "Hexanauplia", "Branchiopoda", "Echinodermata")
coi_tax <- as.data.frame(tax_table(ps_coi), check.names = FALSE)
if (!all(c("Class", "LCA_taxon") %in% names(coi_tax))) stop("COI taxonomy lacks Class/LCA_taxon")
ps_coi_prey <- subset_taxa(ps_coi, Class %in% keep_classes)
ps_coi_prey <- prune_empty_taxa(ps_coi_prey)
ps_coi_prey <- prune_samples(sample_sums(ps_coi_prey) > 0, ps_coi_prey)
if (any(sample_sums(ps_coi_prey) <= 0)) stop("COI prey object contains zero-depth samples")
save_derived(ps_coi_prey, "ps_coi_prey_raw_counts",
             "Original metazoa classes; OS19 excluded; zero-total samples removed")
ps_coi_pct <- transform_sample_counts(ps_coi_prey, function(x) 100 * x / sum(x))
# Unlike the verified MiFish table's >0.1, the original COI code uses >=0.1.
ps_coi_01 <- filter_taxa(ps_coi_pct, function(x) max(x) >= 0.1, prune = TRUE)
coi_tax <- as.data.frame(tax_table(ps_coi_01), check.names = FALSE)
ident_coi <- !is.na(coi_tax$LCA_taxon) &
  tolower(trimws(as.character(coi_tax$LCA_taxon))) != "no identification"
ps_coi_01 <- prune_taxa(rownames(coi_tax)[ident_coi], ps_coi_01)
# Do not delete zero-prey samples after last taxon filter (source behavior).
if (nsamples(ps_coi_01) != 145L || ntaxa(ps_coi_01) != 172L) {
  stop("COI object differs from previously verified 145 samples / 172 ASVs: observed ",
       nsamples(ps_coi_01), " / ", ntaxa(ps_coi_01),
       ". Audit input/source operations before adjusting filters.")
}
md_coi <- as.data.frame(sample_data(ps_coi_01), check.names = FALSE)
stopifnot(sum(as.character(md_coi$species) == "cormorant") == 80L,
          sum(as.character(md_coi$species) == "harbor_seal") == 65L)
save_derived(ps_coi_01, "ps_coi_global_0_1pct_identified",
             "COI max >=0.1%, LCA identified; original 80/65 biological samples")

# ---- 6. Validation only ---------------------------------------------------
summary <- do.call(rbind, records)
write.csv(summary, file.path(VALID_DIR, "03_OBJECT_SUMMARY.csv"), row.names = FALSE)
write.csv(data.frame(source = names(source_files), path = unname(source_files),
                     md5 = unname(tools::md5sum(source_files)), row.names = NULL),
          file.path(VALID_DIR, "03_ORIGINAL_INPUT_HASHES.csv"), row.names = FALSE)
write.csv(data.frame(sample_id = sample_names(ps_mf_identified)[sample_sums(ps_mf_identified) == 0],
                     reason = "No identified MiFish ASV after 0.1% global filter"),
          file.path(VALID_DIR, "03_MIFISH_EMPTY_IDENTIFIED_SAMPLES.csv"), row.names = FALSE)
writeLines(c(
  "03_Phyloseq_Preprocessing.R — METHOD NOTES",
  "Only original source files are inputs; no files in 03_reconstruction or 04_reproduced_analyses are read.",
  "Original raw RDS files are never overwritten.",
  "Two different MiFish thresholds exist in Tarjei's code:",
  "  A. Within-sample threshold: zero each ASV <0.1% / <1% in each sample; includes OS19 in the reloaded ps_1 route.",
  "  B. Global ASV threshold: retain entire ASVs if maximum across samples >0.1%; excludes OS19 and uses 161 biological samples.",
  "For the globally filtered MiFish object, no within-sample re-normalization is performed.",
  "For MiFish identified object, zero-prey samples are retained here and excluded only by analyses explicitly requiring nonzero rows.",
  "COI original >=0.1% threshold differs from MiFish global >0.1%; source code is preserved.",
  "COI 80/65 samples are not the PERMANOVA subgroup, and must not be conflated with PCoA n.",
  "Downstream module 03 should select the intended named object and apply original analysis-specific sample filters.",
  "No manuscript figures or tables are made in this preprocessing module."
), file.path(VALID_DIR, "03_METHOD_NOTES.txt"))
capture.output(sessionInfo(), file = file.path(VALID_DIR, "03_SESSION_INFO.txt"))
cat("\n02 COMPLETE: derived phyloseq objects in\n", OUTPUT_DIR,
    "\nValidation notes in\n", VALID_DIR, "\n")
