# ==============================================================================
# MESOPREDATOR MASTER v1.4 | 2012-2025 | ALL CORMORANT LIBRARIES
# Based on v1.3; 2025 extends diets, fish stocks and predator demography.
# Outputs are isolated under "final result v1.4 2012-2025".
# Cormorant census years remain provisional; 2025 is carried forward
# from the 2024 anchor using stats::approx(rule=2).
# ==============================================================================

# ==============================================================================
# MESOPREDATOR RELEASE — SKAGERRAK | UNIFIED ANNUAL ANALYSIS v1.0
# ==============================================================================
# Full pipeline: library QC -> RRA -> 18 diet models -> ICES and wrasse
# stock proxies -> predator demography -> 5000-draw Monte Carlo -> figures.
# NO SEASON TERM. Age-1 and age-2 cod biomass reconstructed from DATRAS.
#
# Study years: 2012–2025. Output relative to mesopredator_check.
# Important: 142.7 g is a pooled sample estimate, not a population mean
# independently validated for all years. ICES recruitment in thousands
# must be checked against the downloaded assessment.
#
# Cod_Juveniles stock denominator is revised below to AGE 1+2, not ages
# 0–2; the predator juvenile fraction was originally defined
# more broadly. Ratios are therefore scenario PRESSURE INDICES, NOT
# mortality or directly age-matched estimates.
#
# Other assumptions preserved for comparability: predator census,
# annual energetic demand, geographic scaling, diet RRA -> consumption.
# The ICES historic cod assessment may lack 2023–2025. Missing years
# remain NA and should never be silently imputed as zero.
# ==============================================================================

cod_age1_weight_g <- 142.7
cod_age1_weight_kg <- cod_age1_weight_g / 1000
cod_regional_fraction <- 0.06

# ============================================================
# 0. SETTINGS
# ============================================================

setwd("/Users/alexaei/Jottacloud/projects/CEDAP/data/fish/Mesopredator_Collaborator_Package")
# Portable project-root detection. Supports Source, Rscript and interactive
# RStudio Cmd+Enter. All paths below remain inside the extracted package.
.meso_is_root <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x) &&
    file.exists(file.path(x, "Data", "Inputs", "Metadata_TnT - metadata_main.csv")) &&
    file.exists(file.path(x, "Data", "Inputs", "species_table_n.csv"))
}
.meso_sources <- lapply(sys.frames(), function(frame) frame$ofile)
.meso_sources <- Filter(function(x) is.character(x) && length(x) == 1L &&
                          !is.na(x) && nzchar(x) && file.exists(x), .meso_sources)
.meso_candidates <- character(0)
if (length(.meso_sources)) {
  .meso_candidates <- c(.meso_candidates,
    dirname(normalizePath(tail(.meso_sources, 1L)[[1L]])))
}
.meso_command_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(.meso_command_file)) {
  .meso_candidates <- c(.meso_candidates,
    dirname(normalizePath(sub("^--file=", "", .meso_command_file[[1L]]),
                          mustWork = FALSE)))
}
# RStudio active editor: works when running individual lines with Cmd+Enter.
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  .meso_editor <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                           error = function(e) "")
  if (nzchar(.meso_editor)) {
    .meso_candidates <- c(.meso_candidates, dirname(.meso_editor))
  }
}
# RStudio Project / console fallback, if the extracted folder is working dir.
.meso_candidates <- unique(c(.meso_candidates, getwd(),
                             dirname(getwd())))
.meso_matches <- .meso_candidates[vapply(.meso_candidates,
                                        .meso_is_root, logical(1))]
if (!length(.meso_matches)) {
  stop("Cannot locate the Mesopredator package folder. Open the master R file ",
       "in RStudio (save it first), or set the working directory to the ",
       "extracted Mesopredator_Collaborator_Package folder. ",
       "Required: ../../data/raw/ next to the master script.")
}
project_dir <- normalizePath(.meso_matches[[1L]], mustWork = TRUE)
rm(.meso_sources, .meso_candidates, .meso_matches, .meso_command_file)
required_input_paths <- file.path(project_dir, "Data", "Inputs", c(
  "Metadata_TnT - metadata_main.csv",
  "species_table_n.csv",
  "område og art_data.csv"
))
missing_inputs <- required_input_paths[!file.exists(required_input_paths)]
if (length(missing_inputs)) {
  stop("Missing files in extracted package:\n",
       paste(missing_inputs, collapse = "\n"))
}
final_dir <- file.path(project_dir, "final result v1.4 2012-2025")
legacy_final_dir <- file.path(project_dir, "final result")
for (folder in c("Figures", "Tables", "Documentation", "Inputs"))
  dir.create(file.path(final_dir, folder), recursive=TRUE, showWarnings=FALSE)

metadata_file <- file.path(
  project_dir,
  "Data",
  "Inputs",
  "Metadata_TnT - metadata_main.csv"
)

species_file <- file.path(
  project_dir,
  "Data",
  "Inputs",
  "species_table_n.csv"
)

analysis_years <- 2012:2025

random_seed <- 42
set.seed(random_seed)

models_to_run <- c(
  "LM",
  "Quasibinomial",
  "Hurdle"
)

# Before the first observed year, predictions are held
# at the first observed year's fitted value.
#
# This is an explicit modeling assumption.
# It is NOT evidence of observed diets in earlier years.

hold_pre_observation_constant <- TRUE

required_packages <- c(
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Please install: ",
    paste(missing_packages, collapse = ", ")
  )
}

output_dir <- file.path(
  final_dir,
  "Tables",
  "Model_comparison"
)

figure_dir <- file.path(
  final_dir,
  "Figures",
  "Model_comparison"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

for (model_name in models_to_run) {
  
  dir.create(
    file.path(
      final_dir,
      "Figures",
      model_name
    ),
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ============================================================
# 1. LOAD INPUT DATA
# ============================================================

if (!file.exists(metadata_file)) {
  fallback <- file.path(project_dir,"Data","Inputs","Metadata_TnT - metadata_main(1).csv")
  if (file.exists(fallback)) metadata_file <- fallback
}
if (!file.exists(metadata_file)) {
  stop("Missing metadata file: ", metadata_file)
}

if (!file.exists(species_file)) {
  fallback <- file.path(project_dir,"Data","Inputs","species_table_n(1).csv")
  if (file.exists(fallback)) species_file <- fallback
}
if (!file.exists(species_file)) {
  stop("Missing species file: ", species_file)
}

metadata <- readr::read_csv(
  metadata_file,
  show_col_types = FALSE
)

species <- readr::read_csv(
  species_file,
  show_col_types = FALSE
)

# The original collaborator script uses the first column
# of the species table as library_name.

if (!"library_name" %in% names(species)) {
  
  names(species)[1] <- "library_name"
  
}

required_metadata <- c(
  "library_name",
  "sample_year",
  "sample_month",
  "species"
)

missing_metadata <- setdiff(
  required_metadata,
  names(metadata)
)

if (length(missing_metadata) > 0) {
  
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
  
}

if (!"Gadus_morhua" %in% names(species)) {
  stop("Missing Gadus_morhua column.")
}


# ============================================================
# 2. LIBRARY MATCHING
# ============================================================

clean_library_id <- function(x) {
  
  x <- trimws(as.character(x))
  
  x[x %in% c(
    "",
    "NA",
    "N/A",
    "NULL",
    "NaN"
  )] <- NA_character_
  
  x
  
}

metadata <- metadata |>
  dplyr::mutate(
    library_name = clean_library_id(library_name)
  )

species <- species |>
  dplyr::mutate(
    library_name = clean_library_id(library_name)
  )

n_metadata_missing <- sum(
  is.na(metadata$library_name)
)

n_species_missing <- sum(
  is.na(species$library_name)
)

metadata_for_join <- metadata |>
  dplyr::filter(
    !is.na(library_name)
  )

species_for_join <- species |>
  dplyr::filter(
    !is.na(library_name)
  )

metadata_duplicates <- metadata_for_join |>
  dplyr::count(
    library_name,
    name = "n"
  ) |>
  dplyr::filter(n > 1)

species_duplicates <- species_for_join |>
  dplyr::count(
    library_name,
    name = "n"
  ) |>
  dplyr::filter(n > 1)

if (nrow(metadata_duplicates) > 0) {
  
  print(metadata_duplicates)
  
  stop("Duplicate valid metadata library IDs.")
  
}

if (nrow(species_duplicates) > 0) {
  
  print(species_duplicates)
  
  stop("Duplicate valid species library IDs.")
  
}

metadata_unmatched <- dplyr::anti_join(
  metadata_for_join,
  species_for_join,
  by = "library_name",
  na_matches = "never"
)

species_unmatched <- dplyr::anti_join(
  species_for_join,
  metadata_for_join,
  by = "library_name",
  na_matches = "never"
)

merged <- dplyr::inner_join(
  metadata_for_join,
  species_for_join,
  by = "library_name",
  relationship = "one-to-one",
  na_matches = "never"
)

if (nrow(merged) == 0) {
  stop("No matching libraries found.")
}

n_matched_before_filter <- nrow(merged)

library_qc <- tibble::tibble(
  
  Metric = c(
    "Metadata_total",
    "Metadata_missing_ID",
    "Metadata_valid_ID",
    "Species_total",
    "Species_missing_ID",
    "Species_valid_ID",
    "Matched_libraries",
    "Metadata_unmatched",
    "Species_unmatched"
  ),
  
  Value = c(
    nrow(metadata),
    n_metadata_missing,
    nrow(metadata_for_join),
    nrow(species),
    n_species_missing,
    nrow(species_for_join),
    n_matched_before_filter,
    nrow(metadata_unmatched),
    nrow(species_unmatched)
  )
  
)

readr::write_csv(
  library_qc,
  file.path(
    output_dir,
    "library_matching_QC.csv"
  )
)

readr::write_csv(
  metadata_unmatched,
  file.path(
    output_dir,
    "metadata_unmatched.csv"
  )
)

readr::write_csv(
  species_unmatched,
  file.path(
    output_dir,
    "species_unmatched.csv"
  )
)

cat("\n========== LIBRARY MATCHING ==========\n")
print(library_qc)


# ============================================================
# 3. RECONSTRUCT RELATIVE READ ABUNDANCE
# ============================================================

# Original collaborator processing:
#
# 1. Remove samples with >5 Halichoerus_grypus reads.
# 2. Remove homus/unidentified/unknown columns.
# 3. Exclude the three named host species from prey reads.
# 4. Calculate RRA using pmax(Total_Prey_Reads, 1).
#
# Note:
# Halichoerus_grypus is NOT in the original host exclusion
# list after filtering. We retain this baseline choice for
# reproducibility, rather than silently changing it.

if ("Halichoerus_grypus" %in% names(merged)) {
  
  merged <- merged |>
    dplyr::filter(
      !is.na(Halichoerus_grypus),
      Halichoerus_grypus <= 5
    )
  
}

n_after_host_filter <- nrow(merged)

excluded_cols <- grep(
  "homus|unidentified|unknown",
  names(species_for_join),
  ignore.case = TRUE,
  value = TRUE
)

merged <- merged |>
  dplyr::select(
    -dplyr::any_of(excluded_cols)
  )

host_cols <- c(
  "Phoca_vitulina",
  "Phalacrocorax_carbo",
  "Phalacrocorax_aristotelis"
)

prey_cols <- setdiff(
  names(species_for_join),
  c(
    "library_name",
    host_cols,
    excluded_cols
  )
)

prey_cols <- intersect(
  prey_cols,
  names(merged)
)

if (length(prey_cols) == 0) {
  stop("No prey columns found.")
}

# Ensure that all selected read columns are numeric.
# Stop rather than silently converting unexpected text
# values into missing read counts.

for (col in prey_cols) {
  
  original_values <- merged[[col]]
  
  converted_values <- suppressWarnings(
    as.numeric(original_values)
  )
  
  unexpected_na <- (
    !is.na(original_values) &
      is.na(converted_values)
  )
  
  if (any(unexpected_na)) {
    
    stop(
      "Non-numeric read counts found in column: ",
      col
    )
    
  }
  
  merged[[col]] <- converted_values
  
}

if (any(vapply(
  merged[prey_cols],
  function(x) {
    any(!is.na(x) & x < 0)
  },
  logical(1)
))) {
  
  stop("Negative read counts detected.")
  
}

row_sum_existing <- function(df, cols) {
  
  cols <- intersect(
    cols,
    names(df)
  )
  
  if (length(cols) == 0) {
    return(rep(0, nrow(df)))
  }
  
  rowSums(
    as.data.frame(
      df[, cols, drop = FALSE]
    ),
    na.rm = TRUE
  )
  
}

merged$Total_Prey_Reads <- row_sum_existing(
  merged,
  prey_cols
)

# VERSION 0.4:
# Do not remove zero-read samples.
#
# Their RRA is assigned zero to reproduce the original
# collaborator analysis.
#
# A zero-read sample is not necessarily evidence that
# the corresponding prey was biologically absent.

merged <- merged |>
  dplyr::mutate(
    Zero_Prey_Reads = Total_Prey_Reads == 0
  )

n_zero_read <- sum(
  merged$Zero_Prey_Reads,
  na.rm = TRUE
)

rra_denominator <- pmax(
  merged$Total_Prey_Reads,
  1
)

merged$Cod_RRA <- row_sum_existing(
  merged,
  "Gadus_morhua"
) / rra_denominator

merged$Clupeid_RRA <- row_sum_existing(
  merged,
  c(
    "Clupea_harengus",
    "Clupea_pallasii",
    "Sprattus_sprattus"
  )
) / rra_denominator

merged$Wrasse_RRA <- row_sum_existing(
  merged,
  c(
    "Labrus_bergylta",
    "Symphodus_melops"
  )
) / rra_denominator

cat("\n========== READ FILTERING ==========\n")

cat(
  "Matched libraries before filtering:",
  n_matched_before_filter,
  "\n"
)

cat(
  "Libraries after host-read filtering:",
  n_after_host_filter,
  "\n"
)

cat(
  "Zero-prey-read libraries retained:",
  n_zero_read,
  "\n"
)


# ============================================================
# 4. CLASSIFY PREDATOR, YEAR AND SEASON
# ============================================================

map_season <- function(x) {
  
  x <- tolower(
    trimws(as.character(x))
  )
  
  dplyr::case_when(
    
    x %in% c(
      "december", "january", "february",
      "dec", "jan", "feb",
      "12", "1", "2"
    ) ~ "Winter",
    
    x %in% c(
      "march", "april", "may",
      "mar", "apr", "aprl",
      "3", "4", "5"
    ) ~ "Spring",
    
    x %in% c(
      "june", "july", "august",
      "jun", "jul", "aug",
      "6", "7", "8"
    ) ~ "Summer",
    
    x %in% c(
      "september", "october", "november",
      "sep", "oct", "nov",
      "9", "10", "11",
      "august_september"
    ) ~ "Autumn",
    
    TRUE ~ NA_character_
    
  )
  
}

# Preserve every library long enough to document
# why any remaining samples are excluded.

classified_data <- merged |>
  dplyr::mutate(
    
    Year = suppressWarnings(
      as.integer(sample_year)
    ),
    
    Season = "All_year",
    
    Predator = dplyr::case_when(
      
      species == "harbor_seal" ~ "Seal",
      
      species %in% c(
        "cormorant",
        "top_cormorant"
      ) ~ "Cormorant",
      
      TRUE ~ NA_character_
      
    ),
    
    Inclusion_status = dplyr::case_when(
      
      is.na(Predator) ~
        "Excluded: predator classification",
      
      is.na(Year) ~
        "Excluded: missing year",
      
      TRUE ~
        "Included"
      
    )
    
  )

sample_exclusion_qc <- classified_data |>
  dplyr::count(
    Inclusion_status,
    name = "N_libraries"
  )

readr::write_csv(
  sample_exclusion_qc,
  file.path(
    output_dir,
    "sample_exclusion_QC.csv"
  )
)

readr::write_csv(
  classified_data |>
    dplyr::select(
      library_name,
      sample_year,
      sample_month,
      species,
      Predator,
      Year,
      Season,
      Total_Prey_Reads,
      Zero_Prey_Reads,
      Inclusion_status
    ),
  file.path(
    output_dir,
    "sample_inclusion_details.csv"
  )
)

diet_data <- classified_data |>
  dplyr::filter(
    Inclusion_status == "Included"
  ) |>
  dplyr::mutate(
    Season = factor(
      Season,
      levels = "All_year"
    )
  ) |>
  dplyr::select(
    library_name,
    Predator,
    Year,
    Season,
    Total_Prey_Reads,
    Zero_Prey_Reads,
    Cod_RRA,
    Clupeid_RRA,
    Wrasse_RRA
  ) |>
  tidyr::pivot_longer(
    cols = c(
      Cod_RRA,
      Clupeid_RRA,
      Wrasse_RRA
    ),
    names_to = "Prey",
    values_to = "RRA"
  ) |>
  dplyr::mutate(
    Prey = sub(
      "_RRA$",
      "",
      Prey
    )
  )

if (nrow(diet_data) == 0) {
  stop("No diet data after classification.")
}

if (any(
  !is.finite(diet_data$RRA) |
  diet_data$RRA < 0 |
  diet_data$RRA > 1
)) {
  stop("Invalid RRA values detected.")
}



# v1.4: explicitly verify all 2025 cormorant libraries are retained.
qc_cormorant_2025 <- classified_data |>
  dplyr::filter(Year == 2025L, Predator == "Cormorant") |>
  dplyr::select(library_name, species, Year, Inclusion_status, Total_Prey_Reads)
readr::write_csv(qc_cormorant_2025,
  file.path(output_dir, "cormorant_2025_inclusion_QC.csv"))
if (nrow(qc_cormorant_2025) != 3L ||
    any(qc_cormorant_2025$Inclusion_status != "Included") ||
    anyDuplicated(qc_cormorant_2025$library_name)) {
  stop("v1.4 QC: expected three distinct INCLUDED cormorant libraries in 2025")
}
message("v1.4 QC: three 2025 cormorant libraries included: ",
        paste(qc_cormorant_2025$library_name, collapse = ", "))
sample_qc <- diet_data |>
  dplyr::distinct(
    library_name,
    Predator
  ) |>
  dplyr::count(
    Predator,
    name = "N_libraries"
  )

sample_year_qc <- diet_data |>
  dplyr::distinct(
    library_name,
    Predator,
    Year,
    Season,
    Zero_Prey_Reads
  ) |>
  dplyr::group_by(
    Predator,
    Year
  ) |>
  dplyr::summarise(
    N_libraries = dplyr::n(),
    N_zero_prey_reads = sum(Zero_Prey_Reads),
    .groups = "drop"
  )

readr::write_csv(
  sample_qc,
  file.path(
    output_dir,
    "diet_sample_QC.csv"
  )
)

readr::write_csv(
  sample_year_qc,
  file.path(
    output_dir,
    "diet_sample_year_QC.csv"
  )
)

readr::write_csv(
  diet_data,
  file.path(
    output_dir,
    "RRA_input_data.csv"
  )
)

cat("\n========== SAMPLE QC ==========\n")

print(sample_qc)
print(sample_exclusion_qc)

cat("\n========== YEAR COVERAGE ==========\n")

print(
  sample_year_qc,
  n = Inf
)

# Compare with previous QC 2.
#
# These are reference values, not hard-coded sample
# selections or values imposed on the analysis.

expected_counts <- tibble::tibble(
  Predator = c("Seal", "Cormorant"),
  Expected_N = c(60L, 75L)
)

sample_count_comparison <- expected_counts |>
  dplyr::left_join(
    sample_qc,
    by = "Predator"
  ) |>
  dplyr::mutate(
    Difference = N_libraries - Expected_N
  )

readr::write_csv(
  sample_count_comparison,
  file.path(
    output_dir,
    "sample_count_comparison_QC.csv"
  )
)

cat("\n========== QC 2 COMPARISON ==========\n")
print(sample_count_comparison)

if (
  any(
    is.na(sample_count_comparison$Difference) |
    sample_count_comparison$Difference != 0
  )
) {
  
  warning(
    "Sample counts differ from QC 2. ",
    "Inspect sample_inclusion_details.csv ",
    "before interpreting the model results."
  )
  
}


# ============================================================
# 5. FIT DIET MODELS
# ============================================================

fit_one_group <- function(
    df,
    model_name
) {
  
  first_year <- min(df$Year)
  
  df <- df |>
    dplyr::mutate(
      Year_c = Year - first_year,
      Season = droplevels(Season)
    )
  
  has_season <- FALSE  # Annual model: season intentionally excluded
  
  formula_full <- if (has_season) {
    RRA ~ Year_c + Season
  } else {
    RRA ~ Year_c
  }
  
  if (model_name == "LM") {
    
    model <- stats::lm(
      formula_full,
      data = df
    )
    
    return(list(
      type = "LM",
      model = model,
      first_year = first_year,
      has_season = has_season,
      season_levels = levels(df$Season),
      fallback = FALSE
    ))
    
  }
  
  if (model_name == "Quasibinomial") {
    
    model <- stats::glm(
      formula_full,
      family = stats::quasibinomial(
        link = "logit"
      ),
      data = df
    )
    
    return(list(
      type = "Quasibinomial",
      model = model,
      first_year = first_year,
      has_season = has_season,
      season_levels = levels(df$Season),
      fallback = FALSE
    ))
    
  }
  
  if (model_name == "Hurdle") {
    
    df <- df |>
      dplyr::mutate(
        Detected = as.integer(RRA > 0)
      )
    
    positive <- df |>
      dplyr::filter(
        Detected == 1
      )
    
    # For sparse groups, use the constant-RRA fallback.
    # This is recorded explicitly in the model register.
    
    if (
      nrow(positive) < 5 ||
      dplyr::n_distinct(df$Detected) < 2
    ) {
      
      fallback_model <- stats::glm(
        RRA ~ 1,
        family = stats::quasibinomial(),
        data = df
      )
      
      return(list(
        type = "Hurdle",
        model = fallback_model,
        first_year = first_year,
        has_season = FALSE,
        season_levels = NULL,
        fallback = TRUE
      ))
      
    }
    
    detection_model <- stats::glm(
      Detected ~ Year_c,
      family = stats::binomial(),
      data = df
    )
    
    positive_model <- stats::glm(
      RRA ~ 1,
      family = stats::quasibinomial(),
      data = positive
    )
    
    return(list(
      type = "Hurdle",
      detection = detection_model,
      positive = positive_model,
      first_year = first_year,
      has_season = FALSE,
      season_levels = NULL,
      fallback = FALSE
    ))
    
  }
  
  stop("Unknown model: ", model_name)
  
}


# ============================================================
# 6. CALENDAR-YEAR PREDICTIONS
# ============================================================

predict_one_group <- function(
    fit,
    df,
    model_name
) {
  
  first_observed <- min(df$Year)
  last_observed <- max(df$Year)
  
  season_weights <- prop.table(
    table(
      droplevels(df$Season)
    )
  )
  
  predict_at_year <- function(year) {
    
    effective_year <- year
    
    if (
      hold_pre_observation_constant &&
      year < first_observed
    ) {
      
      effective_year <- first_observed
      
    }
    
    if (
      fit$has_season &&
      model_name != "Hurdle"
    ) {
      
      newdata <- data.frame(
        Year_c = rep(
          effective_year - fit$first_year,
          length(season_weights)
        ),
        Season = factor(
          names(season_weights),
          levels = fit$season_levels
        )
      )
      
      predicted_values <- stats::predict(
        fit$model,
        newdata = newdata,
        type = "response"
      )
      
      result <- sum(
        predicted_values *
          as.numeric(season_weights)
      )
      
    } else {
      
      newdata <- data.frame(
        Year_c = effective_year - fit$first_year
      )
      
      if (
        model_name == "Hurdle" &&
        !fit$fallback
      ) {
        
        detection <- stats::predict(
          fit$detection,
          newdata = newdata,
          type = "response"
        )
        
        positive <- stats::predict(
          fit$positive,
          newdata = newdata,
          type = "response"
        )
        
        result <- detection * positive
        
      } else {
        
        result <- stats::predict(
          fit$model,
          newdata = newdata,
          type = "response"
        )
        
      }
      
    }
    
    as.numeric(result)
    
  }
  
  raw_predictions <- vapply(
    analysis_years,
    predict_at_year,
    numeric(1)
  )
  
  bounded_predictions <- pmin(
    pmax(raw_predictions, 0),
    1
  )
  
  tibble::tibble(
    
    Year = analysis_years,
    
    Model = model_name,
    
    Predicted_RRA_raw = raw_predictions,
    
    Predicted_RRA = bounded_predictions,
    
    Prediction_adjusted = abs(
      bounded_predictions -
        raw_predictions
    ) > 1e-10,
    
    Period = dplyr::case_when(
      
      Year < first_observed ~
        "Before observations",
      
      Year > last_observed ~
        "After observations",
      
      TRUE ~
        "Within observed range"
      
    ),
    
    First_observed_year = first_observed,
    
    Last_observed_year = last_observed,
    
    Fallback_used = fit$fallback
    
  )
  
}


# ============================================================
# 7. RUN ALL 18 MODEL COMBINATIONS
# ============================================================

groups <- diet_data |>
  dplyr::group_by(
    Predator,
    Prey
  ) |>
  dplyr::group_split()

all_predictions <- list()
model_register <- list()
fitted_models <- list()

for (df in groups) {
  
  predator <- unique(df$Predator)
  prey <- unique(df$Prey)
  
  for (model_name in models_to_run) {
    
    label <- paste(
      predator,
      prey,
      model_name,
      sep = "_"
    )
    
    message("Fitting: ", label)
    
    result <- tryCatch({
      
      fit <- fit_one_group(
        df,
        model_name
      )
      
      predictions <- predict_one_group(
        fit,
        df,
        model_name
      ) |>
        dplyr::mutate(
          Predator = predator,
          Prey = prey,
          .before = 1
        )
      
      list(
        fit = fit,
        predictions = predictions,
        error = NA_character_
      )
      
    }, error = function(e) {
      
      list(
        fit = NULL,
        predictions = NULL,
        error = conditionMessage(e)
      )
      
    })
    
    model_register[[label]] <- tibble::tibble(
      
      Predator = predator,
      
      Prey = prey,
      
      Model = model_name,
      
      N = nrow(df),
      
      N_positive = sum(
        df$RRA > 0
      ),
      
      N_zero_prey_read_libraries = sum(
        df$Zero_Prey_Reads
      ),
      
      First_year = min(df$Year),
      
      Last_year = max(df$Year),
      
      Status = if (
        is.null(result$fit)
      ) {
        "FAILED"
      } else {
        "SUCCESS"
      },
      
      Fallback = if (
        is.null(result$fit)
      ) {
        NA
      } else {
        result$fit$fallback
      },
      
      Error = result$error
      
    )
    
    if (!is.null(result$fit)) {
      
      fitted_models[[label]] <- result$fit
      
      all_predictions[[label]] <- result$predictions
      
    }
    
  }
  
}

model_register <- dplyr::bind_rows(
  model_register
)

readr::write_csv(
  model_register,
  file.path(
    output_dir,
    "diet_model_register.csv"
  )
)

if (length(all_predictions) == 0) {
  
  stop(
    "All model fits failed. ",
    "Inspect diet_model_register.csv."
  )
  
}

diet_predictions <- dplyr::bind_rows(
  all_predictions
)

readr::write_csv(
  diet_predictions,
  file.path(
    output_dir,
    "diet_model_predictions.csv"
  )
)

saveRDS(
  fitted_models,
  file.path(
    output_dir,
    "fitted_diet_models.rds"
  )
)


# Original plot figures are generated in the final section below.

# ============================================================
# 10. FINAL QC
# ============================================================

cat("\n========================================\n")
cat("MESOPREDATOR DIET ANALYSIS — VERSION 1.0 ANNUAL\n")
cat("========================================\n")

cat(
  "Matched libraries:",
  n_matched_before_filter,
  "\n"
)

cat(
  "Libraries after host-read filtering:",
  n_after_host_filter,
  "\n"
)

cat(
  "Zero-prey-read libraries retained:",
  n_zero_read,
  "\n"
)

cat(
  "Models attempted:",
  nrow(model_register),
  "\n"
)

cat(
  "Successful models:",
  sum(
    model_register$Status == "SUCCESS"
  ),
  "\n"
)

cat(
  "Failed models:",
  sum(
    model_register$Status == "FAILED"
  ),
  "\n"
)

cat(
  "Fallback models:",
  sum(
    model_register$Fallback,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Predictions adjusted to [0, 1]:",
  sum(
    diet_predictions$Prediction_adjusted,
    na.rm = TRUE
  ),
  "\n"
)

cat("\n========== FINAL SAMPLE COUNTS ==========\n")

print(sample_qc)

cat("\n========== QC 2 COMPARISON ==========\n")

print(sample_count_comparison)

cat("\n========== MODEL REGISTER ==========\n")

print(
  model_register,
  n = Inf
)

cat("\n========== OUTPUT LOCATIONS ==========\n")

cat(
  "Tables:",
  output_dir,
  "\n"
)

cat(
  "Figures:",
  figure_dir,
  "\n"
)

cat("\nNEXT STEP:\n")

cat(
  "Validate sample counts and model predictions ",
  "before integrating ICES and Monte Carlo simulations.\n"
)

cat(
  "Diet prediction point estimates alone do not ",
  "represent full parameter uncertainty.\n"
)

cat(
  "RRA is not directly calibrated to biomass ",
  "consumption proportions.\n"
)






# ==============================================================================
# PART II — MONTE CARLO AND PREY BIOMASS INTEGRATION
# ==============================================================================

# ================================================================
# 0. SETTINGS
# ================================================================

# project_dir inherited from unified settings

years <- analysis_years
n_iter <- 5000L
seed <- random_seed

diet_models <- c("LM", "Quasibinomial", "Hurdle")
predators <- c("Seal", "Cormorant")
prey_groups <- c("Cod", "Clupeid", "Wrasse")

input_dir <- file.path(project_dir, "Data", "Inputs")
model_dir <- file.path(final_dir, "Tables", "Model_comparison")

out_dir <- file.path(final_dir, "Tables", "Monte_Carlo")
fig_dir <- file.path(final_dir, "Figures", "Monte_Carlo")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

packages <- c(
  "dplyr", "tidyr", "readr", "ggplot2",
  "tibble", "MASS", "icesSAG", "stringr"
)

missing_packages <- packages[
  !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  stop(
    "Install missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

set.seed(seed)


# ================================================================
# 1. LOAD VALIDATED DIET MODELS
# ================================================================

model_file <- file.path(model_dir, "fitted_diet_models.rds")
diet_file <- file.path(model_dir, "RRA_input_data.csv")

stopifnot(
  file.exists(model_file),
  file.exists(diet_file)
)

fitted_models <- readRDS(model_file)

diet_data <- readr::read_csv(
  diet_file,
  show_col_types = FALSE
)

diet_data <- diet_data |>
  dplyr::mutate(
    Season = factor(
      Season,
      levels = "All_year"
    )
  )

expected_names <- as.vector(
  outer(
    as.vector(
      outer(predators, prey_groups, paste, sep = "_")
    ),
    diet_models,
    paste,
    sep = "_"
  )
)

missing_models <- setdiff(
  expected_names,
  names(fitted_models)
)

if (length(missing_models)) {
  stop(
    "Missing fitted models: ",
    paste(missing_models, collapse = ", ")
  )
}

message("Loaded ", length(fitted_models), " fitted diet models.")


# ================================================================
# 2. MODEL-PARAMETER SIMULATION
# ================================================================
#
# Draw coefficients from approximate multivariate normal
# sampling distributions using coef() and vcov().
#
# For LM: includes residual-scale variability only through
# parameter uncertainty, not extra observation-level noise.
#
# For quasibinomial: vcov() reflects estimated dispersion.
#
# For hurdle: detection and positive-component parameters
# are sampled separately.
#
# Positive-component and detection-model covariance is
# treated as independent in this development version.
# ================================================================

draw_coef <- function(model, n) {
  
  beta <- stats::coef(model)
  V <- stats::vcov(model)
  
  if (
    any(!is.finite(beta)) ||
    any(!is.finite(V))
  ) {
    stop("Non-finite model coefficients or covariance.")
  }
  
  if (length(beta) == 1L) {
    
    sd_beta <- sqrt(max(as.numeric(V), 0))
    
    return(
      matrix(
        stats::rnorm(n, beta, sd_beta),
        ncol = 1,
        dimnames = list(NULL, names(beta))
      )
    )
  }
  
  # Numerical tolerance for near-singular covariance matrices.
  eigenvalues <- eigen(
    (V + t(V)) / 2,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  
  if (min(eigenvalues) < -1e-7) {
    stop("Coefficient covariance matrix is not PSD.")
  }
  
  V <- (V + t(V)) / 2
  
  if (min(eigenvalues) < 0) {
    eig <- eigen(V, symmetric = TRUE)
    eig$values <- pmax(eig$values, 0)
    V <- eig$vectors %*%
      diag(eig$values, nrow = length(eig$values)) %*%
      t(eig$vectors)
  }
  
  result <- MASS::mvrnorm(
    n = n,
    mu = beta,
    Sigma = V
  )
  
  result <- matrix(
    result,
    nrow = n,
    ncol = length(beta)
  )
  
  colnames(result) <- names(beta)
  result
}

# Build model matrices with the model's original terms
# and contrasts, ensuring correct season coefficients.

model_matrix_for <- function(model, newdata) {
  
  terms_no_response <- stats::delete.response(
    stats::terms(model)
  )
  
  stats::model.matrix(
    terms_no_response,
    data = newdata,
    contrasts.arg = model$contrasts,
    xlev = model$xlevels
  )
}

draw_single_diet <- function(
    fit,
    group_data,
    year,
    n
) {
  
  first_observed <- min(group_data$Year)
  
  effective_year <- max(
    year,
    first_observed
  )
  
  year_c <- effective_year - fit$first_year
  
  if (fit$type == "Hurdle" && !fit$fallback) {
    
    newdata <- data.frame(Year_c = year_c)
    
    X_det <- model_matrix_for(
      fit$detection,
      newdata
    )
    
    X_pos <- model_matrix_for(
      fit$positive,
      newdata
    )
    
    beta_det <- draw_coef(fit$detection, n)
    beta_pos <- draw_coef(fit$positive, n)
    
    eta_det <- as.vector(
      beta_det %*% t(X_det)
    )
    
    eta_pos <- as.vector(
      beta_pos %*% t(X_pos)
    )
    
    p_detect <- stats::plogis(eta_det)
    positive_rra <- stats::plogis(eta_pos)
    
    return(
      pmin(pmax(p_detect * positive_rra, 0), 1)
    )
  }
  
  model <- fit$model
  beta <- draw_coef(model, n)
  
  if (fit$has_season && fit$type != "Hurdle") {
    
    season_counts <- table(
      droplevels(group_data$Season)
    )
    
    season_weights <- as.numeric(
      season_counts / sum(season_counts)
    )
    
    newdata <- data.frame(
      Year_c = rep(year_c, length(season_weights)),
      Season = factor(
        names(season_counts),
        levels = fit$season_levels
      )
    )
    
    X <- model_matrix_for(model, newdata)
    
    eta <- beta %*% t(X)
    
    if (fit$type == "LM") {
      predictions <- eta
    } else {
      predictions <- stats::plogis(eta)
    }
    
    result <- as.vector(
      predictions %*% season_weights
    )
    
  } else {
    
    newdata <- data.frame(Year_c = year_c)
    X <- model_matrix_for(model, newdata)
    
    eta <- as.vector(beta %*% t(X))
    
    if (fit$type == "LM") {
      result <- eta
    } else {
      result <- stats::plogis(eta)
    }
  }
  
  pmin(pmax(result, 0), 1)
}

draw_diet_composition <- function(
    predator,
    model_name,
    year,
    n
) {
  
  draws <- lapply(prey_groups, function(prey) {
    
    key <- paste(
      predator,
      prey,
      model_name,
      sep = "_"
    )
    
    group_data <- diet_data |>
      dplyr::filter(
        Predator == predator,
        Prey == prey
      )
    
    draw_single_diet(
      fitted_models[[key]],
      group_data,
      year,
      n
    )
  })
  
  composition <- do.call(cbind, draws)
  colnames(composition) <- prey_groups
  
  # Ensure the three modeled prey groups plus "Other"
  # form a valid composition for each simulation draw.
  #
  # If modeled fractions exceed 1, normalize them to 1.
  # Otherwise, retain the residual as Other.
  #
  # This is a scenario transformation of RRA, not
  # empirical biomass calibration.
  
  modeled_sum <- rowSums(composition)
  
  over_one <- modeled_sum > 1
  
  if (any(over_one)) {
    composition[over_one, ] <-
      composition[over_one, , drop = FALSE] /
      modeled_sum[over_one]
  }
  
  other <- pmax(
    1 - rowSums(composition),
    0
  )
  
  composition <- cbind(
    composition,
    Other = other
  )
  
  if (any(
    abs(rowSums(composition) - 1) > 1e-8
  )) {
    stop("Invalid diet composition.")
  }
  
  composition
}


# ================================================================
# 3. ICES DATA — DIRECT DOWNLOAD
# ================================================================
#
# Retain original stock codes.
# Save assessment provenance and raw downloaded tables.
#
# Assessment availability and stock coverage can change.
# The original cod assessment is not silently replaced
# with a different stock definition.
# ================================================================

get_latest_sag <- function(
    stock_code,
    target_year,
    lookback = 3L
) {
  
  for (assessment_year in seq(
    target_year,
    target_year - lookback,
    by = -1L
  )) {
    
    keys <- tryCatch(
      icesSAG::findAssessmentKey(
        stock_code,
        year = assessment_year
      ),
      error = function(e) integer(0)
    )
    
    if (!length(keys)) next
    
    for (key in keys) {
      
      tbl <- tryCatch(
        icesSAG::getSummaryTable(key),
        error = function(e) NULL
      )
      
      if (is.data.frame(tbl) && nrow(tbl)) {
        
        return(list(
          data = tbl,
          stock_code = stock_code,
          assessment_year = assessment_year,
          assessment_key = as.character(key)
        ))
      }
    }
  }
  
  stop(
    "No ICES assessment found for ",
    stock_code,
    ". Check stock code and API availability."
  )
}

stock_specs <- tibble::tribble(
  ~Stock,          ~Code,           ~Target_year, ~Ratio, ~Ratio_CV, ~Prefer_total,
  "Sprat",         "spr.27.3a4",     2026L,        0.15,   0.30,      TRUE,
  "NSAS_Herring",  "her.27.3a47d",   2026L,        0.10,   0.20,      TRUE,
  "WBSS_Herring",  "her.27.20-24",   2026L,        0.40,   0.20,      TRUE,
  "Cod_Adults",    "cod.27.47d20",   2025L,        0.06,   0.25,      FALSE
)

ices_results <- vector(
  "list",
  nrow(stock_specs)
)

names(ices_results) <- stock_specs$Stock

for (i in seq_len(nrow(stock_specs))) {
  
  spec <- stock_specs[i, ]
  
  message(
    "Downloading ICES: ",
    spec$Stock,
    " (",
    spec$Code,
    ")"
  )
  
  downloaded <- get_latest_sag(
    stock_code = spec$Code,
    target_year = spec$Target_year
  )
  
  ices_results[[spec$Stock]] <- downloaded
  
  readr::write_csv(
    downloaded$data,
    file.path(
      out_dir,
      paste0("ICES_raw_", spec$Stock, ".csv")
    )
  )
}

ices_provenance <- dplyr::bind_rows(
  lapply(names(ices_results), function(stock_name) {
    x <- ices_results[[stock_name]]
    tibble::tibble(
      Stock = stock_name,
      Stock_code = x$stock_code,
      Assessment_year = x$assessment_year,
      Assessment_key = x$assessment_key,
      Rows = nrow(x$data)
    )
  })
)

readr::write_csv(
  ices_provenance,
  file.path(out_dir, "ICES_assessment_provenance.csv")
)


# ================================================================
# 4. ICES BIOMASS EXTRACTION
# ================================================================

get_numeric_column <- function(df, candidates) {
  
  matches <- candidates[
    tolower(candidates) %in% tolower(names(df))
  ]
  
  if (!length(matches)) {
    return(rep(NA_real_, nrow(df)))
  }
  
  actual_name <- names(df)[
    match(
      tolower(matches[1]),
      tolower(names(df))
    )
  ]
  
  suppressWarnings(
    as.numeric(df[[actual_name]])
  )
}

extract_stock <- function(
    downloaded,
    spec
) {
  
  df <- downloaded$data
  
  year <- get_numeric_column(df, "Year")
  
  total <- get_numeric_column(
    df,
    c("StockSize", "TotalBiomass")
  )
  
  total_low <- get_numeric_column(
    df,
    c("low_StockSize", "low_TotalBiomass")
  )
  
  total_high <- get_numeric_column(
    df,
    c("high_StockSize", "high_TotalBiomass")
  )
  
  ssb <- get_numeric_column(df, "SSB")
  ssb_low <- get_numeric_column(df, "low_SSB")
  ssb_high <- get_numeric_column(df, "high_SSB")
  
  use_total <- isTRUE(spec$Prefer_total) &&
    sum(is.finite(total)) > 3
  
  if (use_total) {
    
    biomass <- total
    low <- total_low
    high <- total_high
    metric <- "TSB"
    
  } else {
    
    biomass <- ssb
    low <- ssb_low
    high <- ssb_high
    metric <- "SSB"
    
  }
  
  # Original approximation: (upper - lower) / 3.29
  # where both interval limits exist.
  #
  # Otherwise use a 20% biomass CV assumption.
  
  sd_biomass <- dplyr::if_else(
    is.finite(low) &
      is.finite(high) &
      high > low,
    (high - low) / 3.29,
    biomass * 0.20
  )
  
  tibble::tibble(
    Year = as.integer(year),
    Stock = spec$Stock,
    Biomass_metric = metric,
    Assessment_key = downloaded$assessment_key,
    Combined_biomass_t = biomass,
    Combined_biomass_sd_t = sd_biomass,
    Regional_ratio = spec$Ratio,
    Regional_ratio_cv = spec$Ratio_CV
  ) |>
    dplyr::filter(
      Year %in% years
    )
}

stock_components <- dplyr::bind_rows(
  lapply(seq_len(nrow(stock_specs)), function(i) {
    
    spec <- stock_specs[i, ]
    
    extract_stock(
      ices_results[[spec$Stock]],
      spec
    )
  })
)

# Age-1 cod recruitment biomass proxy: DATRAS mean weight.
#
# Recruitment is assumed to represent thousands of fish,
# with a constant DATRAS-derived mean mass of 0.1427 kg.
#
# Original equation:
# Recruit (thousands) * 0.1427 kg = tonnes.
#
# Confirm recruitment units against the chosen assessment.

cod_df <- ices_results[["Cod_Adults"]]$data

cod_year <- get_numeric_column(
  cod_df,
  "Year"
)

cod_recruit <- get_numeric_column(
  cod_df,
  c("Recruitment", "Recruit")
)

cod_juvenile <- tibble::tibble(
  Year = as.integer(cod_year),
  Stock = "Cod_Juveniles",
  Biomass_metric = "Age1_recruitment_derived_142.7g",
  Assessment_key =
    ices_results[["Cod_Adults"]]$assessment_key,
  Combined_biomass_t = cod_recruit * cod_age1_weight_kg,
  Combined_biomass_sd_t = cod_recruit * cod_age1_weight_kg * 0.30,
  Regional_ratio = ifelse(as.integer(cod_year) >= 2023L, 0.06 * 679300 / 916060, cod_regional_fraction),
  Regional_ratio_cv = 0.25
) |>
  dplyr::filter(
    Year %in% years
  )

stock_components <- dplyr::bind_rows(
  stock_components,
  cod_juvenile
)

# Historical sensitivity comparison removed from final pipeline.
# ICES recruitment must be in THOUSANDS of age-1 fish; verify in
# assessment metadata before interpreting these biomass values.
# Do not fill missing ICES years automatically.
# Missing stock biomass must remain missing.

readr::write_csv(
  stock_components,
  file.path(out_dir, "ICES_stock_components.csv")
)


# ================================================================
# 4B. DOWNLOAD RECRUITMENT FROM ICES FOR 2011--2025
# ================================================================
# Historical assessment: 17652 (cod.27.47d20)
# New substock assessments: 19663 Viking, 19661 Northwestern,
# 19662 Southern. These are different assessment definitions.
# We combine three new substocks for 2023--25 ONLY as a scenario.
# Raw ICES tables and provenance are saved; API failure uses supplied
# archived summary tables ONLY if their AssessmentKey matches.
assessment_keys <- c(Historical=17652L,Viking=19663L,
                     Northwestern=19661L,Southern=19662L)
read_ices_assessment <- function(name, key) {
  target <- file.path(final_dir,'Tables','ICES',paste0('ICES_',name,'_',key,'.csv'))
  dir.create(dirname(target),recursive=TRUE,showWarnings=FALSE)
  dat <- tryCatch(icesSAG::getSummaryTable(key),error=function(e) {
    message('ICES API failed for ',name,': ',conditionMessage(e)); NULL
  })
  if (is.null(dat) || !is.data.frame(dat) || !nrow(dat)) {
    archive <- file.path(legacy_final_dir,'Inputs',paste0('summary_',key,'.csv'))
    if (file.exists(archive)) {
      message('Using archived ICES assessment ',archive)
      dat <- readr::read_csv(archive,show_col_types=FALSE)
    } else stop('Cannot retrieve ICES assessment ',key,
                '. Check internet, or supply archived ',archive)
  }
  if (!all(c('Year','recruitment','recruitment_age','AssessmentKey') %in% names(dat)))
    stop('Missing required ICES fields for ',name)
  if (any(unique(stats::na.omit(as.integer(dat$recruitment_age))) != 1L))
    stop('Recruitment is not age 1 for ',name)
  if (!all(unique(stats::na.omit(as.integer(dat$AssessmentKey))) == key))
    stop('Assessment key mismatch for ',name)
  dat$Year <- as.integer(dat$Year)
  dat$recruitment <- as.numeric(dat$recruitment)
  if (anyDuplicated(dat$Year)) stop('Duplicate ICES year in ',name)
  readr::write_csv(dat,target)
  dat
}
ices_cod <- Map(read_ices_assessment,names(assessment_keys),assessment_keys)
names(ices_cod) <- names(assessment_keys)
recruit_rows <- list()
for (y in 2011:2025) {
  components <- if (y <= 2022) 'Historical' else c('Viking','Northwestern','Southern')
  for (component in components) {
    row <- ices_cod[[component]][ices_cod[[component]]$Year == y,,drop=FALSE]
    if (nrow(row)!=1L || !is.finite(row$recruitment) || row$recruitment < 0)
      stop('Missing/invalid ICES recruitment for ',component,' ',y)
    recruit_rows[[length(recruit_rows)+1L]] <- data.frame(
      Year=y, Stock_component=component,
      Recruitment_thousands=row$recruitment,
      Assessment_key=assessment_keys[[component]],
      Assessment_year=if ('AssessmentYear' %in% names(row)) row$AssessmentYear else NA_integer_)
  }
}
verified <- dplyr::bind_rows(recruit_rows)
readr::write_csv(verified,file.path(final_dir,'Tables','ICES','cod_recruitment_components.csv'))
recruit_annual <- verified |>
  dplyr::group_by(Year) |>
  dplyr::summarise(Recruitment_thousands=sum(Recruitment_thousands),
                   Components=paste(Stock_component,collapse='+'),.groups='drop')
# Spatial ratio is an AREA proxy, not a validated fish-density fraction.
old_area_km2 <- 679300
additional_area_km2 <- 236759
f_old <- 0.06
f_new <- f_old * old_area_km2/(old_area_km2+additional_area_km2)
regional_fraction <- function(y) ifelse(y>=2023,f_new,f_old)
# Literature-informed age-1 to age-2 cohort survival (Risør/Skagerrak).
# Assumed total instantaneous mortality Z = 0.55 yr^-1;
# S = exp(-Z) = 0.57695. See Documentation for citation and limitations.
survival <- exp(-0.55)
W1 <- 0.1427
W2 <- 0.6938
cod_cohorts <- tibble::tibble(Year=years) |>
  dplyr::left_join(recruit_annual |> dplyr::select(Year,R1=Recruitment_thousands),by='Year') |>
  dplyr::left_join(recruit_annual |>
     dplyr::transmute(Year=Year+1L,R_prev=Recruitment_thousands),by='Year') |>
  dplyr::mutate(
    Spatial_fraction_age1=regional_fraction(Year),
    Spatial_fraction_age2_cohort=regional_fraction(Year-1L),
    Survival_age1_to_age2=survival,
    Age1_t=R1*W1*Spatial_fraction_age1,
    Age2_t=R_prev*survival*W2*Spatial_fraction_age2_cohort,
    Age1_2_t=Age1_t+Age2_t,
    Recruitment_source=dplyr::if_else(Year<=2022,'Historical assessment',
                                      'Sum of three new substocks'))
if (any(!is.finite(cod_cohorts$Age1_2_t))) stop('Incomplete cod cohort estimates')
readr::write_csv(cod_cohorts,file.path(final_dir,'Tables','ICES','cod_age1_age2_2012_2025.csv'))
readr::write_csv(data.frame(Original_fraction=f_old,Old_area_km2=old_area_km2,
 Additional_area_km2=additional_area_km2,New_fraction=f_new),
 file.path(final_dir,'Tables','ICES','spatial_fraction_provenance.csv'))
# Replace original age-1-only stock proxy for ALL MC models.
stock_components <- stock_components |> dplyr::filter(Stock!='Cod_Juveniles')
new_cod <- cod_cohorts |> dplyr::transmute(
 Year,Stock='Cod_Juveniles',Biomass_metric='Age1_plus_age2_cohort_scenario',
 Assessment_key='17652_to_19661_19662_19663',
 Combined_biomass_t=Age1_2_t,
 Combined_biomass_sd_t=Age1_2_t*0.30,
 Regional_ratio=1,Regional_ratio_cv=0)
stock_components <- dplyr::bind_rows(stock_components,new_cod)
# Age-specific biomass components retained for transparent QC.
readr::write_csv(cod_cohorts |> dplyr::select(Year,Age1_t,Age2_t,Age1_2_t),
 file.path(final_dir,'Tables','ICES','cod_age_components_QC.csv'))


# ================================================================
# 4C. ADULT COD SSB — ASSESSMENT TRANSITION (NEW IN v1.3)
# ================================================================
# PREVIOUS: Adult cod used only cod.27.47d20 (original ICES assessment,
# key 16999 in the checked output), with no SSB after 2022.
# NOW: Preserve the original assessment SSB for 2012--2022 and use the
# SUM of Viking (19663), Northwestern (19661) and Southern (19662) SSB
# for 2023--2025. This parallels the juvenile assessment transition.
# IMPORTANT: The combined series is a SCENARIO, not a homogeneous ICES
# assessment. The area factor is a proxy, not observed adult-cod density.
# ================================================================

adult_old_factor <- cod_regional_fraction  # 0.06, unchanged
adult_new_factor <- f_new                  # 0.04449277..., same as juveniles
adult_old <- stock_components |>
  dplyr::filter(Stock == "Cod_Adults", Year %in% years) |>
  dplyr::select(Year, Combined_biomass_t, Combined_biomass_sd_t,
                Regional_ratio, Assessment_key)

# Check assessment metadata and extract each component. All stockSizeUnits
# must be tonnes; SSB and year must exist in each downloaded assessment.
read_adult_component <- function(name) {
  z <- ices_cod[[name]]
  if (!all(c("Year", "SSB", "low_SSB", "high_SSB", "stockSizeUnits") %in% names(z)))
    stop("Adult SSB: missing ICES fields for ", name)
  units <- unique(stats::na.omit(as.character(z$stockSizeUnits)))
  if (length(units) != 1L || tolower(units) != "tonnes")
    stop("Adult SSB: units are not tonnes for ", name)
  tibble::tibble(
    Year = as.integer(z$Year), Component = name,
    SSB_t = as.numeric(z$SSB),
    Low_t = as.numeric(z$low_SSB),
    High_t = as.numeric(z$high_SSB),
    Assessment_key = as.character(assessment_keys[[name]])
  ) |>
    dplyr::filter(Year %in% 2020:2025)
}

adult_components <- dplyr::bind_rows(lapply(
  c("Viking", "Northwestern", "Southern"), read_adult_component
))
if (anyDuplicated(adult_components[c("Year", "Component")]))
  stop("Adult SSB: duplicate year/component")

# The ICES component intervals are NOT added to form a 95% CI.
# Instead derive approximate component SD from interval width / 3.29,
# consistent with the pre-existing extract_stock() convention, and
# combine SDs in quadrature. Independence is an explicit assumption.
# This is a working uncertainty scenario, not a validated joint CI.
adult_new <- adult_components |>
  dplyr::group_by(Year) |>
  dplyr::summarise(
    N_components = dplyr::n(),
    All_valid = N_components == 3L &&
      all(is.finite(SSB_t) & SSB_t > 0 &
          is.finite(Low_t) & is.finite(High_t) & High_t > Low_t),
    SSB_t = if (All_valid) sum(SSB_t) else NA_real_,
    SSB_sd_t = if (All_valid)
      sqrt(sum(((High_t - Low_t) / 3.29)^2)) else NA_real_,
    .groups = "drop"
  )

adult_2325 <- adult_new |> dplyr::filter(Year %in% 2023:2025)
if (nrow(adult_2325) != 3L || !all(adult_2325$All_valid))
  stop("Adult SSB: incomplete new assessment data for 2023--2025")
if (nrow(adult_old |> dplyr::filter(Year %in% 2012:2022,
                                  is.finite(Combined_biomass_t))) != 11L)
  stop("Adult SSB: original historical series is incomplete for 2012--2022")

# Overlap sensitivity (2020--2022): compare both spatially scaled
# approaches in the SAME calendar year, rather than comparing adjacent
# years. The difference quantifies method sensitivity, not ecological change.
adult_overlap <- adult_old |>
  dplyr::filter(Year %in% 2020:2022) |>
  dplyr::left_join(adult_new, by = "Year") |>
  dplyr::transmute(
    Year, Old_assessment_key = as.character(Assessment_key),
    Old_SSB_t = Combined_biomass_t,
    New_SSB_t = SSB_t,
    Old_factor = adult_old_factor,
    New_factor = adult_new_factor,
    Old_regional_t = Old_SSB_t * Old_factor,
    New_regional_t = New_SSB_t * New_factor,
    Difference_pct = 100 * (New_regional_t / Old_regional_t - 1)
  )
if (nrow(adult_overlap) != 3L ||
    any(!is.finite(adult_overlap$Difference_pct)))
  stop("Adult SSB: overlap sensitivity check failed")

ices_output_dir <- file.path(final_dir, "Tables", "ICES")
readr::write_csv(adult_components,
  file.path(ices_output_dir, "adult_cod_ssb_new_components.csv"))
readr::write_csv(adult_overlap,
  file.path(ices_output_dir, "adult_cod_ssb_overlap_sensitivity.csv"))

# Replace only 2023--2025 adult stock rows. Historical adult and ALL
# juvenile stock rows are unchanged. Regional_ratio is applied once in
# draw_stock_year(), as in the original Monte Carlo workflow.
adult_replacement <- adult_2325 |>
  dplyr::transmute(
    Year, Stock = "Cod_Adults",
    Biomass_metric = "SSB_three_substocks_scenario",
    Assessment_key = "19661+19662+19663",
    Combined_biomass_t = SSB_t,
    Combined_biomass_sd_t = SSB_sd_t,
    Regional_ratio = adult_new_factor,
    Regional_ratio_cv = 0.25
  )
stock_components <- stock_components |>
  dplyr::filter(!(Stock == "Cod_Adults" & Year %in% 2023:2025)) |>
  dplyr::bind_rows(adult_replacement)

# Plot: one separate sensitivity figure, without modifying any of the
# seven original plots or their styling. Display overlap and the actual
# spliced time series; do not imply that overlap is a continuous series.
adult_sensitivity_plot_data <- adult_overlap |>
  dplyr::select(Year, Old_regional_t, New_regional_t) |>
  tidyr::pivot_longer(-Year, names_to = "Method", values_to = "SSB_t") |>
  dplyr::mutate(Method = dplyr::recode(Method,
    Old_regional_t = "Historical SSB x 0.06",
    New_regional_t = "Three substocks x 0.04449"))
adult_sensitivity_plot <- ggplot2::ggplot(
    adult_sensitivity_plot_data,
    ggplot2::aes(x = Year, y = SSB_t, color = Method, group = Method)
  ) +
  ggplot2::geom_line(linewidth = 1.1) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_x_continuous(breaks = 2020:2022) +
  ggplot2::labs(
    title = "Adult cod SSB: assessment-transition sensitivity",
    subtitle = "Overlapping years 2020–2022; area-scaled scenario estimates",
    x = "Year", y = "Estimated regional SSB (tonnes)", color = "Method"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(legend.position = "bottom")

ggplot2::ggsave(
  file.path(final_dir, "Figures", "Model_comparison",
            "08_Adult_cod_SSB_overlap_sensitivity.png"),
  adult_sensitivity_plot, width = 10, height = 6, dpi = 300
)
cat("\n========== ADULT COD SSB TRANSITION QC ==========\n")
print(adult_overlap, n = Inf)
print(adult_2325, n = Inf)

# ================================================================
# 5. WRASSE BIOMASS — OPTION A (IMR BEACH SEINE CPUE ANCHORED)
# ================================================================
# Choose "beach_seine_point" (interannual beach seine points, Option A2)
# or "beach_seine_smooth" (GAM/LOESS trend from IMR Fig. 1, Option A1)
wrasse_cpue_mode <- "beach_seine_point"

wrasse_file <- file.path(
  input_dir,
  "område og art_data.csv"
)

if (!file.exists(wrasse_file)) {
  stop("Missing wrasse file: ", wrasse_file)
}

wrasse_raw <- readr::read_delim(
  wrasse_file,
  delim = "\t",
  locale = readr::locale(
    encoding = "UTF-16",
    decimal_mark = ","
  ),
  show_col_types = FALSE
)

names(wrasse_raw) <- trimws(
  names(wrasse_raw)
)

wrasse_required <- c(
  "Art - gruppe",
  "Områdegruppering",
  "Fangstår",
  "Measure Values"
)

missing_wrasse <- setdiff(
  wrasse_required,
  names(wrasse_raw)
)

if (length(missing_wrasse)) {
  stop(
    "Unexpected wrasse file columns: ",
    paste(missing_wrasse, collapse = ", "),
    ". Inspect names(wrasse_raw)."
  )
}

wrasse_catches <- wrasse_raw |>
  dplyr::filter(
    stringr::str_detect(
      `Art - gruppe`,
      stringr::regex("leppefisk", ignore_case = TRUE)
    ),
    stringr::str_detect(
      Områdegruppering,
      stringr::regex("Skagerrak", ignore_case = TRUE)
    )
  ) |>
  dplyr::group_by(
    Year = as.integer(Fangstår)
  ) |>
  dplyr::summarise(
    Catch_t = sum(
      as.numeric(`Measure Values`),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::filter(Year %in% years)

# Digitized IMR Flødevigen Beach Seine Survey (Strandnotundersøkelsen)
# from Halvorsen et al. (2026), Rapport fra havforskningen 2026-21, Figure 1.
# Values represent mean combined Labridae (berggylt + bergnebb + grønngylt)
# individuals per beach-seine haul across the 5 Skagerrak subregions
# (Søgne-Grimstad, Arendal-Risør, Kragerø-Larvik, Ytre Oslofjord, Indre Oslofjord).
# 2025 is carried forward from 2024 (rule = 2), matching cormorant census handling.
wrasse_cpue_table <- tibble::tribble(
  ~Year, ~Berggylt_Smooth, ~Bergnebb_Smooth, ~Gronngylt_Smooth, ~CPUE_Smooth_Mean, ~CPUE_Smooth_SE, ~CPUE_Point_Mean, ~CPUE_Point_SE,
  2012L, 1.00,             13.44,            13.18,             27.61,             1.05,            19.97,            5.00,
  2013L, 1.02,             13.44,            13.34,             27.80,             1.02,            27.19,            6.27,
  2014L, 1.04,             13.38,            13.60,             28.02,             1.01,            37.24,            9.71,
  2015L, 1.06,             13.32,            13.86,             28.23,             1.00,            24.44,            5.08,
  2016L, 1.08,             13.23,            14.10,             28.41,             1.00,            37.19,            6.42,
  2017L, 1.11,             13.10,            14.30,             28.51,             1.02,            27.23,            2.89,
  2018L, 1.13,             12.94,            14.58,             28.65,             1.05,            30.80,            5.07,
  2019L, 1.14,             12.74,            14.81,             28.69,             1.10,            20.80,            3.90,
  2020L, 1.15,             12.54,            14.97,             28.67,             1.16,            33.18,            8.29,
  2021L, 1.17,             12.42,            15.08,             28.67,             1.24,            30.57,            6.50,
  2022L, 1.20,             12.39,            15.12,             28.71,             1.36,            26.33,            3.85,
  2023L, 1.23,             12.51,            15.08,             28.82,             1.51,            28.11,            5.00,
  2024L, 1.27,             12.81,            14.96,             29.03,             1.70,            25.89,            5.65,
  2025L, 1.27,             12.81,            14.96,             29.03,             1.70,            25.89,            5.65
)

# Save digitized CPUE table to Data/Inputs and Monte_Carlo output directory for provenance
readr::write_csv(
  wrasse_cpue_table,
  file.path(input_dir, "wrasse_cpue_skagerrak.csv")
)
readr::write_csv(
  wrasse_cpue_table,
  file.path(out_dir, "IMR_wrasse_beach_seine_CPUE_2012_2025.csv")
)

# Anchor unscaled biomass to the peak active commercial fishery window (2015-2018)
wrasse_anchor_years <- 2015:2018

wrasse_anchor_unscaled_t <- wrasse_catches |>
  dplyr::filter(Year %in% wrasse_anchor_years) |>
  dplyr::summarise(
    Mean_unscaled_t = mean(Catch_t / (1 - exp(-0.2)), na.rm = TRUE)
  ) |>
  dplyr::pull(Mean_unscaled_t)

wrasse_stock <- wrasse_cpue_table |>
  dplyr::filter(Year %in% years) |>
  dplyr::left_join(wrasse_catches, by = "Year") |>
  dplyr::mutate(
    Active_CPUE = if (wrasse_cpue_mode == "beach_seine_smooth") {
      CPUE_Smooth_Mean
    } else {
      CPUE_Point_Mean
    },
    Active_SE = if (wrasse_cpue_mode == "beach_seine_smooth") {
      CPUE_Smooth_SE
    } else {
      CPUE_Point_SE
    },
    Ref_CPUE_2015_2018 = mean(
      Active_CPUE[Year %in% wrasse_anchor_years],
      na.rm = TRUE
    ),
    Relative_CPUE_index = Active_CPUE / Ref_CPUE_2015_2018,
    Relative_CPUE_cv = pmax(Active_SE / Active_CPUE, 0.20),
    Stock = "Wrasse",
    Biomass_metric = paste0("IMR_Strandnot_", wrasse_cpue_mode, "_anchored_2015_2018"),
    Assessment_key = "IMR_2026_21_Fig1",
    Combined_biomass_t = wrasse_anchor_unscaled_t * Relative_CPUE_index,
    Combined_biomass_sd_t = Combined_biomass_t * Relative_CPUE_cv,
    Regional_ratio = 20.0, # ONE application only: original scenario multiplier
    Regional_ratio_cv = 0.30
  ) |>
  dplyr::select(
    Year,
    Stock,
    Biomass_metric,
    Assessment_key,
    Combined_biomass_t,
    Combined_biomass_sd_t,
    Regional_ratio,
    Regional_ratio_cv
  )

stock_components <- dplyr::bind_rows(
  stock_components,
  wrasse_stock
)


# ================================================================
# 6. PREDATOR DEMOGRAPHICS
# ================================================================

seal_K <- 13965
seal_N0 <- 5930
seal_mu <- 0.28
seal_skagerrak_ratio <- 1.0

# Corrected calendar alignment must be checked against the source census.
# The original annual script hardcoded 2012,2016,2020,2024.
# Set CORM_CENSUS_YEARS to VERIFIED census calendar years (comma separated).
corm_year_env <- Sys.getenv("CORM_CENSUS_YEARS", unset="")
# Preserve original calendar mapping as a provisional scenario unless an
# independently verified census-year mapping is supplied by the user.
if (!nzchar(corm_year_env)) {
  warning("Cormorant census calendar years are not independently verified; using original script mapping 2012, 2016, 2020, 2024 as a scenario.")
  corm_census_years <- c(2012L, 2016L, 2020L, 2024L)
} else {
  corm_census_years <- as.integer(strsplit(corm_year_env, ",")[[1]])
}
if (length(corm_census_years)!=4L || anyNA(corm_census_years) || anyDuplicated(corm_census_years) || is.unsorted(corm_census_years)) stop("Invalid CORM_CENSUS_YEARS")

corm_census_counts <- c(
  28500, 33750, 44137, 42875
)

predator_demographics <- tibble::tibble(
  Year = years
) |>
  dplyr::mutate(
    
    Year_c = Year - 2012,
    
    Seal_N_mean =
      (
        seal_K /
          (
            1 +
              ((seal_K - seal_N0) / seal_N0) *
              exp(-seal_mu * Year_c)
          )
      ) * seal_skagerrak_ratio,
    
    Seal_N_sd = Seal_N_mean * 0.15,
    
    Corm_N_mean = stats::approx(
      x = corm_census_years,
      y = corm_census_counts,
      xout = Year,
      rule = 2
    )$y,
    
    Corm_N_sd = Corm_N_mean * 0.15
    
  )

# 2025 cormorants: approx(rule=2) carries the 2024 anchor forward.
qc_predator_2025 <- predator_demographics |>
  dplyr::filter(Year == 2025L)
if (nrow(qc_predator_2025) != 1L ||
    !isTRUE(all.equal(qc_predator_2025$Corm_N_mean,
      stats::approx(corm_census_years, corm_census_counts,
                    xout = 2025, rule = 2)$y)) ||
    any(!is.finite(unlist(qc_predator_2025)))) {
  stop("v1.4 QC: 2025 predator population mismatch")
}
readr::write_csv(qc_predator_2025,
  file.path(out_dir, "predator_demographics_2025_QC.csv"))

readr::write_csv(
  predator_demographics,
  file.path(out_dir, "predator_demographics.csv")
)


# ================================================================
# 7. COMMON MONTE CARLO RANDOM DRAWS
# ================================================================
#
# Positive lognormal draws with specified arithmetic mean
# and standard deviation.
#
# Biological and stock draws are generated once per year
# and reused across all three diet models.
# ================================================================

rlnorm_mean_sd <- function(n, mean, sd) {
  
  if (
    !is.finite(mean) ||
    !is.finite(sd) ||
    mean <= 0 ||
    sd < 0
  ) {
    return(rep(NA_real_, n))
  }
  
  sdlog <- sqrt(
    log1p((sd / mean)^2)
  )
  
  meanlog <- log(mean) -
    0.5 * sdlog^2
  
  stats::rlnorm(
    n,
    meanlog = meanlog,
    sdlog = sdlog
  )
}

draw_stock_year <- function(year, n) {
  
  rows <- stock_components |>
    dplyr::filter(Year == year)
  
  stock_names <- c(
    "Sprat",
    "NSAS_Herring",
    "WBSS_Herring",
    "Cod_Juveniles",
    "Cod_Adults",
    "Wrasse"
  )
  
  stock_draws <- setNames(
    lapply(stock_names, function(x) {
      rep(NA_real_, n)
    }),
    stock_names
  )
  
  for (stock_name in stock_names) {
    
    row <- rows |>
      dplyr::filter(Stock == stock_name)
    
    if (nrow(row) == 0) {
      next
    }
    
    biomass <- rlnorm_mean_sd(
      n,
      row$Combined_biomass_t[1],
      row$Combined_biomass_sd_t[1]
    )
    
    ratio <- rlnorm_mean_sd(
      n,
      row$Regional_ratio[1],
      row$Regional_ratio[1] *
        row$Regional_ratio_cv[1]
    )
    
    stock_draws[[stock_name]] <- biomass * ratio
  }
  
  # All three components must be available.
  # Missing components propagate NA, never zero.
  
  clupeid <- Reduce(
    `+`,
    stock_draws[c(
      "Sprat",
      "NSAS_Herring",
      "WBSS_Herring"
    )]
  )
  
  list(
    Cod_Juveniles = stock_draws[["Cod_Juveniles"]],
    Cod_Adults = stock_draws[["Cod_Adults"]],
    Clupeid = clupeid,
    Wrasse = stock_draws[["Wrasse"]]
  )
}

draw_common_year <- function(year, n) {
  
  demo <- predator_demographics |>
    dplyr::filter(Year == year)
  
  if (nrow(demo) != 1) {
    stop("Invalid predator demographics year.")
  }
  
  seal_N <- rlnorm_mean_sd(
    n,
    demo$Seal_N_mean,
    demo$Seal_N_sd
  )
  
  corm_N <- rlnorm_mean_sd(
    n,
    demo$Corm_N_mean,
    demo$Corm_N_sd
  )
  
  seal_D <- rlnorm_mean_sd(
    n, 4.0, 0.5
  )
  
  corm_D <- rlnorm_mean_sd(
    n, 0.4, 0.05
  )
  
  seal_T <- stats::runif(
    n, 240, 290
  )
  
  corm_T <- stats::runif(
    n, 200, 265
  )
  
  # Consumption in kg per year.
  
  seal_total_kg <- seal_N * seal_D * seal_T
  corm_total_kg <- corm_N * corm_D * corm_T
  
  seal_juvenile_fraction <- pmin(
    pmax(stats::rnorm(n, 0.65, 0.10), 0),
    1
  )
  
  corm_juvenile_fraction <- pmin(
    pmax(stats::rnorm(n, 0.95, 0.05), 0),
    1
  )
  
  list(
    Seal_total_t = seal_total_kg / 1000,
    Cormorant_total_t = corm_total_kg / 1000,
    Seal_juvenile_fraction = seal_juvenile_fraction,
    Cormorant_juvenile_fraction =
      corm_juvenile_fraction,
    Stock = draw_stock_year(year, n)
  )
}

set.seed(seed)

common_draws <- setNames(
  lapply(years, function(year) {
    
    message("Generating common draws: ", year)
    
    draw_common_year(year, n_iter)
    
  }),
  as.character(years)
)


# ================================================================
# 8. MONTE CARLO SIMULATION
# ================================================================

run_model_year <- function(
    model_name,
    year,
    common,
    n
) {
  
  seal_diet <- draw_diet_composition(
    "Seal",
    model_name,
    year,
    n
  )
  
  corm_diet <- draw_diet_composition(
    "Cormorant",
    model_name,
    year,
    n
  )
  
  seal_total <- common$Seal_total_t
  corm_total <- common$Cormorant_total_t
  
  # Draw-level consumption.
  
  seal_cod <- seal_total * seal_diet[, "Cod"]
  corm_cod <- corm_total * corm_diet[, "Cod"]
  
  seal_clupeid <- seal_total * seal_diet[, "Clupeid"]
  corm_clupeid <- corm_total * corm_diet[, "Clupeid"]
  
  seal_wrasse <- seal_total * seal_diet[, "Wrasse"]
  corm_wrasse <- corm_total * corm_diet[, "Wrasse"]
  
  seal_cod_juv <- seal_cod *
    common$Seal_juvenile_fraction
  
  seal_cod_adult <- seal_cod *
    (1 - common$Seal_juvenile_fraction)
  
  corm_cod_juv <- corm_cod *
    common$Cormorant_juvenile_fraction
  
  corm_cod_adult <- corm_cod *
    (1 - common$Cormorant_juvenile_fraction)
  
  # IMPORTANT:
  # Sum predator consumption within each draw,
  # not after calculating separate medians.
  
  tibble::tibble(
    
    Model = model_name,
    Year = year,
    Iteration = seq_len(n),
    
    Seal_total_t = seal_total,
    Cormorant_total_t = corm_total,
    
    Seal_Cod_Juveniles_t = seal_cod_juv,
    Cormorant_Cod_Juveniles_t = corm_cod_juv,
    
    Seal_Cod_Adults_t = seal_cod_adult,
    Cormorant_Cod_Adults_t = corm_cod_adult,
    
    Seal_Clupeid_t = seal_clupeid,
    Cormorant_Clupeid_t = corm_clupeid,
    
    Seal_Wrasse_t = seal_wrasse,
    Cormorant_Wrasse_t = corm_wrasse,
    
    Cod_Juveniles_stock_t =
      common$Stock$Cod_Juveniles,
    
    Cod_Adults_stock_t =
      common$Stock$Cod_Adults,
    
    Clupeid_stock_t =
      common$Stock$Clupeid,
    
    Wrasse_stock_t =
      common$Stock$Wrasse
    
  ) |>
    dplyr::mutate(
      
      Total_Cod_Juveniles_t =
        Seal_Cod_Juveniles_t +
        Cormorant_Cod_Juveniles_t,
      
      Total_Cod_Adults_t =
        Seal_Cod_Adults_t +
        Cormorant_Cod_Adults_t,
      
      Total_Clupeid_t =
        Seal_Clupeid_t +
        Cormorant_Clupeid_t,
      
      Total_Wrasse_t =
        Seal_Wrasse_t +
        Cormorant_Wrasse_t,
      
      Total_Cod_t =
        Total_Cod_Juveniles_t +
        Total_Cod_Adults_t,
      
      Cod_Juveniles_ratio =
        Total_Cod_Juveniles_t /
        Cod_Juveniles_stock_t,
      
      Cod_Adults_ratio =
        Total_Cod_Adults_t /
        Cod_Adults_stock_t,
      
      Clupeid_ratio =
        Total_Clupeid_t /
        Clupeid_stock_t,
      
      Wrasse_ratio =
        Total_Wrasse_t /
        Wrasse_stock_t
      
    )
}

all_mc <- list()

for (model_name in diet_models) {
  
  message("\nRunning Monte Carlo: ", model_name)
  
  model_results <- lapply(years, function(year) {
    
    message("  Year: ", year)
    
    run_model_year(
      model_name = model_name,
      year = year,
      common = common_draws[[as.character(year)]],
      n = n_iter
    )
  })
  
  all_mc[[model_name]] <- dplyr::bind_rows(
    model_results
  )
  
  saveRDS(
    all_mc[[model_name]],
    file.path(
      out_dir,
      paste0("MC_draws_", model_name, ".rds")
    )
  )
}

mc_draws <- dplyr::bind_rows(all_mc)


# ================================================================
# 9. SUMMARY STATISTICS
# ================================================================

summary_variables <- c(
  "Seal_total_t",
  "Cormorant_total_t",
  "Total_Cod_Juveniles_t",
  "Total_Cod_Adults_t",
  "Total_Clupeid_t",
  "Total_Wrasse_t",
  "Total_Cod_t",
  "Cod_Juveniles_ratio",
  "Cod_Adults_ratio",
  "Clupeid_ratio",
  "Wrasse_ratio"
)

summarise_mc <- function(x) {
  
  valid <- x[is.finite(x)]
  
  if (!length(valid)) {
    
    return(tibble::tibble(
      Median = NA_real_,
      Mean = NA_real_,
      LCI = NA_real_,
      UCI = NA_real_,
      N_valid = 0L
    ))
  }
  
  tibble::tibble(
    Median = stats::median(valid),
    Mean = mean(valid),
    LCI = as.numeric(
      stats::quantile(valid, 0.025)
    ),
    UCI = as.numeric(
      stats::quantile(valid, 0.975)
    ),
    N_valid = length(valid)
  )
}

mc_summary <- mc_draws |>
  dplyr::select(
    Model,
    Year,
    dplyr::all_of(summary_variables)
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::all_of(summary_variables),
    names_to = "Metric",
    values_to = "Value"
  ) |>
  dplyr::group_by(
    Model,
    Year,
    Metric
  ) |>
  dplyr::group_modify(
    ~summarise_mc(.x$Value)
  ) |>
  dplyr::ungroup()

readr::write_csv(
  mc_summary,
  file.path(
    out_dir,
    "MC_summary_all_models.csv"
  )
)

for (model_name in diet_models) {
  
  readr::write_csv(
    mc_summary |>
      dplyr::filter(Model == model_name),
    file.path(
      out_dir,
      paste0(
        "MC_summary_",
        model_name,
        ".csv"
      )
    )
  )
}


# Original plot figures are generated in the final section below.

# ================================================================
# 12. VALIDATION TABLES
# ================================================================

stock_coverage_qc <- tidyr::expand_grid(
  Stock = unique(stock_components$Stock),
  Year = years
) |>
  dplyr::left_join(
    stock_components,
    by = c("Stock", "Year")
  ) |>
  dplyr::mutate(
    Biomass_available =
      is.finite(Combined_biomass_t) &
      Combined_biomass_t > 0
  )

readr::write_csv(
  stock_coverage_qc,
  file.path(
    out_dir,
    "stock_coverage_QC.csv"
  )
)

mc_validity_qc <- mc_summary |>
  dplyr::mutate(
    Expected_iterations = n_iter,
    Missing_draws = n_iter - N_valid
  ) |>
  dplyr::filter(
    Missing_draws > 0 |
      !is.finite(Median)
  )

readr::write_csv(
  mc_validity_qc,
  file.path(
    out_dir,
    "MC_missing_draws_QC.csv"
  )
)

model_draw_counts <- mc_draws |>
  dplyr::count(
    Model,
    Year,
    name = "N_draws"
  )

readr::write_csv(
  model_draw_counts,
  file.path(
    out_dir,
    "MC_iteration_QC.csv"
  )
)

if (any(model_draw_counts$N_draws != n_iter)) {
  stop("Incorrect Monte Carlo iteration count.")
}


# ================================================================
# 13. FINAL REPORT
# ================================================================

cat("\n========================================\n")
cat("MONTE CARLO INTEGRATION — VERSION 0.5\n")
cat("========================================\n")

cat(
  "Diet models:",
  paste(diet_models, collapse = ", "),
  "\n"
)

cat(
  "Years:",
  min(years),
  "to",
  max(years),
  "\n"
)

cat(
  "Iterations per year/model:",
  n_iter,
  "\n"
)

cat(
  "Total simulation rows:",
  nrow(mc_draws),
  "\n"
)

cat(
  "Missing summary combinations:",
  nrow(mc_validity_qc),
  "\n"
)

cat("\nICES assessments:\n")
print(ices_provenance)

cat("\nStock coverage issues:\n")

print(
  stock_coverage_qc |>
    dplyr::filter(!Biomass_available),
  n = Inf
)

cat("\nMonte Carlo missing-draw issues:\n")

print(mc_validity_qc, n = Inf)

cat("\nOutput tables:\n", out_dir, "\n")
cat("Output figures:\n", fig_dir, "\n")

cat("\nIMPORTANT INTERPRETATION NOTES:\n")

cat(
  "- Diet RRA is not calibrated prey biomass.\n",
  "- Spatial ratios are model assumptions.\n",
  "- Cod adult biomass uses SSB, not total biomass.\n",
  "- Juvenile cod biomass is recruitment-derived.\n",
  "- Wrasse biomass is a catch-derived proxy.\n",
  "- Pre-observation cormorant diet is held at its first observed year.\n",
  "- Missing ICES years are not imputed.\n",
  "- No causal predation effect is estimated.\n",
  sep = ""
)

cat("\nDevelopment run completed.\n")


# Original plot figures are generated in the final section below.
# Source-verified stock calculations are identical across diet models;
# only diet model and its uncertainty simulation differ.
readr::write_csv(data.frame(Model=diet_models,Draws=n_iter,
  Cod_stock_scenario='Age1+Age2',Cod_survival=survival,
  ICES_old_key=17652L,ICES_new_keys='19661+19662+19663',
  Spatial_fraction_2023_2025=f_new,
  Adult_SSB_old_assessment='Original cod.27.47d20 (2012-2022)',
  Adult_SSB_new_assessments='19661+19662+19663 (2023-2025)',
  Adult_SSB_old_factor=adult_old_factor,
  Adult_SSB_new_factor=adult_new_factor,
  Adult_SSB_uncertainty='Component SDs from interval width / 3.29; independent quadrature'),
  file.path(final_dir,'Documentation','model_run_manifest.csv'))
message('FINAL RESULT folder: ',final_dir)

# IMPORTANT: MC ratio = annual consumption / standing biomass,
# not a direct estimate of annual predation mortality. RRA is not
# calibrated to biomass fractions. All three models share the same
# ICES stock scenario and demographic assumptions.


# ============================================================================
# ORIGINAL FIGURES — exact original ggplot layers/themes/labels/colours
# Data are mapped from the revised LM / quasibinomial / hurdle MC outputs.
# The original plot code is preserved below; no new plot styles are added.
# ============================================================================
if (!requireNamespace("patchwork", quietly=TRUE)) stop("Install patchwork: install.packages('patchwork')")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# The same original plots are produced separately for each diet model.
# NEW v1.3: sensitivity figure 08 is generated once in Section 4C above.
# All MC summaries use medians across iterations (not medians summed afterward).
for (model_name in diet_models) {
  original_dir <- file.path(final_dir, "Figures", model_name)
  dir.create(original_dir, recursive=TRUE, showWarnings=FALSE)
  d <- mc_draws |> dplyr::filter(Model == model_name)
  agg <- d |> dplyr::group_by(Year) |> dplyr::summarise(
    Seal_CodJuv_Impact_Median=median(Seal_Cod_Juveniles_t / pmax(Cod_Juveniles_stock_t, .1),na.rm=TRUE),
    Corm_CodJuv_Impact_Median=median(Cormorant_Cod_Juveniles_t / pmax(Cod_Juveniles_stock_t, .1),na.rm=TRUE),
    Seal_CodAdult_Impact_Median=if(all(is.na(Cod_Adults_stock_t))) NA_real_ else median(Seal_Cod_Adults_t / pmax(Cod_Adults_stock_t,.1),na.rm=TRUE),
    Corm_CodAdult_Impact_Median=if(all(is.na(Cod_Adults_stock_t))) NA_real_ else median(Cormorant_Cod_Adults_t / pmax(Cod_Adults_stock_t,.1),na.rm=TRUE),
    Seal_Clupeid_Impact_Median=median(Seal_Clupeid_t / pmax(Clupeid_stock_t,.1),na.rm=TRUE),
    Corm_Clupeid_Impact_Median=median(Cormorant_Clupeid_t / pmax(Clupeid_stock_t,.1),na.rm=TRUE),
    Seal_Wrasse_Impact_Median=median(Seal_Wrasse_t / pmax(Wrasse_stock_t,.1),na.rm=TRUE),
    Corm_Wrasse_Impact_Median=median(Cormorant_Wrasse_t / pmax(Wrasse_stock_t,.1),na.rm=TRUE),
    Seal_CodJuv_C_Tonnes=median(Seal_Cod_Juveniles_t,na.rm=TRUE),
    Corm_CodJuv_C_Tonnes=median(Cormorant_Cod_Juveniles_t,na.rm=TRUE),
    Seal_CodAdult_C_Tonnes=median(Seal_Cod_Adults_t,na.rm=TRUE),
    Corm_CodAdult_C_Tonnes=median(Cormorant_Cod_Adults_t,na.rm=TRUE),
    Seal_Clupeid_C_Tonnes=median(Seal_Clupeid_t,na.rm=TRUE),
    Corm_Clupeid_C_Tonnes=median(Cormorant_Clupeid_t,na.rm=TRUE),
    Seal_Wrasse_C_Tonnes=median(Seal_Wrasse_t,na.rm=TRUE),
    Corm_Wrasse_C_Tonnes=median(Cormorant_Wrasse_t,na.rm=TRUE),
    Seal_Total_C_Tonnes=median(Seal_total_t,na.rm=TRUE),
    Corm_Total_C_Tonnes=median(Cormorant_total_t,na.rm=TRUE),
    .groups="drop")
  mc_results <- agg
  impact_plot_data <- mc_results %>%
    dplyr::select(Year, ends_with("_Impact_Median")) %>%
    pivot_longer(cols=-Year,names_to=c("Predator","Prey"),names_pattern="(.*)_(.*)_Impact_Median") %>%
    mutate(Predator=if_else(Predator=="Seal","Harbour Seals","Great Cormorants"),
           Prey=case_when(Prey=="CodJuv"~"Cod Juveniles",Prey=="CodAdult"~"Cod Adults",Prey=="Clupeid"~"Clupeidae",Prey=="Wrasse"~"Wrasses"),
           Prey=factor(Prey,levels=c("Cod Juveniles","Cod Adults","Clupeidae","Wrasses")))
  consumption_plot_data <- mc_results %>%
    dplyr::select(Year, matches("^(Seal|Corm)_(CodJuv|CodAdult|Clupeid|Wrasse)_C_Tonnes$")) %>%
    pivot_longer(cols=-Year,names_to=c("Predator","Prey"),names_pattern="(.*)_(.*)_C_Tonnes") %>%
    mutate(Predator=if_else(Predator=="Seal","Harbour Seals","Great Cormorants"),
           Prey=case_when(Prey=="CodJuv"~"Cod Juveniles",Prey=="CodAdult"~"Cod Adults",Prey=="Clupeid"~"Clupeidae",Prey=="Wrasse"~"Wrasses"),
           Prey=factor(Prey,levels=c("Cod Juveniles","Cod Adults","Clupeidae","Wrasses")))
p1 <- ggplot(impact_plot_data, aes(x = Year, y = value, color = Predator)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Prey, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("Harbour Seals" = "#8c510a", "Great Cormorants" = "#01665e")) +
  scale_x_continuous(breaks = seq(2012, 2025, by = 3)) +
  theme_minimal() +
  labs(
    title = "Top-Down Impact Ratio by Predator and Prey Stock",
    y = "Impact Ratio (Consumption / Biomass)",
    x = "Year"
  ) +
  theme(legend.position = "bottom", legend.title = element_blank(), strip.text = element_text(face = "bold"))

p2 <- ggplot(consumption_plot_data, aes(x = Year, y = value, fill = Prey)) +
  geom_area(alpha = 0.8, color = "white", linewidth = 0.2) +
  facet_wrap(~ Predator, scales = "free_y") +
  # Use light blue for juveniles, dark blue for adults
  scale_fill_manual(values = c("Cod Juveniles" = "#74add1", "Cod Adults" = "#313695", 
                               "Clupeidae" = "darkorange", "Wrasses" = "seagreen")) +
  scale_x_continuous(breaks = seq(2012, 2025, by = 3)) +
  theme_minimal() +
  labs(
    title = "Total Prey Consumption by Species and Life Stage",
    y = "Consumption (Metric Tonnes)",
    x = "Year"
  ) +
  theme(legend.position = "bottom", legend.title = element_blank(), strip.text = element_text(face = "bold"))


  all_prey_preds <- diet_predictions %>% filter(Model==model_name) %>%
    transmute(x=Year-2012,predicted=Predicted_RRA,Predator=if_else(Predator=="Seal","Harbour Seals","Great Cormorants"),Prey=Prey)
  # Original RRA age partitions are assumptions, not age-specific DNA observations.
  all_prey_preds <- bind_rows(
    all_prey_preds %>% filter(Prey=="Cod",Predator=="Harbour Seals") %>% mutate(predicted=predicted*.65,Prey="Cod Juveniles"),
    all_prey_preds %>% filter(Prey=="Cod",Predator=="Harbour Seals") %>% mutate(predicted=predicted*.35,Prey="Cod Adults"),
    all_prey_preds %>% filter(Prey=="Cod",Predator=="Great Cormorants") %>% mutate(predicted=predicted*.95,Prey="Cod Juveniles"),
    all_prey_preds %>% filter(Prey=="Cod",Predator=="Great Cormorants") %>% mutate(predicted=predicted*.05,Prey="Cod Adults"),
    all_prey_preds %>% filter(Prey!="Cod") %>% mutate(Prey=recode(Prey,Clupeid="Clupeidae",Wrasse="Wrasses"))
  ) %>% mutate(Prey=factor(Prey,levels=c("Cod Juveniles","Cod Adults","Clupeidae","Wrasses")),
              conf.low=NA_real_,conf.high=NA_real_)
  # Uncertainty bands are deliberately NA rather than misusing CI bounds as SE.
  # The original geom_ribbon is retained and will be populated only with valid CIs.

  prey_plot_data <- stock_components %>%
    filter(Year>=2012,Year<=2025) %>%
    mutate(Stock=recode(Stock,Cod_Juveniles="Cod Juveniles",Cod_Adults="Cod Adults",Wrasse="Wrasses",NSAS_Herring="NSAS Herring",WBSS_Herring="WBSS Herring"),
           Skagerrak_Biomass_Median=Combined_biomass_t*Regional_ratio,
           Skagerrak_Biomass_LCI=pmax(Skagerrak_Biomass_Median-1.96*Combined_biomass_sd_t*Regional_ratio,0),
           Skagerrak_Biomass_UCI=Skagerrak_Biomass_Median+1.96*Combined_biomass_sd_t*Regional_ratio)
  predator_plot_data <- predator_demographics %>%
    mutate(Seal_LCI=Seal_N_mean-1.96*Seal_N_sd,Seal_UCI=Seal_N_mean+1.96*Seal_N_sd,
           Corm_LCI=Corm_N_mean-1.96*Corm_N_sd,Corm_UCI=Corm_N_mean+1.96*Corm_N_sd)
diet_trend_plot <- ggplot(all_prey_preds, aes(x = x + 2012, y = predicted, color = Predator, fill = Predator)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  facet_wrap(~ Prey, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("Harbour Seals" = "#8c510a", "Great Cormorants" = "#01665e")) +
  scale_fill_manual(values = c("Harbour Seals" = "#8c510a", "Great Cormorants" = "#01665e")) +
  scale_x_continuous(breaks = seq(2012, 2025, by = 3)) + 
  theme_minimal() +
  labs(
    title = "Season-Adjusted Dietary RRA Trends Over Time",
    x = "Year",
    y = "Predicted Relative Read Abundance (RRA)"
  ) +
  theme(legend.position = "bottom", legend.title = element_blank(), strip.text = element_text(face = "bold"))

prey_dynamics_plot <- ggplot(prey_plot_data, aes(x = Year, y = Skagerrak_Biomass_Median, color = Stock, fill = Stock)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = Skagerrak_Biomass_LCI, ymax = Skagerrak_Biomass_UCI), alpha = 0.2, color = NA) +
  facet_wrap(~ Stock, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = seq(2012, 2025, by = 3)) +
  scale_color_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option = "turbo") +
  theme_minimal() +
  labs(
    title = "Estimated Prey Population Dynamics in the Skagerrak",
    x = "Year", 
    y = "Standing Stock Biomass (Metric Tonnes)"
  ) +
  theme(
    legend.position = "none", 
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )

seal_pop_plot <- ggplot(predator_plot_data, aes(x = Year, y = Seal_N_mean)) +
  geom_line(color = "#8c510a", linewidth = 1.2) +
  geom_ribbon(aes(ymin = Seal_LCI, ymax = Seal_UCI), fill = "#8c510a", alpha = 0.2) +
  scale_x_continuous(breaks = seq(2012, 2025, by = 3)) +
  theme_minimal() +
  labs(
    title = "Harbour Seal Population Dynamics",
    x = "Year", 
    y = "Number of Individuals"
  )

cormorant_pop_plot <- ggplot(predator_plot_data, aes(x = Year, y = Corm_N_mean)) +
  geom_line(color = "#01665e", linewidth = 1.2) +
  geom_ribbon(aes(ymin = Corm_LCI, ymax = Corm_UCI), fill = "#01665e", alpha = 0.2) +
  scale_x_continuous(breaks = seq(2012, 2025, by = 3)) +
  theme_minimal() +
  labs(
    title = "Great Cormorant Population Dynamics",
    x = "Year", 
    y = "Number of Individuals"
  )


  # Fisheries extraction comparison requires catch series. Use original ICES raw
  # assessment CSVs when they include catches/landings; never replace NA with zero.
  get_annual_catch_safe <- function(filename,ratio) {
    f <- file.path(final_dir,"Tables","Monte_Carlo",filename)
    if (!file.exists(f)) return(tibble(Year=2012:2025,Catch_Val=NA_real_))
    z <- readr::read_csv(f,show_col_types=FALSE)
    y <- intersect(c("Year","year"),names(z))
    c <- intersect(c("catches","Catches","landings","Landings"),names(z))
    if (!length(y) || !length(c)) return(tibble(Year=2012:2025,Catch_Val=NA_real_))
    z %>% transmute(Year=as.integer(.data[[y[1]]]),Catch_Val=as.numeric(.data[[c[1]]])*ratio) %>% filter(Year %in% 2012:2025) %>% distinct(Year,.keep_all=TRUE)
  }
  annual_comm_cod <- get_annual_catch_safe("ICES_raw_Cod_Adults.csv",.06) %>% rename(Comm_Cod=Catch_Val)
  annual_sprat <- get_annual_catch_safe("ICES_raw_Sprat.csv",.15) %>% rename(Sprat=Catch_Val)
  annual_nsas <- get_annual_catch_safe("ICES_raw_NSAS_Herring.csv",.10) %>% rename(NSAS=Catch_Val)
  annual_wbss <- get_annual_catch_safe("ICES_raw_WBSS_Herring.csv",.40) %>% rename(WBSS=Catch_Val)
  annual_human <- tibble(Year=2012:2025) %>% left_join(annual_comm_cod,by="Year") %>%
    left_join(annual_sprat,by="Year") %>% left_join(annual_nsas,by="Year") %>% left_join(annual_wbss,by="Year") %>%
    mutate(Comm_Pelagic=Sprat+NSAS+WBSS,Rec_Cod=if_else(Year>2021,NA_real_,Comm_Cod*3))
  annual_predator <- mc_results %>% transmute(Year,
    Corm_Total=Corm_Total_C_Tonnes,Seal_Total=Seal_Total_C_Tonnes,
    Corm_Cod=Corm_CodJuv_C_Tonnes+Corm_CodAdult_C_Tonnes,
    Seal_Cod=Seal_CodJuv_C_Tonnes+Seal_CodAdult_C_Tonnes)
  plot_data_time <- annual_human %>% left_join(annual_predator,by="Year") %>%
    pivot_longer(cols=-Year,names_to="Variable",values_to="Biomass_tonnes") %>%
    filter(Variable %in% c("Corm_Total","Corm_Cod","Seal_Total","Seal_Cod","Comm_Pelagic","Comm_Cod","Rec_Cod")) %>%
    mutate(Source=case_when(Variable %in% c("Corm_Total","Corm_Cod")~"Great Cormorant",
      Variable %in% c("Seal_Total","Seal_Cod")~"Harbour Seal",
      Variable %in% c("Comm_Pelagic","Comm_Cod")~"Commercial Fleet",Variable=="Rec_Cod"~"Recreational Fishery"),
      Category=case_when(Variable %in% c("Corm_Total","Seal_Total","Comm_Pelagic")~"Total Community Extraction",TRUE~"Targeted Atlantic Cod Extraction"),
      Type=case_when(Variable %in% c("Corm_Total","Corm_Cod")~"Avian Predation",
        Variable %in% c("Seal_Total","Seal_Cod")~"Mammalian Predation",TRUE~"Human Extraction"),
      Source=factor(Source,levels=c("Great Cormorant","Harbour Seal","Commercial Fleet","Recreational Fishery")))
extraction_time_plot <- ggplot(plot_data_time, aes(x = Year, y = Biomass_tonnes, color = Source, linetype = Type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~ Category, scales = "free_y") + 
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_color_manual(values = c("Great Cormorant" = "#01665e",        
                                "Harbour Seal" = "#8c510a",    
                                "Commercial Fleet" = "grey20",     
                                "Recreational Fishery" = "grey50")) +
  scale_linetype_manual(values = c("Avian Predation" = "solid",
                                   "Mammalian Predation" = "solid",
                                   "Human Extraction" = "dashed")) +
  scale_x_continuous(breaks = seq(2012, 2025, by = 2)) +
  labs(
    title = "Temporal Dynamics of Biomass Extraction",
    y = "Extracted Biomass (Tonnes)",
    x = "Year",
    color = "Extraction Agent",
    linetype = "Extraction Type"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.box = "vertical"
  )

# 6. Display the final plot

  ggplot2::ggsave(file.path(original_dir,"01_Top_down_impact.png"),p1,width=15,height=5,dpi=300)
  ggplot2::ggsave(file.path(original_dir,"02_Total_prey_consumption.png"),p2,width=12,height=6,dpi=300)
  ggplot2::ggsave(file.path(original_dir,"03_Original_A_B_combined.png"),p1 / p2,width=16,height=12,dpi=300)
  ggplot2::ggsave(file.path(original_dir,"04_Dietary_RRA_trends.png"),diet_trend_plot,width=15,height=5,dpi=300)
  ggplot2::ggsave(file.path(original_dir,"05_Prey_population_dynamics.png"),prey_dynamics_plot,width=13,height=10,dpi=300)
  ggplot2::ggsave(file.path(original_dir,"06_Predator_population_dynamics.png"),seal_pop_plot+cormorant_pop_plot,width=12,height=5,dpi=300)
  ggplot2::ggsave(file.path(original_dir,"07_Temporal_biomass_extraction.png"),extraction_time_plot,width=12,height=7,dpi=300)
}
message("Original plotting layouts saved in: ",file.path(final_dir,"Figures"))

# ============================================================================
# ORIGINAL FIGURES — exact original ggplot layers/themes/labels/colours
# Data are mapped from the revised LM / quasibinomial / hurdle MC outputs.
# The original plot code is preserved below; no new plot styles are added.
# ============================================================================
if (!requireNamespace("patchwork", quietly=TRUE)) stop("Install patchwork: install.packages('patchwork')")
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# The same original plots are produced separately for each diet model.
# NEW v1.3: sensitivity figure 08 is generated once in Section 4C above.
# All MC summaries use medians across iterations (not medians summed afterward).
for (model_name in diet_models) {
  original_dir <- file.path(final_dir, "Figures", model_name)
  dir.create(original_dir, recursive=TRUE, showWarnings=FALSE)
  d <- mc_draws |> dplyr::filter(Model == model_name)
  agg <- d |> dplyr::group_by(Year) |> dplyr::summarise(
    Seal_CodJuv_Impact_Median=median(Seal_Cod_Juveniles_t / pmax(Cod_Juveniles_stock_t, .1),na.rm=TRUE),
    Corm_CodJuv_Impact_Median=median(Cormorant_Cod_Juveniles_t / pmax(Cod_Juveniles_stock_t, .1),na.rm=TRUE),
    Seal_CodAdult_Impact_Median=if(all(is.na(Cod_Adults_stock_t))) NA_real_ else median(Seal_Cod_Adults_t / pmax(Cod_Adults_stock_t,.1),na.rm=TRUE),
    Corm_CodAdult_Impact_Median=if(all(is.na(Cod_Adults_stock_t))) NA_real_ else median(Cormorant_Cod_Adults_t / pmax(Cod_Adults_stock_t,.1),na.rm=TRUE),
    Seal_Clupeid_Impact_Median=median(Seal_Clupeid_t / pmax(Clupeid_stock_t,.1),na.rm=TRUE),
    Corm_Clupeid_Impact_Median=median(Cormorant_Clupeid_t / pmax(Clupeid_stock_t,.1),na.rm=TRUE),
    Seal_Wrasse_Impact_Median=median(Seal_Wrasse_t / pmax(Wrasse_stock_t,.1),na.rm=TRUE),
    Corm_Wrasse_Impact_Median=median(Cormorant_Wrasse_t / pmax(Wrasse_stock_t,.1),na.rm=TRUE),
    Seal_CodJuv_C_Tonnes=median(Seal_Cod_Juveniles_t,na.rm=TRUE),
    Corm_CodJuv_C_Tonnes=median(Cormorant_Cod_Juveniles_t,na.rm=TRUE),
    Seal_CodAdult_C_Tonnes=median(Seal_Cod_Adults_t,na.rm=TRUE),
    Corm_CodAdult_C_Tonnes=median(Cormorant_Cod_Adults_t,na.rm=TRUE),
    Seal_Clupeid_C_Tonnes=median(Seal_Clupeid_t,na.rm=TRUE),
    Corm_Clupeid_C_Tonnes=median(Cormorant_Clupeid_t,na.rm=TRUE),
    Seal_Wrasse_C_Tonnes=median(Seal_Wrasse_t,na.rm=TRUE),
    Corm_Wrasse_C_Tonnes=median(Cormorant_Wrasse_t,na.rm=TRUE),
    Seal_Total_C_Tonnes=median(Seal_total_t,na.rm=TRUE),
    Corm_Total_C_Tonnes=median(Cormorant_total_t,na.rm=TRUE),
    .groups="drop")
  mc_results <- agg
  impact_plot_data <- mc_results %>%
    dplyr::select(Year, ends_with("_Impact_Median")) %>%
    pivot_longer(cols=-Year,names_to=c("Predator","Prey"),names_pattern="(.*)_(.*)_Impact_Median") %>%
    mutate(Predator=if_else(Predator=="Seal","Harbour Seals","Great Cormorants"),
           Prey=case_when(Prey=="CodJuv"~"Cod Juveniles",Prey=="CodAdult"~"Cod Adults",Prey=="Clupeid"~"Clupeidae",Prey=="Wrasse"~"Wrasses"),
           Prey=factor(Prey,levels=c("Cod Juveniles","Cod Adults","Clupeidae","Wrasses")))
  consumption_plot_data <- mc_results %>%
    dplyr::select(Year, matches("^(Seal|Corm)_(CodJuv|CodAdult|Clupeid|Wrasse)_C_Tonnes$")) %>%
    pivot_longer(cols=-Year,names_to=c("Predator","Prey"),names_pattern="(.*)_(.*)_C_Tonnes") %>%
    mutate(Predator=if_else(Predator=="Seal","Harbour Seals","Great Cormorants"),
           Prey=case_when(Prey=="CodJuv"~"Cod Juveniles",Prey=="CodAdult"~"Cod Adults",Prey=="Clupeid"~"Clupeidae",Prey=="Wrasse"~"Wrasses"),
           Prey=factor(Prey,levels=c("Cod Juveniles","Cod Adults","Clupeidae","Wrasses")))
p1 <- ggplot(impact_plot_data, aes(x = Year, y = value, color = Predator)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ Prey, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("Harbour Seals" = "#8c510a", "Great Cormorants" = "#01665e")) +
  scale_x_continuous(breaks = seq(2012, 2025, by = 3)) +
  theme_minimal() +
  labs(
    title = "A - Top-Down Impact Ratio",
    y = "Impact Ratio (Consumption / Biomass)",
    x = "Year"
  ) +
  theme(legend.position = "bottom", legend.title = element_blank(), strip.text = element_text(face = "bold"))

total_demand_data <- mc_results %>%
    dplyr::select(Year, Seal_Total_C_Tonnes, Corm_Total_C_Tonnes) %>%
    tidyr::pivot_longer(
      cols = c(Seal_Total_C_Tonnes, Corm_Total_C_Tonnes),
      names_to = "Predator",
      values_to = "Total_Demand"
    ) %>%
    dplyr::mutate(
      Predator = dplyr::if_else(Predator == "Seal_Total_C_Tonnes",
                                "Harbour Seals",
                                "Great Cormorants")
    )
p2 <- ggplot() +
    geom_area(data = consumption_plot_data, aes(x = Year, y = value, fill = Prey),
              alpha = 0.8, color = "white", linewidth = 0.2) +
    geom_line(data = total_demand_data, aes(x = Year, y = Total_Demand),
              color = "black", linetype = "dashed", linewidth = 1) +
    facet_wrap(~ Predator, scales = "free_y") +
    scale_fill_manual(values = c("Cod Juveniles" = "#74add1", "Cod Adults" = "#313695",
                                 "Clupeidae" = "darkorange", "Wrasses" = "seagreen")) +
    scale_x_continuous(breaks = seq(2012, 2025, by = 3)) +
    theme_minimal() +
    labs(
      title = "B - Total Prey Consumption ",
      y = "Consumption (Metric Tonnes)",
      x = "Year"
    ) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          strip.text = element_text(face = "bold"))

  all_prey_preds <- diet_predictions %>% filter(Model==model_name) %>%
    transmute(x=Year-2012,predicted=Predicted_RRA,Predator=if_else(Predator=="Seal","Harbour Seals","Great Cormorants"),Prey=Prey)

  all_prey_preds <- bind_rows(
    all_prey_preds %>% filter(Prey=="Cod",Predator=="Harbour Seals") %>% mutate(predicted=predicted*.65,Prey="Cod Juveniles"),
    all_prey_preds %>% filter(Prey=="Cod",Predator=="Harbour Seals") %>% mutate(predicted=predicted*.35,Prey="Cod Adults"),
    all_prey_preds %>% filter(Prey=="Cod",Predator=="Great Cormorants") %>% mutate(predicted=predicted*.95,Prey="Cod Juveniles"),
    all_prey_preds %>% filter(Prey=="Cod",Predator=="Great Cormorants") %>% mutate(predicted=predicted*.05,Prey="Cod Adults"),
    all_prey_preds %>% filter(Prey!="Cod") %>% mutate(Prey=recode(Prey,Clupeid="Clupeidae",Wrasse="Wrasses"))
  ) %>% mutate(Prey=factor(Prey,levels=c("Cod Juveniles","Cod Adults","Clupeidae","Wrasses")),
              conf.low=NA_real_,conf.high=NA_real_)

  prey_plot_data <- stock_components %>%
    filter(Year>=2012,Year<=2025) %>%
    mutate(Stock=recode(Stock,Cod_Juveniles="Cod Juveniles",Cod_Adults="Cod Adults",Wrasse="Wrasses",NSAS_Herring="NSAS Herring",WBSS_Herring="WBSS Herring"),
           Skagerrak_Biomass_Median=Combined_biomass_t*Regional_ratio,
           Skagerrak_Biomass_LCI=pmax(Skagerrak_Biomass_Median-1.96*Combined_biomass_sd_t*Regional_ratio,0),
           Skagerrak_Biomass_UCI=Skagerrak_Biomass_Median+1.96*Combined_biomass_sd_t*Regional_ratio)
  predator_plot_data <- predator_demographics %>%
    mutate(Seal_LCI=Seal_N_mean-1.96*Seal_N_sd,Seal_UCI=Seal_N_mean+1.96*Seal_N_sd,
           Corm_LCI=Corm_N_mean-1.96*Corm_N_sd,Corm_UCI=Corm_N_mean+1.96*Corm_N_sd)

diet_trend_plot <- ggplot(all_prey_preds, aes(x = x + 2012, y = predicted, color = Predator, fill = Predator)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, color = NA) +
  facet_wrap(~ Prey, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("Harbour Seals" = "#8c510a", "Great Cormorants" = "#01665e")) +
  scale_fill_manual(values = c("Harbour Seals" = "#8c510a", "Great Cormorants" = "#01665e")) +
  scale_x_continuous(breaks = seq(2012, 2025, by = 3)) +
  theme_minimal() +
  labs(
    title = "Season-Adjusted Dietary RRA Trends Over Time",
    x = "Year",
    y = "Predicted Relative Read Abundance (RRA)"
  ) +
  theme(legend.position = "bottom", legend.title = element_blank(), strip.text = element_text(face = "bold"))

prey_dynamics_plot <- ggplot(prey_plot_data, aes(x = Year, y = Skagerrak_Biomass_Median, color = Stock, fill = Stock)) +
  geom_line(linewidth = 1.1) +
  geom_ribbon(aes(ymin = Skagerrak_Biomass_LCI, ymax = Skagerrak_Biomass_UCI), alpha = 0.2, color = NA) +
  facet_wrap(~ Stock, scales = "free_y", ncol = 3) +
  scale_x_continuous(breaks = seq(2012, 2025, by = 3)) +
  scale_color_viridis_d(option = "turbo") +
  scale_fill_viridis_d(option = "turbo") +
  theme_minimal() +
  labs(
    title = "Estimated Prey Population Dynamics in the Skagerrak",
    x = "Year",
    y = "Standing Stock Biomass (Metric Tonnes)"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines")
  )

seal_pop_plot <- ggplot(predator_plot_data, aes(x = Year, y = Seal_N_mean)) +
  geom_line(color = "#8c510a", linewidth = 1.2) +
  geom_ribbon(aes(ymin = Seal_LCI, ymax = Seal_UCI), fill = "#8c510a", alpha = 0.2) +
  scale_x_continuous(breaks = seq(2012, 2025, by = 3)) +
  theme_minimal() +
  labs(
    title = "Harbour Seal Population Dynamics",
    x = "Year",
    y = "Number of Individuals"
  )

cormorant_pop_plot <- ggplot(predator_plot_data, aes(x = Year, y = Corm_N_mean)) +
  geom_line(color = "#01665e", linewidth = 1.2) +
  geom_ribbon(aes(ymin = Corm_LCI, ymax = Corm_UCI), fill = "#01665e", alpha = 0.2) +
  scale_x_continuous(breaks = seq(2012, 2025, by = 3)) +
  theme_minimal() +
  labs(
    title = "Great Cormorant Population Dynamics",
    x = "Year",
    y = "Number of Individuals"
  )

  # Fisheries extraction comparison
  get_annual_catch_safe <- function(filename, ratio) {
    f <- file.path(final_dir, "Tables", "Monte_Carlo", filename)
    if (!file.exists(f)) return(tibble(Year = 2012:2025, Catch_Val = NA_real_))
    z <- readr::read_csv(f, show_col_types = FALSE)
    y <- intersect(c("Year", "year"), names(z))
    c <- intersect(c("catches", "Catches", "landings", "Landings"), names(z))
    if (!length(y) || !length(c)) return(tibble(Year = 2012:2025, Catch_Val = NA_real_))
    z %>% transmute(Year = as.integer(.data[[y[1]]]), Catch_Val = as.numeric(.data[[c[1]]]) * ratio) %>%
      filter(Year %in% 2012:2025) %>% distinct(Year, .keep_all = TRUE)
  }

  annual_comm_cod <- get_annual_catch_safe("ICES_raw_Cod_Adults.csv", .06) %>% rename(Comm_Cod = Catch_Val)
  annual_sprat <- get_annual_catch_safe("ICES_raw_Sprat.csv", .15) %>% rename(Sprat = Catch_Val)
  annual_nsas <- get_annual_catch_safe("ICES_raw_NSAS_Herring.csv", .10) %>% rename(NSAS = Catch_Val)
  annual_wbss <- get_annual_catch_safe("ICES_raw_WBSS_Herring.csv", .40) %>% rename(WBSS = Catch_Val)

  annual_human <- tibble(Year = 2012:2025) %>%
    left_join(annual_comm_cod, by = "Year") %>%
    left_join(annual_sprat, by = "Year") %>%
    left_join(annual_nsas, by = "Year") %>%
    left_join(annual_wbss, by = "Year") %>%
    mutate(
      Comm_Pelagic = Sprat + NSAS + WBSS,
      # FIX: Drop recreational cod line entirely after 2022 (turns into NA so line stops)
      Rec_Cod = if_else(Year > 2022, NA_real_, Comm_Cod * 3)
    )

  annual_predator <- mc_results %>% transmute(
    Year,
    Corm_Total = Corm_Total_C_Tonnes,
    Seal_Total = Seal_Total_C_Tonnes,
    Corm_Cod = Corm_CodJuv_C_Tonnes + Corm_CodAdult_C_Tonnes,
    Seal_Cod = Seal_CodJuv_C_Tonnes + Seal_CodAdult_C_Tonnes
  )

  plot_data_time <- annual_human %>%
    left_join(annual_predator, by = "Year") %>%
    pivot_longer(cols = -Year, names_to = "Variable", values_to = "Biomass_tonnes") %>%
    filter(Variable %in% c("Corm_Total", "Corm_Cod", "Seal_Total", "Seal_Cod", "Comm_Pelagic", "Comm_Cod", "Rec_Cod")) %>%
    mutate(
      Source = case_when(
        Variable %in% c("Corm_Total", "Corm_Cod") ~ "Great Cormorant",
        Variable %in% c("Seal_Total", "Seal_Cod") ~ "Harbour Seal",
        Variable %in% c("Comm_Pelagic", "Comm_Cod") ~ "Commercial Fleet",
        Variable == "Rec_Cod" ~ "Recreational Fishery"
      ),
      Category = case_when(
        Variable %in% c("Corm_Total", "Seal_Total", "Comm_Pelagic") ~ "Total Community Extraction",
        TRUE ~ "Targeted Atlantic Cod Extraction"
      ),
      Type = case_when(
        Variable %in% c("Corm_Total", "Corm_Cod") ~ "Avian Predation",
        Variable %in% c("Seal_Total", "Seal_Cod") ~ "Mammalian Predation",
        TRUE ~ "Human Extraction"
      ),
      Source = factor(Source, levels = c("Great Cormorant", "Harbour Seal", "Commercial Fleet", "Recreational Fishery"))
    )

  extraction_time_plot <- ggplot(plot_data_time, aes(x = Year, y = Biomass_tonnes, color = Source, linetype = Type)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    facet_wrap(~ Category, scales = "free_y") +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_color_manual(values = c("Great Cormorant" = "#01665e",
                                  "Harbour Seal" = "#8c510a",
                                  "Commercial Fleet" = "grey20",
                                  "Recreational Fishery" = "grey50")) +
    scale_linetype_manual(values = c("Avian Predation" = "solid",
                                     "Mammalian Predation" = "solid",
                                     "Human Extraction" = "dashed")) +
    scale_x_continuous(breaks = seq(2012, 2025, by = 2)) +
    labs(
      title = "Temporal Dynamics of Biomass Extraction",
      y = "Extracted Biomass (Tonnes)",
      x = "Year",
      color = "Extraction Agent",
      linetype = "Extraction Type"
    ) +
    theme_bw(base_size = 14) +
    theme(
      strip.background = element_rect(fill = "grey90", color = "black"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical"
    )

  ggplot2::ggsave(file.path(original_dir, "01_Top_down_impact.png"), p1, width = 15, height = 5, dpi = 300)
  ggplot2::ggsave(file.path(original_dir, "02_Total_prey_consumption.png"), p2, width = 12, height = 6, dpi = 300)
  ggplot2::ggsave(file.path(original_dir, "03_Original_A_B_combined.png"), p1 / p2, width = 16, height = 12, dpi = 300)
  ggplot2::ggsave(file.path(original_dir, "04_Dietary_RRA_trends.png"), diet_trend_plot, width = 15, height = 5, dpi = 300)
  ggplot2::ggsave(file.path(original_dir, "05_Prey_population_dynamics.png"), prey_dynamics_plot, width = 13, height = 10, dpi = 300)
  ggplot2::ggsave(file.path(original_dir, "06_Predator_population_dynamics.png"), seal_pop_plot + cormorant_pop_plot, width = 12, height = 5, dpi = 300)
  ggplot2::ggsave(file.path(original_dir, "07_Temporal_biomass_extraction.png"), extraction_time_plot, width = 12, height = 7, dpi = 300)
}
message("Original plotting layouts saved in: ", file.path(final_dir, "Figures"))
