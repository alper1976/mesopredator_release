# PORTABLE STANDALONE RECONSTRUCTION. Original calculations retained.
# ============================================================
# 07_Trophic_Network.R
# Mesopredator_check_2 / standalone manuscript reconstruction
# Reconstructs original Tarjei three-level predator-fish-invertebrate network.
# IMPORTANT: Co-occurrence is not proof of direct predation.
# Inputs: validated script 02 phyloseq RDS files only.
# Manuscript output: one PNG in outputs/07_Trophic_Network/figures.
# Method/source audits: ../../docs/validation/07_Trophic_Network.
# Not run in this environment; validate against original Figure 4 locally.
# ============================================================

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
INPUT <- file.path(BASE, "processed_data", "03_Phyloseq_Preprocessing")
MIFISH <- file.path(INPUT, "ps_mifish_global_0_1pct_identified.rds")
COI <- file.path(INPUT, "ps_coi_global_0_1pct_identified.rds")
OUT <- file.path(BASE, "processed_data", "07_Trophic_Network")
FIG <- file.path(OUT, "figures")
TAB <- file.path(OUT, "tables")
VAL <- file.path(BASE, "validation", "07_Trophic_Network")
for (dir in c(FIG, TAB, VAL)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
if (!all(file.exists(c(MIFISH, COI)))) stop("Missing script 02 phyloseq inputs")

# ------------------------------------------------------------
# 2. PACKAGES
# ------------------------------------------------------------

for (pkg in c("phyloseq", "dplyr", "ggplot2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing package: ", pkg)
  }
}

library(phyloseq)
library(dplyr)
library(ggplot2)


# ------------------------------------------------------------
# 3. LOAD ORIGINAL RECONSTRUCTED DATA
# ------------------------------------------------------------

ps_mf <- readRDS(MIFISH)
ps_coi <- readRDS(COI)

stopifnot(
  inherits(ps_mf, "phyloseq"),
  inherits(ps_coi, "phyloseq")
)

# Exclude unidentified MiFish ASVs as in original ps_01.
# Detection-based calculations do not require re-normalization.

if (!"LCA_taxon" %in%
    colnames(as(tax_table(ps_mf), "matrix"))) {
  stop("MiFish LCA_taxon not available")
}

ps_mf <- subset_taxa(
  ps_mf,
  LCA_taxon != "no identification"
)


# ------------------------------------------------------------
# 4. EXTRACT AND AGGREGATE LCA DATA
# ------------------------------------------------------------

extract_data <- function(ps) {
  
  m <- as(otu_table(ps), "matrix")
  
  if (taxa_are_rows(ps)) {
    m <- t(m)
  }
  
  tx <- as.data.frame(
    tax_table(ps),
    check.names = FALSE
  )
  
  md <- as.data.frame(
    sample_data(ps),
    check.names = FALSE
  )
  
  stopifnot(
    all(colnames(m) %in% rownames(tx)),
    all(rownames(m) %in% rownames(md))
  )
  
  tx <- tx[colnames(m), , drop = FALSE]
  md <- md[rownames(m), , drop = FALSE]
  
  stopifnot(
    all(c("LCA_taxon", "Class") %in% names(tx)),
    "species" %in% names(md)
  )
  
  labels <- gsub(
    "_",
    " ",
    as.character(tx$LCA_taxon),
    fixed = TRUE
  )
  
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("Missing LCA labels")
  }
  
  aggregated <- t(
    rowsum(
      t(m),
      group = labels,
      reorder = FALSE
    )
  )
  
  list(
    otu = aggregated,
    tx = tx,
    md = md,
    labels = labels
  )
}

mi <- extract_data(ps_mf)
co <- extract_data(ps_coi)


# ------------------------------------------------------------
# 5. FIND SHARED SAMPLES
# ------------------------------------------------------------

shared <- intersect(
  rownames(mi$otu),
  rownames(co$otu)
)

if (!length(shared)) {
  stop("No shared sample IDs")
}

cat("\nShared samples:", length(shared), "\n")


# ------------------------------------------------------------
# 6. SELECT COI INVERTEBRATES
# ------------------------------------------------------------

class_by_label <- split(
  as.character(co$tx$Class),
  co$labels
)

class_collision <- vapply(
  class_by_label,
  function(x) length(unique(x)) > 1L,
  logical(1)
)

if (any(class_collision)) {
  stop("Conflicting classes for aggregated COI LCA")
}

invertebrates <- names(class_by_label)[
  vapply(
    class_by_label,
    function(x) unique(x) != "Actinopteri",
    logical(1)
  )
]

invertebrates <- intersect(
  colnames(co$otu),
  invertebrates
)

m <- mi$otu[shared, , drop = FALSE]

v <- co$otu[
  shared,
  invertebrates,
  drop = FALSE
]


# ------------------------------------------------------------
# 7. VALIDATE PREDATOR METADATA - FIXED
# ------------------------------------------------------------

# Extract a single predator label per sample.
# The previous error was caused by subsetting the full
# metadata data frame and converting it to character.

md <- as.data.frame(
  mi$md[shared, , drop = FALSE]
)

co_md <- as.data.frame(
  co$md[shared, , drop = FALSE]
)

mifish_species <- trimws(
  as.character(md[["species"]])
)

coi_species <- trimws(
  as.character(co_md[["species"]])
)

# Normalize only the documented cormorant naming variant.

mifish_species[
  mifish_species == "top_cormorant"
] <- "cormorant"

coi_species[
  coi_species == "top_cormorant"
] <- "cormorant"

stopifnot(
  length(mifish_species) == length(shared),
  length(coi_species) == length(shared),
  !anyNA(mifish_species),
  !anyNA(coi_species),
  identical(mifish_species, coi_species),
  all(mifish_species %in% c("cormorant", "harbor_seal"))
)

md$species <- mifish_species

cat("\n=== VALIDATED PREDATOR METADATA ===\n")
print(table(md$species))

metadata_audit <- data.frame(
  sample_id = shared,
  MiFish_original = as.character(
    mi$md[shared, , drop = FALSE][["species"]]
  ),
  COI_original = as.character(
    co$md[shared, , drop = FALSE][["species"]]
  ),
  validated_predator = mifish_species,
  stringsAsFactors = FALSE
)

write.csv(
  metadata_audit,
  file.path(VAL, "07_PREDATOR_METADATA_VALIDATION.csv"),
  row.names = FALSE
)

cat("\n=== ABUNDANCE VALIDATION ===\n")

cat("Shared samples:", nrow(m), "\n")
cat("Samples with zero MiFish abundance:",
    sum(rowSums(m) == 0), "\n")
cat("Samples with zero COI invertebrate abundance:",
    sum(rowSums(v) == 0), "\n")

if (
  anyNA(m) ||
  anyNA(v) ||
  any(!is.finite(m)) ||
  any(!is.finite(v)) ||
  any(m < 0) ||
  any(v < 0)
) {
  stop("Invalid abundance values: NA, Inf or negative counts")
}

cat("\nAbundance validation passed.\n")


# ------------------------------------------------------------
# 8. ORIGINAL TAXON NAME SHORTENING
# ------------------------------------------------------------

short_fish <- function(x) {
  
  vapply(x, function(z) {
    
    if (grepl("Myoxocephalus scorpius Myoxocephalus", z))
      return("Myoxocephalus scorpius")
    
    if (grepl("Pleuronectes platessa Platichthys", z))
      return("Pleuronectes platessa")
    
    if (grepl("Platichthys flesus Pleuronectes platessa Hippoglossoides", z))
      return("Platichthys flesus")
    
    if (grepl("Pollachius virens Pollachius|Pollachius pollachius", z))
      return("Pollachius spp.")
    
    if (grepl("Limanda limanda Hippoglossus", z))
      return("Limanda limanda")
    
    if (grepl("Clupea harengus.*Clupea pallasii|Clupea pallasii.*Clupea harengus", z))
      return("Clupea harengus")
    
    if (grepl("Hyperoplus.*Ammodytes|Ammodytes.*Hyperoplus", z))
      return("Ammodytes spp.")
    
    if (grepl("Sprattus sprattus Clupea|Clupea pallasii Sprattus", z))
      return("Sprattus sprattus")
    
    paste(
      head(strsplit(z, " ", fixed = TRUE)[[1]], 2),
      collapse = " "
    )
    
  }, character(1), USE.NAMES = FALSE)
}


short_inv <- function(x) {
  
  vapply(x, function(z) {
    
    if (grepl("Evadne", z))
      return("Evadne sp.")
    
    if (grepl("Alitta|Hediste", z))
      return("Alitta / Hediste")
    
    if (grepl("Echinogammarus.*Gammarus", z))
      return("Echinogammarus sp.")
    
    if (grepl("Heteromastus|Capitellidae", z))
      return("Heteromastus sp.")
    
    if (grepl("Palaemon", z))
      return("Palaemon adspersus")
    
    if (grepl("Harmothoe", z))
      return("Harmothoe imbricata")
    
    if (grepl("Liocarcinus", z))
      return("Liocarcinus holsatus")
    
    if (grepl("Penilia", z))
      return("Penilia avirostris")
    
    if (grepl("Acartia", z))
      return("Acartia spp.")
    
    if (grepl("Cyclopoida", z))
      return("Cyclopoida sp.")
    
    paste(
      head(strsplit(z, " ", fixed = TRUE)[[1]], 2),
      collapse = " "
    )
    
  }, character(1), USE.NAMES = FALSE)
}


# ------------------------------------------------------------
# 9. PREDATOR-SPECIFIC MiFish FOO
# ------------------------------------------------------------

corm <- md$species == "cormorant"
seal <- md$species == "harbor_seal"

stopifnot(any(corm), any(seal))

foo_corm <- colSums(
  m[corm, , drop = FALSE] > 0
) / sum(corm) * 100

foo_seal <- colSums(
  m[seal, , drop = FALSE] > 0
) / sum(seal) * 100

fish_data <- data.frame(
  taxon = colnames(m),
  taxon_short = short_fish(colnames(m)),
  foo_corm = as.numeric(foo_corm),
  foo_seal = as.numeric(foo_seal),
  mean_foo = (
    as.numeric(foo_corm) +
      as.numeric(foo_seal)
  ) / 2
)


# ------------------------------------------------------------
# 10. TOP EIGHT FISH
# ------------------------------------------------------------

# Original operation: aggregate shortened names with max().

fish_foo <- fish_data %>%
  group_by(taxon_short) %>%
  summarise(
    foo_corm = max(foo_corm),
    foo_seal = max(foo_seal),
    mean_foo = max(mean_foo),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_foo)) %>%
  slice_head(n = 8)

if (nrow(fish_foo) != 8L) {
  stop("Fewer than eight fish taxa")
}

cat("\n=== TOP EIGHT MiFish FISH ===\n")
print(fish_foo)


# ------------------------------------------------------------
# 11. FISH-INVERTEBRATE CO-OCCURRENCE
# ------------------------------------------------------------

edges <- list()

for (fish in fish_foo$taxon_short) {
  
  cols <- colnames(m)[
    short_fish(colnames(m)) == fish
  ]
  
  fish_present <- rowSums(
    m[, cols, drop = FALSE]
  ) > 0
  
  denominator <- sum(fish_present)
  
  if (denominator == 0L) next
  
  for (inv_name in colnames(v)) {
    
    co_count <- sum(
      fish_present &
        v[, inv_name] > 0
    )
    
    if (co_count > 0L) {
      
      edges[[length(edges) + 1L]] <- data.frame(
        fish = fish,
        inv = short_inv(inv_name),
        source_invertebrate = inv_name,
        cooccur = co_count,
        pct = 100 * co_count / denominator,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (!length(edges)) {
  stop("No fish-invertebrate co-occurrences")
}

raw_edges <- dplyr::bind_rows(edges)

# Preserve original grouping and max() logic.

fish_inv_edges <- raw_edges %>%
  group_by(fish, inv) %>%
  summarise(
    cooccur = max(cooccur),
    pct = max(pct),
    .groups = "drop"
  ) %>%
  filter(pct >= 10)

cat("\nFiltered fish-invertebrate edges:",
    nrow(fish_inv_edges), "\n")


# ------------------------------------------------------------
# 12. PREDATOR-FISH EDGES
# ------------------------------------------------------------

# Original strict threshold: FOO > 5%.

pred_fish <- dplyr::bind_rows(
  
  data.frame(
    from = "Great cormorant",
    to = fish_foo$taxon_short,
    weight = fish_foo$foo_corm,
    predator = "Great cormorant"
  ),
  
  data.frame(
    from = "Harbour seal",
    to = fish_foo$taxon_short,
    weight = fish_foo$foo_seal,
    predator = "Harbour seal"
  )
  
) %>%
  filter(weight > 5)


# ------------------------------------------------------------
# 13. EXPORT NETWORK DATA
# ------------------------------------------------------------

write.csv(
  fish_foo,
  file.path(TAB, "07_TOP8_FISH_FOO.csv"),
  row.names = FALSE
)

write.csv(
  fish_data,
  file.path(VAL, "07_ALL_FISH_FOO_SOURCE.csv"),
  row.names = FALSE
)

write.csv(
  raw_edges,
  file.path(VAL, "07_RAW_FISH_INVERTEBRATE_COOCCURRENCE.csv"),
  row.names = FALSE
)

write.csv(
  fish_inv_edges,
  file.path(TAB, "07_FILTERED_FISH_INVERTEBRATE_EDGES.csv"),
  row.names = FALSE
)

write.csv(
  pred_fish,
  file.path(TAB, "07_PREDATOR_FISH_EDGES.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    sample_id = shared,
    predator = md$species
  ),
  file.path(VAL, "07_SHARED_SAMPLE_AUDIT.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------
# 14. NETWORK NODES
# ------------------------------------------------------------

pred_nodes <- data.frame(
  name = c("Great cormorant", "Harbour seal"),
  y = 3,
  x = c(2.5, 6.5),
  type = "predator"
)

fish_nodes <- data.frame(
  name = fish_foo$taxon_short,
  y = 2,
  x = seq(1, 9, length.out = nrow(fish_foo)),
  type = "fish"
)

inv_names <- unique(fish_inv_edges$inv)

inv_nodes <- data.frame(
  name = inv_names,
  y = 1,
  x = if (length(inv_names)) {
    seq(1, 9, length.out = length(inv_names))
  } else {
    numeric()
  },
  type = "invertebrate"
)

nodes <- dplyr::bind_rows(
  pred_nodes,
  fish_nodes,
  inv_nodes
)


# ------------------------------------------------------------
# 15. LINK EDGES TO NODE COORDINATES
# ------------------------------------------------------------

coords <- function(df) {
  
  df %>%
    left_join(
      nodes %>% select(name, x, y),
      by = c("from" = "name")
    ) %>%
    rename(
      x_from = x,
      y_from = y
    ) %>%
    left_join(
      nodes %>% select(name, x, y),
      by = c("to" = "name")
    ) %>%
    rename(
      x_to = x,
      y_to = y
    )
}

pred_plot <- coords(pred_fish)

inv_plot <- coords(
  fish_inv_edges %>%
    transmute(
      from = fish,
      to = inv,
      weight = pct
    )
)

if (
  anyNA(pred_plot$x_from) ||
  anyNA(pred_plot$x_to) ||
  anyNA(inv_plot$x_from) ||
  anyNA(inv_plot$x_to)
) {
  stop("Node/edge name mismatch")
}


# ------------------------------------------------------------
# 16. RECONSTRUCT NETWORK FIGURE
# ------------------------------------------------------------

p <- ggplot() +
  
  geom_segment(
    data = pred_plot,
    aes(
      x = x_from,
      y = y_from,
      xend = x_to,
      yend = y_to,
      colour = predator,
      linewidth = weight
    ),
    alpha = 0.4
  ) +
  
  geom_segment(
    data = inv_plot,
    aes(
      x = x_from,
      y = y_from,
      xend = x_to,
      yend = y_to,
      linewidth = weight
    ),
    colour = "grey40",
    alpha = 0.6
  ) +
  
  geom_point(
    data = nodes,
    aes(x = x, y = y, fill = type),
    shape = 21,
    size = 10,
    colour = "grey30",
    stroke = 0.8
  ) +
  
  geom_text(
    data = pred_nodes,
    aes(x = x, y = y + 0.18, label = name),
    size = 3.2,
    fontface = "bold"
  ) +
  
  geom_text(
    data = fish_nodes,
    aes(x = x, y = y - 0.15, label = name),
    size = 2.5,
    fontface = "italic",
    angle = 35,
    hjust = 1,
    vjust = 1
  ) +
  
  geom_text(
    data = inv_nodes,
    aes(x = x, y = y - 0.15, label = name),
    size = 2.5,
    fontface = "italic",
    angle = 35,
    hjust = 1,
    vjust = 1
  ) +
  
  annotate(
    "text",
    x = -0.2, y = 3,
    label = "Predator",
    hjust = 1,
    size = 3.5,
    fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = -0.2, y = 2,
    label = "Primary\nprey",
    hjust = 1,
    size = 3.5,
    fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = -0.2, y = 1,
    label = "Secondary\nprey",
    hjust = 1,
    size = 3.5,
    fontface = "bold"
  ) +
  
  scale_colour_manual(
    values = c(
      "Great cormorant" = "#d6604d",
      "Harbour seal" = "#2166ac"
    ),
    name = "Predator"
  ) +
  
  scale_fill_manual(
    values = c(
      predator = "white",
      fish = "#d9f0d3",
      invertebrate = "#e8d5f5"
    ),
    guide = "none"
  ) +
  
  scale_linewidth_continuous(
    range = c(0.3, 3),
    name = "Frequency (%)",
    breaks = c(10, 20, 30, 40)
  ) +
  
  coord_cartesian(
    xlim = c(-0.5, 10.5),
    ylim = c(0.5, 3.5)
  ) +
  
  theme_void(base_size = 11) +
  
  theme(
    legend.position = "right",
    plot.margin = margin(20, 20, 80, 60),
    plot.title = element_text(
      face = "bold",
      size = 12
    )
  ) +
  
  labs(
    title = "Trophic network - primary and secondary prey",
    subtitle = paste(
      "Edge width: predator-fish = FOO (%);",
      "fish-invertebrate = co-occurrence (%)"
    )
  )


# ------------------------------------------------------------
# 17. SAVE FIGURES
# ------------------------------------------------------------

ggsave(
  file.path(
    FIG,
    "trophic_network.png"
  ),
  p,
  width = 14,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(
    VAL,
    "07_TROPHIC_NETWORK_VALIDATION.pdf"
  ),
  p,
  width = 14,
  height = 8,
  bg = "white"
)


# ------------------------------------------------------------
# 18. FINAL RECONSTRUCTION SUMMARY
# ------------------------------------------------------------

summary <- data.frame(
  item = c(
    "MiFish input ASVs",
    "COI input ASVs",
    "Shared samples",
    "Shared cormorant",
    "Shared harbor seal",
    "Top fish nodes",
    "Secondary invertebrate nodes",
    "Predator-fish edges",
    "Fish-invertebrate edges"
  ),
  value = c(
    ntaxa(ps_mf),
    ntaxa(ps_coi),
    length(shared),
    sum(corm),
    sum(seal),
    nrow(fish_foo),
    length(inv_names),
    nrow(pred_fish),
    nrow(fish_inv_edges)
  )
)

write.csv(
  summary,
  file.path(VAL, "07_RECONSTRUCTION_SUMMARY.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "Source: Tarjei Network.analysis.Rmd and independent script 28 reconstruction.",
    "Original three-level trophic network reconstruction.",
    "Matched samples across MiFish and COI by exact sample ID.",
    "Predator metadata validated for all shared samples.",
    "Top eight MiFish fish taxa selected by mean FOO.",
    "Shortened taxon names aggregated using max().",
    "Predator-fish threshold: FOO > 5%.",
    "Fish-invertebrate threshold: co-occurrence >= 10%.",
    "Fish-invertebrate associations represent co-detection.",
    "Co-detection does not establish direct consumption.",
    "MiFish and COI inputs are identified global 0.1 percent phyloseq objects from script 02.",
    "Final visual fidelity requires comparison with Tarjei original trophic network Figure 4."
  ),
  file.path(VAL, "07_METHOD_AND_LIMITATIONS.txt")
)



# ------------------------------------------------------------
# 18B. PROVENANCE, FINAL OUTPUT CHECK, AND QA
# ------------------------------------------------------------
stopifnot(file.exists(file.path(FIG, "trophic_network.png")))
stopifnot(all(file.exists(file.path(TAB, c(
  "07_TOP8_FISH_FOO.csv",
  "07_FILTERED_FISH_INVERTEBRATE_EDGES.csv",
  "07_PREDATOR_FISH_EDGES.csv"
)))))
write.csv(data.frame(
  input = c("MiFish", "COI"),
  path = c(MIFISH, COI),
  md5 = unname(tools::md5sum(c(MIFISH, COI)))
), file.path(VAL, "07_INPUT_PROVENANCE.csv"), row.names = FALSE)
write.csv(data.frame(
  file = c("trophic_network.png", "07_TOP8_FISH_FOO.csv",
           "07_FILTERED_FISH_INVERTEBRATE_EDGES.csv",
           "07_PREDATOR_FISH_EDGES.csv"),
  kind = c("figure", "table", "table", "table")
), file.path(VAL, "07_OUTPUT_MANIFEST.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(VAL, "07_SESSION_INFO.txt"))

# ------------------------------------------------------------
# 19. CONSOLE REPORT
# ------------------------------------------------------------

cat("\n======================================\n")
cat("SCRIPT 06 COMPLETE\n")
cat("======================================\n")

cat("\nReconstruction summary:\n")
print(summary, row.names = FALSE)

cat("\nTop eight fish taxa:\n")
print(fish_foo)

cat("\nOutput directory:\n")
cat(FIG, "\n")

cat("\nOriginal phyloseq objects were not modified.\n")
