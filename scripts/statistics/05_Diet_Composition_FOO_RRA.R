# PORTABLE STANDALONE RECONSTRUCTION. Original calculations retained.
# 05_Diet_Composition_FOO_RRA.R
# Mesopredator_paper / Mesopredator_check_2
# Faithful reproduction candidate from Tarjei Rostum:
# prevalence.distribution.coi.mifish.Rmd, original R chunks 4,6-12.
# Raw original Rmd code is preserved below in sequence where possible.
# Final figure: diet_FOO_four_panels.png (A-D combined).
# Individual panel previews are validation only; three original taxa_list_*.csv.
# Important: original source SAVES diet_composition_mifish.png twice;
# its LAST plot supersedes the first. We export only the final version.
# No previously reconstructed data are inputs: only script 02 outputs.
# Not executed here: validate numerical/visual fidelity in local RStudio.

required <- c('phyloseq','dplyr','tidyr','stringr','ggplot2','patchwork')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Install required R packages: ',paste(missing,collapse=', '))
suppressPackageStartupMessages({
 library(phyloseq);library(dplyr);library(tidyr);library(stringr)
 library(ggplot2);library(patchwork)
})
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
OUT <- file.path(BASE,'processed_data','05_Diet_Composition_FOO_RRA')
figs_dir <- file.path(OUT,'figures')
home_dir <- file.path(OUT,'tables')
VAL <- file.path(BASE,'validation','05_Diet_Composition_FOO_RRA')
for(p in c(figs_dir,home_dir,VAL)) dir.create(p,recursive=TRUE,showWarnings=FALSE)
mi_path <- file.path(INPUT,'ps_mifish_global_0_1pct_identified.rds')
coi_path <- file.path(INPUT,'ps_coi_global_0_1pct_identified.rds')
if(!all(file.exists(c(mi_path,coi_path)))) stop('Missing script 02 identified phyloseq outputs')
ps_01 <- readRDS(mi_path)
ps_coi_01 <- readRDS(coi_path)
stopifnot(inherits(ps_01,'phyloseq'),inherits(ps_coi_01,'phyloseq'))
stopifnot(nsamples(ps_01)==161L,ntaxa(ps_01)==69L,
          nsamples(ps_coi_01)==145L,ntaxa(ps_coi_01)==172L)
# Ensure ASV order is identical in otu_table and tax_table because original
# code renames columns positionally and assumes order correspondence.
check_alignment <- function(ps,tag) {
 mat <- as(otu_table(ps),'matrix')
 if(taxa_are_rows(ps)) mat <- t(mat)
 tx <- as.data.frame(tax_table(ps))
 if(!identical(colnames(mat),rownames(tx))) {
   stop(tag,': OTU/taxonomy column order differs; original positional renaming unsafe')
 }
 md <- as.data.frame(sample_data(ps))
 if(!identical(rownames(mat),rownames(md))) {
   stop(tag,': OTU/metadata row order differs')
 }
 invisible(TRUE)
}
check_alignment(ps_01,'MiFish')
check_alignment(ps_coi_01,'COI')
mdmi <- as.data.frame(sample_data(ps_01))
mdcoi <- as.data.frame(sample_data(ps_coi_01))
stopifnot(sum(as.character(mdmi$species)=='cormorant')==86L,
 sum(as.character(mdmi$species)=='harbor_seal')==75L,
 sum(as.character(mdcoi$species)=='cormorant')==80L,
 sum(as.character(mdcoi$species)=='harbor_seal')==65L)
write.csv(data.frame(source=c('MiFish','COI'),file=c(mi_path,coi_path),
  md5=unname(tools::md5sum(c(mi_path,coi_path)))),
  file.path(VAL,'05_INPUT_PROVENANCE.csv'),row.names=FALSE)
cat('04 input validated: MiFish 161x69; COI 145x172\\n')


# ---- Original two-panel prey detection, chunk 4 ----
# ============================================================
# PREY SPECIES DETECTED — two-panel FOO lollipop plot
# Panel A: MiFish top 20 fish taxa
# Panel B: COI top 20 taxa (fish + invertebrates mixed)
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ----------------------------------------------------------
# SHARED SETTINGS
# ----------------------------------------------------------

pred_colours <- c("Great cormorant" = "#d6604d",
                  "Harbour seal"    = "#2166ac")

shorten_taxon <- function(x) {
  sapply(strsplit(x, " "), function(w) {
    if (length(w) >= 2) paste(w[1:2], collapse = " ") else w[1]
  })
}

# ----------------------------------------------------------
# PANEL A: MiFish — top 20 fish taxa
# ----------------------------------------------------------

tt_mf <- as.data.frame(tax_table(ps_01))
tt_mf$LCA_taxon_clean <- gsub("_", " ", tt_mf$LCA_taxon)
sd_mf  <- as.data.frame(sample_data(ps_01))
otu_mf <- as.data.frame(otu_table(ps_01))
if (taxa_are_rows(ps_01)) otu_mf <- as.data.frame(t(otu_mf))
colnames(otu_mf) <- tt_mf$LCA_taxon_clean

corm_rows_mf <- rownames(sd_mf)[sd_mf$species == "cormorant"]
seal_rows_mf <- rownames(sd_mf)[sd_mf$species == "harbor_seal"]
n_corm_mf    <- length(corm_rows_mf)
n_seal_mf    <- length(seal_rows_mf)

foo_corm_mf <- colSums(otu_mf[rownames(otu_mf) %in% corm_rows_mf, ] > 0) / n_corm_mf * 100
foo_seal_mf <- colSums(otu_mf[rownames(otu_mf) %in% seal_rows_mf, ] > 0) / n_seal_mf * 100

foo_mf_df <- data.frame(
  taxon             = tt_mf$LCA_taxon_clean,
  taxon_short       = shorten_taxon(tt_mf$LCA_taxon_clean),
  `Great cormorant` = as.numeric(foo_corm_mf),
  `Harbour seal`    = as.numeric(foo_seal_mf),
  check.names       = FALSE
) %>%
  mutate(
    taxon_short = ifelse(
      grepl("Clupea harengus.*Clupea pallasii|Clupea pallasii.*Clupea harengus", taxon),
      "Clupea harengus", taxon_short),
    mean_foo = (`Great cormorant` + `Harbour seal`) / 2
  ) %>%
  group_by(taxon_short) %>%
  summarise(
    `Great cormorant` = max(`Great cormorant`),
    `Harbour seal`    = max(`Harbour seal`),
    mean_foo          = max(mean_foo),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_foo)) %>%
  slice_head(n = 20) %>%
  pivot_longer(cols      = c("Great cormorant", "Harbour seal"),
               names_to  = "predator",
               values_to = "FOO")

taxon_order_mf <- foo_mf_df %>%
  distinct(taxon_short, mean_foo) %>%
  arrange(mean_foo) %>%
  pull(taxon_short)

foo_mf_df$taxon_short <- factor(foo_mf_df$taxon_short, levels = taxon_order_mf)
foo_mf_df$predator    <- factor(foo_mf_df$predator,
                                levels = c("Harbour seal", "Great cormorant"))

p_mf <- ggplot(foo_mf_df,
               aes(x = FOO, y = taxon_short, colour = predator)) +
  geom_linerange(aes(xmin = 0, xmax = FOO),
                 position = position_dodge(width = 0.6),
                 linewidth = 0.5, alpha = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5, shape = 16) +
  scale_colour_manual(values = pred_colours, name = NULL,
                      guide  = guide_legend(override.aes = list(shape = 16,
                                                                linetype = 0))) +
  scale_x_continuous(limits = c(0, 55),
                     breaks = c(0, 10, 20, 30, 40, 50)) +
  labs(x = "Frequency of occurrence (%)", y = NULL,
       title = "A   MiFish — fish prey (top 20 by FOO)") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y        = element_text(face = "italic", size = 8),
    plot.title         = element_text(face = "bold", size = 10),
    legend.position    = "none",
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3)
  )

# ----------------------------------------------------------
# PANEL B: COI — top 20 taxa (fish + invertebrates mixed)
# ----------------------------------------------------------

tt_coi <- as.data.frame(tax_table(ps_coi_01))
tt_coi$LCA_taxon_clean <- gsub("_", " ", tt_coi$LCA_taxon)
tt_coi$is_invertebrate <- tt_coi$Class != "Actinopteri"

sd_coi  <- as.data.frame(sample_data(ps_coi_01))
otu_coi <- as.data.frame(otu_table(ps_coi_01))
if (taxa_are_rows(ps_coi_01)) otu_coi <- as.data.frame(t(otu_coi))
colnames(otu_coi) <- tt_coi$LCA_taxon_clean

otu_coi_agg <- t(rowsum(t(otu_coi), tt_coi$LCA_taxon_clean))

tax_agg <- tt_coi %>%
  group_by(LCA_taxon_clean) %>%
  summarise(is_invertebrate = any(is_invertebrate),
            Class           = first(Class),
            .groups         = "drop")

corm_rows_coi <- rownames(sd_coi)[sd_coi$species == "cormorant"]
seal_rows_coi <- rownames(sd_coi)[sd_coi$species == "harbor_seal"]
n_corm_coi    <- length(corm_rows_coi)
n_seal_coi    <- length(seal_rows_coi)

foo_corm_coi <- colSums(otu_coi_agg[rownames(otu_coi_agg) %in% corm_rows_coi, ] > 0) / n_corm_coi * 100
foo_seal_coi <- colSums(otu_coi_agg[rownames(otu_coi_agg) %in% seal_rows_coi, ] > 0) / n_seal_coi * 100

foo_coi_df <- data.frame(
  taxon             = colnames(otu_coi_agg),
  `Great cormorant` = as.numeric(foo_corm_coi),
  `Harbour seal`    = as.numeric(foo_seal_coi),
  check.names       = FALSE
) %>%
  left_join(tax_agg, by = c("taxon" = "LCA_taxon_clean")) %>%
  mutate(
    mean_foo    = (`Great cormorant` + `Harbour seal`) / 2,
    taxon_short = shorten_taxon(taxon),
    taxon_short = case_when(
      grepl("Evadne", taxon)                                                             ~ "Evadne sp.",
      grepl("Alitta|Hediste", taxon)                                                     ~ "Alitta / Hediste",
      grepl("Echinogammarus.*Gammarus", taxon)                                           ~ "Echinogammarus / Gammarus",
      grepl("Heteromastus|Capitellidae", taxon)                                          ~ "Heteromastus / Capitellidae",
      grepl("Palaemon", taxon)                                                           ~ "Palaemon adspersus",
      grepl("Harmothoe", taxon)                                                          ~ "Harmothoe imbricata",
      grepl("Liocarcinus", taxon)                                                        ~ "Liocarcinus holsatus",
      grepl("Penilia", taxon)                                                            ~ "Penilia avirostris",
      grepl("Acartia", taxon)                                                            ~ "Acartia spp.",
      grepl("Cyclopoida", taxon)                                                         ~ "Cyclopoida sp.",
      grepl("Clupea harengus.*Clupea pallasii|Clupea pallasii.*Clupea harengus", taxon) ~ "Clupea harengus",
      grepl("Pollachius virens.*environmental", taxon)                                   ~ "Pollachius virens",
      grepl("Melanogrammus", taxon)                                                      ~ "Melanogrammus aeglefinus",
      grepl("Micromesistius.*environmental", taxon)                                      ~ "Micromesistius poutassou",
      grepl("Enchelyopus.*environmental", taxon)                                         ~ "Enchelyopus cimbrius",
      grepl("Belone svetovidovi", taxon)                                                 ~ "Belone belone",
      grepl("Pleuronectes platessa.*Pleuronectes sp", taxon)                             ~ "Pleuronectes platessa",
      grepl("Pollachius pollachius.*Pollachius virens", taxon)                           ~ "Pollachius spp.",
      TRUE ~ taxon_short
    )
  ) %>%
  filter(mean_foo > 0) %>%
  group_by(taxon_short, is_invertebrate) %>%
  summarise(
    `Great cormorant` = max(`Great cormorant`),
    `Harbour seal`    = max(`Harbour seal`),
    mean_foo          = max(mean_foo),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_foo)) %>%
  slice_head(n = 20) %>%
  pivot_longer(cols      = c("Great cormorant", "Harbour seal"),
               names_to  = "predator",
               values_to = "FOO")

taxon_order_coi <- foo_coi_df %>%
  distinct(taxon_short, mean_foo) %>%
  arrange(mean_foo) %>%
  pull(taxon_short)

foo_coi_df$taxon_short    <- factor(foo_coi_df$taxon_short, levels = taxon_order_coi)
foo_coi_df$predator       <- factor(foo_coi_df$predator,
                                    levels = c("Harbour seal", "Great cormorant"))
foo_coi_df$is_invertebrate <- as.logical(foo_coi_df$is_invertebrate)

p_coi <- ggplot(foo_coi_df,
                aes(x = FOO, y = taxon_short,
                    colour = predator, shape = is_invertebrate)) +
  geom_linerange(aes(xmin = 0, xmax = FOO),
                 position = position_dodge(width = 0.6),
                 linewidth = 0.5, alpha = 0.5) +
  geom_point(position = position_dodge(width = 0.6), size = 2.5) +
  scale_colour_manual(values = pred_colours, name = NULL,
                      guide  = guide_legend(override.aes = list(shape  = 16,
                                                                linetype = 0))) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                     labels = c("FALSE" = "Fish",
                                "TRUE"  = "Invertebrate (secondary prey)"),
                     name   = "Prey type") +
  scale_x_continuous(limits = c(0, 55),
                     breaks = c(0, 10, 20, 30, 40, 50)) +
  labs(x = "Frequency of occurrence (%)", y = NULL,
       title = "B   COI Leray-XT — all prey (top 20 by FOO)") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y        = element_text(face = "italic", size = 8),
    plot.title         = element_text(face = "bold", size = 10),
    legend.position    = "right",
    axis.line.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3)
  )

# ----------------------------------------------------------
# COMBINE AND SAVE
# ----------------------------------------------------------

p_combined <- p_mf + p_coi +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position  = "right",
        legend.box       = "vertical",
        legend.spacing.y = unit(0.5, "cm"))

ggsave(file.path(VAL, "prey_species_detected_individual.png"),
       p_combined, width = 16, height = 8, dpi = 300, bg = "white")

cat("Saved: prey_species_detected.png\n")


# ---- Original MiFish functional groups and initial diet, chunk 6 ----
# ============================================================
# MIFISH FUNCTIONAL GROUP LOOKUP
# ============================================================

fg_mifish <- c(
  # Pelagic
  "Clupea harengus"                                                                         = "Pelagic",
  "Clupea harengus Clupea pallasii"                                                         = "Pelagic",
  "Clupea pallasii Sprattus sprattus Clupea harengus"                                       = "Pelagic",
  "Sprattus sprattus"                                                                        = "Pelagic",
  "Sprattus sprattus Clupea pallasii"                                                        = "Pelagic",
  "Sprattus sprattus Clupea harengus"                                                        = "Pelagic",
  "Clupea pallasii"                                                                          = "Pelagic",
  "Scomber scombrus"                                                                         = "Pelagic",
  "Belone belone"                                                                            = "Pelagic",
  "Gasterosteus aculeatus"                                                                   = "Pelagic",
  "Spinachia spinachia"                                                                      = "Pelagic",
  "Micromesistius poutassou"                                                                 = "Pelagic",
  "Hyperoplus lanceolatus"                                                                   = "Pelagic",
  "Ammodytes tobianus"                                                                       = "Pelagic",
  "Ammodytes marinus"                                                                        = "Pelagic",
  "Ammodytes marinus Hyperoplus lanceolatus Ammodytes tobianus"                              = "Pelagic",
  "Hyperoplus lanceolatus Ammodytes tobianus Ammodytes marinus Hyperoplus immaculatus"       = "Pelagic",
  "Argentina silus"                                                                          = "Pelagic",
  "Argentina"                                                                                = "Pelagic",
  "Salmo trutta"                                                                             = "Pelagic",
  "Mallotus villosus"                                                                        = "Pelagic",

  # Demersal
  "Gadus morhua"                                                                             = "Demersal",
  "Merlangius merlangus"                                                                     = "Demersal",
  "Merlangius merlangus Gadus morhua"                                                        = "Demersal",
  "Pollachius virens"                                                                        = "Demersal",
  "Pollachius pollachius"                                                                    = "Demersal",
  "Pollachius virens Pollachius pollachius"                                                  = "Demersal",
  "Pollachius pollachius Pollachius virens"                                                  = "Demersal",
  "Merluccius merluccius"                                                                    = "Demersal",
  "Molva molva"                                                                              = "Demersal",
  "Trisopterus minutus"                                                                      = "Demersal",
  "Enchelyopus cimbrius"                                                                     = "Demersal",
  "Raniceps raninus"                                                                         = "Demersal",
  "Aphia minuta"                                                                             = "Demersal",
  "Symphodus melops"                                                                         = "Demersal",
  "Labrus bergylta"                                                                          = "Demersal",
  "Ctenolabrus rupestris"                                                                    = "Demersal",
  "Gobiusculus flavescens"                                                                   = "Demersal",
  "Callionymus lyra"                                                                         = "Demersal",
  "Melanogrammus aeglefinus"                                                                 = "Demersal",

  # Benthic
  "Myoxocephalus scorpius"                                                                   = "Benthic",
  "Myoxocephalus scorpius Myoxocephalus quadricornis"                                        = "Benthic",
  "Myoxocephalus scorpius Myoxocephalus scorpioides Myoxocephalus aenaeus Myoxocephalus verrucosus" = "Benthic",
  "Myoxocephalus"                                                                            = "Benthic",
  "Taurulus bubalis"                                                                         = "Benthic",
  "Pomatoschistus minutus"                                                                   = "Benthic",
  "Pomatoschistus microps"                                                                   = "Benthic",
  "Gobius niger"                                                                             = "Benthic",
  "Pholis gunnellus"                                                                         = "Benthic",
  "Zoarces viviparus"                                                                        = "Benthic",
  "Zoarces viviparus Zoarces fedorovi Pholidapus dybowskii"                                 = "Benthic",
  "Anguilla anguilla"                                                                        = "Benthic",
  "Pleuronectes platessa"                                                                    = "Benthic",
  "Pleuronectes platessa Platichthys flesus"                                                 = "Benthic",
  "Pleuronectes platessa Pleuronectes sp. SH-2018"                                          = "Benthic",
  "Platichthys flesus"                                                                       = "Benthic",
  "Platichthys flesus Pleuronectes platessa Hippoglossoides platessoides"                   = "Benthic",
  "Hippoglossoides platessoides"                                                             = "Benthic",
  "Limanda limanda"                                                                          = "Benthic",
  "Limanda limanda Hippoglossus hippoglossus"                                                = "Benthic",
  "Scophthalmus rhombus"                                                                     = "Benthic",
  "Solea solea"                                                                              = "Benthic",
  "Glyptocephalus cynoglossus"                                                               = "Benthic",
  "Microstomus kitt"                                                                         = "Benthic",
  "Leptoclinus maculatus"                                                                    = "Benthic",
  "Lumpenus lampretaeformis"                                                                 = "Benthic",
  "Lumpenus lampretaeformis Lumpenus fabricii"                                               = "Benthic",
  "Lycodes vahlii Lycodes mcallisteri Lycodes terraenovae"                                   = "Benthic",
  "Lycodes vahlii Lycodes mcallisteri"                                                       = "Benthic"
)

# ============================================================
# DIET COMPOSITION — overall stacked bar per predator
# ============================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

tt_mf <- as.data.frame(tax_table(ps_01))
tt_mf$LCA_taxon_clean <- gsub("_", " ", tt_mf$LCA_taxon)
sd_mf  <- as.data.frame(sample_data(ps_01))
otu_mf <- as.data.frame(otu_table(ps_01))
if (taxa_are_rows(ps_01)) otu_mf <- as.data.frame(t(otu_mf))
colnames(otu_mf) <- tt_mf$LCA_taxon_clean

# Assign functional groups
tt_mf$fg <- fg_mifish[tt_mf$LCA_taxon_clean]
tt_mf$fg[is.na(tt_mf$fg)] <- "Other"
cat("Unassigned taxa:\n")
print(tt_mf$LCA_taxon_clean[tt_mf$fg == "Other"])

# Collapse Clupea harengus/pallasii amalgam
colnames(otu_mf) <- ifelse(
  grepl("Clupea harengus.*Clupea pallasii|Clupea pallasii.*Clupea harengus",
        colnames(otu_mf)),
  "Clupea harengus", colnames(otu_mf))
tt_mf$LCA_taxon_clean <- ifelse(
  grepl("Clupea harengus.*Clupea pallasii|Clupea pallasii.*Clupea harengus",
        tt_mf$LCA_taxon_clean),
  "Clupea harengus", tt_mf$LCA_taxon_clean)

# Aggregate duplicate columns after renaming
otu_mf <- t(rowsum(t(otu_mf), colnames(otu_mf)))

# Per-predator mean relative abundance
corm_rows <- rownames(sd_mf)[sd_mf$species == "cormorant"]
seal_rows  <- rownames(sd_mf)[sd_mf$species == "harbor_seal"]

mean_corm <- colMeans(otu_mf[rownames(otu_mf) %in% corm_rows, ])
mean_seal <- colMeans(otu_mf[rownames(otu_mf) %in% seal_rows, ])

# Rebuild fg lookup after aggregation
fg_lookup <- setNames(tt_mf$fg, tt_mf$LCA_taxon_clean)
fg_lookup <- fg_lookup[!duplicated(names(fg_lookup))]

# Top 15 taxa by mean across both predators
top15_names <- names(sort((mean_corm + mean_seal) / 2, decreasing = TRUE)[1:15])

# Build long data frame
build_df <- function(mean_vec, predator_label, top_names, fg_lkp) {
  top  <- mean_vec[top_names]
  other <- sum(mean_vec[!names(mean_vec) %in% top_names])
  df <- data.frame(
    taxon    = c(names(top), "Other"),
    relabund = c(as.numeric(top), other),
    predator = predator_label
  )
  df$fg <- c(fg_lkp[names(top)], "Other")
  df
}

df_corm <- build_df(mean_corm, "Great cormorant", top15_names, fg_lookup)
df_seal <- build_df(mean_seal, "Harbour seal",    top15_names, fg_lookup)
df_all  <- rbind(df_corm, df_seal)

# Shorten display names
df_all$taxon_short <- sapply(strsplit(df_all$taxon, " "), function(w) {
  if (length(w) >= 2) paste(w[1:2], collapse = " ") else w[1]
})
df_all$taxon_short[df_all$taxon == "Other"] <- "Other"

# Fix known amalgam shortenings
df_all$taxon_short <- case_when(
  grepl("Pollachius virens Pollachius", df_all$taxon) ~ "Pollachius spp.",
  grepl("Platichthys flesus Pleuronectes", df_all$taxon) ~ "Platichthys / Pleuronectes",
  grepl("Hyperoplus.*Ammodytes|Ammodytes.*Hyperoplus", df_all$taxon) ~ "Ammodytes / Hyperoplus",
  grepl("Limanda limanda Hippoglossus", df_all$taxon) ~ "Limanda limanda",
  grepl("Sprattus sprattus Clupea", df_all$taxon) ~ "Sprattus sprattus",
  TRUE ~ df_all$taxon_short
)

# Factor ordering — by mean relative abundance, Other at bottom
taxon_order <- df_all %>%
  group_by(taxon_short) %>%
  summarise(mean_ra = mean(relabund), .groups = "drop") %>%
  arrange(mean_ra) %>%
  pull(taxon_short)
taxon_order <- c("Other", taxon_order[taxon_order != "Other"])

df_all$taxon_short <- factor(df_all$taxon_short, levels = taxon_order)
df_all$predator    <- factor(df_all$predator,
                             levels = c("Great cormorant", "Harbour seal"))

# Colours — by functional group
fg_colours <- c(
  "Pelagic"  = "#4393c3",
  "Demersal" = "#8B4513",
  "Benthic"  = "#8B008B",
  "Other"    = "#D3D3D3"
)

# Map taxon to fg for fill
df_all$fg <- factor(df_all$fg, levels = names(fg_colours))

p_diet <- ggplot(df_all,
                 aes(x = predator, y = relabund,
                     fill = fg, group = taxon_short)) +
  geom_bar(stat = "identity", position = "stack",
           colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = fg_colours,
                    name   = "Functional group") +
  scale_y_continuous(limits = c(0, 105),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100")) +
  labs(x = NULL, y = "Mean relative abundance (%)",
       title = "Diet composition — MiFish 12S") +
  theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5),
    legend.position  = "right",
    axis.text.x      = element_text(size = 10),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
  )

# Add taxon labels inside bars
df_all_labels <- df_all %>%
  filter(relabund >= 1) %>%  # only label slices >= 1%
  group_by(predator) %>%
  arrange(predator, fg, desc(relabund)) %>%
  mutate(ypos = cumsum(relabund) - relabund / 2)

p_diet <- p_diet +
  geom_text(data = df_all_labels,
            aes(x = predator, y = ypos, label = taxon_short),
            inherit.aes = FALSE,
            size = 2.5, fontface = "italic", colour = "white")

ggsave(file.path(VAL, "intermediate_diet_composition_mifish_NOT_FINAL.png"),
       p_diet, width = 8, height = 7, dpi = 300, bg = "white")

cat("Saved: diet_composition_mifish.png\n")


# ---- MiFish renormalisation diagnostic, chunk 7 ----
# Fix 1: renormalise each sample to 100% after threshold
otu_mf_norm <- otu_mf / rowSums(otu_mf) * 100
otu_mf_norm[is.nan(otu_mf_norm)] <- 0  # handle zero-sum samples

# Recalculate means
mean_corm <- colMeans(otu_mf_norm[rownames(otu_mf_norm) %in% corm_rows, ])
mean_seal  <- colMeans(otu_mf_norm[rownames(otu_mf_norm) %in% seal_rows, ])

cat("Cormorant total after renorm:", round(sum(mean_corm), 1), "\n")
cat("Seal total after renorm:     ", round(sum(mean_seal), 1), "\n")

# Fix 2: add missing amalgam names to fg lookup
fg_mifish_extra <- c(
  "Myoxocephalus scorpius Myoxocephalus quadricornis"                                                    = "Benthic",
  "Myoxocephalus scorpius Myoxocephalus scorpioides Myoxocephalus aenaeus Myoxocephalus verrucosus"      = "Benthic",
  "Pleuronectes platessa Platichthys flesus"                                                             = "Benthic",
  "Platichthys flesus Pleuronectes platessa Hippoglossoides platessoides"                                = "Benthic",
  "Zoarces viviparus Zoarces fedorovi Pholidapus dybowskii"                                              = "Benthic",
  "Lumpenus lampretaeformis Lumpenus fabricii"                                                           = "Benthic",
  "Lycodes vahlii Lycodes mcallisteri"                                                                   = "Benthic",
  "Lycodes vahlii Lycodes mcallisteri Lycodes terraenovae"                                               = "Benthic",
  "Pleuronectes platessa Pleuronectes sp. SH-2018"                                                       = "Benthic",
  "Limanda limanda Hippoglossus hippoglossus"                                                             = "Benthic",
  "Pollachius virens Pollachius pollachius"                                                               = "Demersal",
  "Pollachius pollachius Pollachius virens"                                                               = "Demersal",
  "Merlangius merlangus Gadus morhua"                                                                    = "Demersal",
  "Ammodytes marinus Hyperoplus lanceolatus Ammodytes tobianus"                                          = "Pelagic",
  "Hyperoplus lanceolatus Ammodytes tobianus Ammodytes marinus Hyperoplus immaculatus"                   = "Pelagic",
  "Clupea pallasii Sprattus sprattus Clupea harengus"                                                    = "Pelagic",
  "Sprattus sprattus Clupea harengus"                                                                    = "Pelagic",
  "Sprattus sprattus Clupea pallasii"                                                                    = "Pelagic",
  "Belone svetovidovi Belone belone"                                                                     = "Pelagic",
  "Enchelyopus cimbrius Actinopterygii environmental sample"                                             = "Demersal",
  "Pollachius virens Pollachius environmental sample Actinopterygii environmental sample"                = "Demersal"
)

# Merge into main lookup
fg_mifish_full <- c(fg_mifish, fg_mifish_extra)

# Rebuild fg_lookup after aggregation
fg_lookup <- fg_mifish_full[colnames(otu_mf_norm)]
fg_lookup[is.na(fg_lookup)] <- "Other"

cat("\nStill unassigned after fix:\n")
print(names(fg_lookup)[fg_lookup == "Other"])

# Recalculate top15
top15_names <- names(sort((mean_corm + mean_seal) / 2, decreasing = TRUE)[1:15])

cat("\nTop 15 sums after renorm:\n")
cat("Cormorant:", round(sum(mean_corm[top15_names]), 1), "\n")
cat("Seal:     ", round(sum(mean_seal[top15_names]), 1), "\n")


# ---- Corrected MiFish final plot: rebuild data after renormalisation ----------
# The original sequential code recalculated mean_corm/mean_seal above but
# continued to plot df_all built BEFORE that calculation. Rebuild explicitly.
# This is a documented reproduction correction, not a new ecological method.

species_colours <- c(
  "Gasterosteus aculeatus" = "#9ecae1",
  "Sprattus sprattus" = "#4292c6",
  "Clupea harengus" = "#08519c",
  "Merlangius merlangus" = "#fdae6b",
  "Symphodus melops" = "#f16913",
  "Pollachius spp." = "#d94801",
  "Gadus morhua" = "#8c2d04",
  "Platichthys flesus" = "#d4b9da",
  "Taurulus bubalis" = "#c994c7",
  "Pomatoschistus minutus" = "#df65b0",
  "Pholis gunnellus" = "#e7298a",
  "Pleuronectes platessa" = "#980043",
  "Myoxocephalus scorpius" = "#49006a",
  "Other" = "#d9d9d9"
)

# Restrict separate coloured segments to the taxa actually shown in Tarjei's
# final figure; all remaining taxa are included in Other rather than dropped.
prey_shown <- c(
  "Clupea harengus", "Myoxocephalus scorpius",
  "Pleuronectes platessa", "Gadus morhua", "Symphodus melops",
  "Pollachius spp.", "Gasterosteus aculeatus", "Sprattus sprattus",
  "Platichthys flesus", "Merlangius merlangus"
)
stopifnot(length(prey_shown) == 10L,
          all(prey_shown %in% names(species_colours)))

short_mifish <- function(x) {
  x <- as.character(x)
  out <- vapply(strsplit(x, " ", fixed = TRUE), function(w) {
    if (length(w) >= 2L) paste(w[1:2], collapse = " ") else w[1]
  }, character(1))
  out[grepl("^Pollachius (virens|pollachius).*Pollachius", x)] <- "Pollachius spp."
  out
}

# Normalize each sample to 100%, as in Tarjei's preceding diagnostic block.
# Note: zero-prey samples remain zero, so predator means are rescaled after
# averaging to make each diet bar sum to 100% (conditional prey composition).
make_mifish_plot_rows <- function(mean_vec, predator_label) {
  if (anyNA(mean_vec) || any(mean_vec < 0)) stop("Invalid MiFish mean abundance")
  display <- short_mifish(names(mean_vec))
  by_species <- tapply(as.numeric(mean_vec), display, sum)
  vals <- setNames(numeric(length(prey_shown)), prey_shown)
  matched <- intersect(names(by_species), prey_shown)
  vals[matched] <- by_species[matched]
  other <- sum(by_species[!names(by_species) %in% prey_shown])
  total <- sum(vals) + other
  if (!is.finite(total) || total <= 0) stop("No MiFish prey for ", predator_label)
  values <- 100 * c(vals, Other = other) / total
  data.frame(predator = predator_label,
             taxon = names(values), relabund = unname(values),
             stringsAsFactors = FALSE)
}

df_all <- rbind(
  make_mifish_plot_rows(mean_corm, "Great cormorant"),
  make_mifish_plot_rows(mean_seal, "Harbour seal")
)
stopifnot(all(abs(tapply(df_all$relabund,df_all$predator,sum)-100)<1e-8))

# Original final figure's legend order (top to bottom). ggplot stacks
# the last factor level at the bottom, hence Other is last here.
taxon_order <- c(
  "Merlangius merlangus", "Platichthys flesus",
  "Sprattus sprattus", "Gasterosteus aculeatus",
  "Pollachius spp.", "Symphodus melops", "Gadus morhua",
  "Pleuronectes platessa", "Myoxocephalus scorpius",
  "Clupea harengus", "Other"
)
stopifnot(setequal(as.character(unique(df_all$taxon)), taxon_order),
          nrow(df_all) == 22L)
df_all$taxon <- factor(df_all$taxon, levels = taxon_order)
df_all$predator <- factor(df_all$predator,
                           levels=c("Great cormorant","Harbour seal"))

p_diet <- ggplot(df_all, aes(x=predator,y=relabund,fill=taxon)) +
  geom_col(colour="white",linewidth=0.3) +
  scale_fill_manual(values=species_colours, name="Prey species",drop=TRUE) +
  scale_y_continuous(limits=c(0,101),breaks=c(0,25,50,75,100),
                     expand=c(0,0)) +
  labs(x=NULL,y="Mean relative abundance (%)",
       title="Diet composition — MiFish 12S") +
  theme_classic(base_size=11) +
  theme(plot.title=element_text(face="bold",size=11,hjust=0.5),
        legend.position="right",legend.text=element_text(face="italic",size=9),
        legend.title=element_text(size=10),axis.text.x=element_text(size=11),
        panel.grid.major.y=element_line(colour="grey90",linewidth=0.3))

ggsave(file.path(VAL,"diet_composition_mifish_individual.png"),
       p_diet,width=9,height=7,dpi=300,bg="white")
write.csv(df_all,file.path(VAL,"05_MIFISH_FINAL_PLOT_DATA.csv"),row.names=FALSE)
cat("Saved corrected MiFish composition: each predator bar totals 100%\n")

# ---- Original COI diagnostic, chunk 9 ----
# ============================================================
# COI diet composition — diagnostic
# ============================================================

tt_coi <- as.data.frame(tax_table(ps_coi_01))
tt_coi$LCA_taxon_clean <- gsub("_", " ", tt_coi$LCA_taxon)
tt_coi$is_invertebrate <- tt_coi$Class != "Actinopteri"

sd_coi  <- as.data.frame(sample_data(ps_coi_01))
otu_coi <- as.data.frame(otu_table(ps_coi_01))
if (taxa_are_rows(ps_coi_01)) otu_coi <- as.data.frame(t(otu_coi))
colnames(otu_coi) <- tt_coi$LCA_taxon_clean

# Aggregate duplicate ASV names
otu_coi_agg <- t(rowsum(t(otu_coi), tt_coi$LCA_taxon_clean))

# Renormalise to 100% per sample
otu_coi_norm <- otu_coi_agg / rowSums(otu_coi_agg) * 100
otu_coi_norm[is.nan(otu_coi_norm)] <- 0

# Per-predator means
corm_rows_coi <- rownames(sd_coi)[sd_coi$species == "cormorant"]
seal_rows_coi <- rownames(sd_coi)[sd_coi$species == "harbor_seal"]

mean_corm_coi <- colMeans(otu_coi_norm[rownames(otu_coi_norm) %in% corm_rows_coi, ])
mean_seal_coi <- colMeans(otu_coi_norm[rownames(otu_coi_norm) %in% seal_rows_coi, ])

# Top 10
top10_coi <- names(sort((mean_corm_coi + mean_seal_coi) / 2, decreasing = TRUE)[1:10])

cat("=== COI top 10 taxa ===\n")
for (tx in top10_coi) {
  cat(sprintf("%-70s corm: %5.2f  seal: %5.2f\n",
              tx, mean_corm_coi[tx], mean_seal_coi[tx]))
}

cat("\nCormorant total:", round(sum(mean_corm_coi), 1), "\n")
cat("Seal total:     ", round(sum(mean_seal_coi), 1), "\n")

# Functional group for COI taxa
tax_agg_coi <- tt_coi %>%
  group_by(LCA_taxon_clean) %>%
  summarise(is_invertebrate = any(is_invertebrate),
            Class           = first(Class),
            .groups         = "drop")

cat("\n=== Class of top 10 ===\n")
for (tx in top10_coi) {
  cl <- tax_agg_coi$Class[tax_agg_coi$LCA_taxon_clean == tx]
  inv <- tax_agg_coi$is_invertebrate[tax_agg_coi$LCA_taxon_clean == tx]
  cat(sprintf("%-70s %s (invert: %s)\n", tx, cl, inv))
}


# ---- Original COI diet figure, chunk 10 ----
# ============================================================
# DIET COMPOSITION — COI
# ============================================================

# Functional group lookup for COI top 10
fg_coi_lookup <- c(
  "Clupea harengus"                                                                                        = "Pelagic",
  "Sprattus sprattus"                                                                                      = "Pelagic",
  "Micromesistius poutassou fish environmental sample"                                                     = "Pelagic",
  "Pollachius virens Pollachius environmental sample Actinopterygii environmental sample"                  = "Demersal",
  "Merlangius merlangus"                                                                                   = "Demersal",
  "Ctenolabrus rupestris"                                                                                  = "Demersal",
  "Pleuronectes platessa Pleuronectes sp. SH-2018"                                                        = "Benthic",
  "Myoxocephalus scorpius Myoxocephalus scorpioides Myoxocephalus aenaeus Myoxocephalus verrucosus"       = "Benthic",
  "Pomatoschistus minutus Pomatoschistus microps Neogobius melanostomus Gobius niger"                     = "Benthic",
  "Pontophilus norvegicus"                                                                                 = "Invertebrate"
)

# Shorten display names
shorten_coi <- function(x) {
  case_when(
    grepl("Pleuronectes platessa Pleuronectes sp", x)                                    ~ "Pleuronectes platessa",
    grepl("Myoxocephalus scorpius Myoxocephalus scorpioides", x)                         ~ "Myoxocephalus scorpius",
    grepl("Pollachius virens Pollachius environmental", x)                                ~ "Pollachius virens",
    grepl("Micromesistius poutassou fish", x)                                            ~ "Micromesistius poutassou",
    grepl("Pomatoschistus minutus Pomatoschistus microps", x)                            ~ "Pomatoschistus spp.",
    TRUE ~ sapply(strsplit(x, " "), function(w)
      if (length(w) >= 2) paste(w[1:2], collapse = " ") else w[1])
  )
}

# Build data frame
build_df_coi <- function(mean_vec, predator_label, top_names) {
  top   <- mean_vec[top_names]
  other <- sum(mean_vec[!names(mean_vec) %in% top_names], na.rm = TRUE)
  data.frame(
    taxon    = c(shorten_coi(names(top)), "Other"),
    relabund = c(as.numeric(top), other),
    predator = predator_label,
    fg       = c(as.character(fg_coi_lookup[names(top)]), "Other"),
    stringsAsFactors = FALSE
  )
}

df_corm_coi <- build_df_coi(mean_corm_coi, "Great cormorant", top10_coi)
df_seal_coi <- build_df_coi(mean_seal_coi, "Harbour seal",    top10_coi)
df_coi      <- rbind(df_corm_coi, df_seal_coi)

# Renormalise to 100%
df_coi <- df_coi %>%
  group_by(predator) %>%
  mutate(relabund = relabund / sum(relabund) * 100) %>%
  ungroup()

# Stack order — most abundant at bottom, Other at top
taxon_order_coi <- df_coi %>%
  group_by(taxon) %>%
  summarise(mean_ra = mean(relabund), .groups = "drop") %>%
  filter(taxon != "Other") %>%
  arrange(mean_ra) %>%
  pull(taxon)
taxon_order_coi <- c(taxon_order_coi, "Other")

df_coi$taxon    <- factor(df_coi$taxon,    levels = taxon_order_coi)
df_coi$predator <- factor(df_coi$predator,
                          levels = c("Great cormorant", "Harbour seal"))

# Colours — same family logic as MiFish
species_colours_coi <- c(
  # Pelagic — blues
  "Clupea harengus"        = "#08519c",
  "Sprattus sprattus"      = "#2171b5",
  "Micromesistius poutassou" = "#6baed6",

  # Demersal — oranges
  "Merlangius merlangus"   = "#fdae6b",
  "Pollachius virens"      = "#d94801",
  "Ctenolabrus rupestris"  = "#8c2d04",

  # Benthic — purples
  "Pleuronectes platessa"  = "#980043",
  "Myoxocephalus scorpius" = "#49006a",
  "Pomatoschistus spp."    = "#df65b0",

  # Invertebrate — greens
  "Pontophilus norvegicus" = "#238b45",

  # Other
  "Other"                  = "#d9d9d9"
)

p_diet_coi <- ggplot(df_coi,
                     aes(x = predator, y = relabund, fill = taxon)) +
  geom_bar(stat = "identity", position = "stack",
           colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = species_colours_coi,
                    name   = "Prey taxa",
                    guide  = guide_legend(reverse = FALSE)) +
  scale_y_continuous(limits = c(0, 101),
                     breaks = c(0, 25, 50, 75, 100),
                     labels = c("0", "25", "50", "75", "100"),
                     expand = c(0, 0)) +
  labs(x     = NULL,
       y     = "Mean relative abundance (%)",
       title = "Diet composition — COI Leray-XT") +
  theme_classic(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 11, hjust = 0.5),
    legend.position    = "right",
    legend.text        = element_text(face = "italic", size = 9),
    legend.title       = element_text(size = 10),
    axis.text.x        = element_text(size = 11),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
  )

ggsave(file.path(VAL, "diet_composition_coi_individual.png"),
       p_diet_coi, width = 9, height = 7, dpi = 300, bg = "white")

cat("Saved: diet_composition_coi.png (validation panel)\n")

# ---- Final manuscript figure: original four-panel A-D layout -----------------
# A: MiFish FOO; B: COI FOO; C: MiFish diet; D: COI diet.
# The underlying FOO and COI calculations remain unchanged.
# Keep individual plots in validation; only this combined image belongs to
# the final figure set.
pC_final <- p_diet + labs(title = "C   Diet composition — MiFish")
pD_final <- p_diet_coi + labs(title = "D   Diet composition — COI Leray-XT")

p_four <- (p_mf + p_coi) / (pC_final + pD_final) +
  patchwork::plot_layout(heights = c(1, 1))

final_four_path <- file.path(figs_dir, "diet_FOO_four_panels.png")
ggsave(final_four_path, plot = p_four,
       width = 18, height = 13, units = "in", dpi = 300, bg = "white")
if (!file.exists(final_four_path)) stop("Combined figure was not created")

# Remove stale previous manuscript figures; individual panels are recreated
# under validation above. These names are exclusive to module 04.
unlink(file.path(figs_dir, c(
  "prey_species_detected.png",
  "diet_composition_mifish.png",
  "diet_composition_coi.png"
)))
cat("Saved final four-panel figure: ", final_four_path, "\n", sep = "")


# ---- Original MiFish and COI invertebrate taxa lists, chunk 11 ----
# ============================================================
# FULL TAXA LISTS — MiFish and COI
# ============================================================

library(dplyr)

# ----------------------------------------------------------
# MIFISH TAXA TABLE
# ----------------------------------------------------------

tt_mf <- as.data.frame(tax_table(ps_01))
tt_mf$LCA_taxon_clean <- gsub("_", " ", tt_mf$LCA_taxon)

sd_mf  <- as.data.frame(sample_data(ps_01))
otu_mf <- as.data.frame(otu_table(ps_01))
if (taxa_are_rows(ps_01)) otu_mf <- as.data.frame(t(otu_mf))
colnames(otu_mf) <- tt_mf$LCA_taxon_clean

# Collapse Clupea amalgam
colnames(otu_mf) <- ifelse(
  grepl("Clupea harengus.*Clupea pallasii|Clupea pallasii.*Clupea harengus",
        colnames(otu_mf)),
  "Clupea harengus", colnames(otu_mf))
tt_mf$LCA_taxon_clean <- ifelse(
  grepl("Clupea harengus.*Clupea pallasii|Clupea pallasii.*Clupea harengus",
        tt_mf$LCA_taxon_clean),
  "Clupea harengus", tt_mf$LCA_taxon_clean)

# Aggregate duplicates
otu_mf <- t(rowsum(t(otu_mf), colnames(otu_mf)))

# Per-predator FOO
corm_rows <- rownames(sd_mf)[sd_mf$species == "cormorant"]
seal_rows  <- rownames(sd_mf)[sd_mf$species == "harbor_seal"]
n_corm <- length(corm_rows)
n_seal <- length(seal_rows)

foo_corm <- colSums(otu_mf[rownames(otu_mf) %in% corm_rows, ] > 0) / n_corm * 100
foo_seal <- colSums(otu_mf[rownames(otu_mf) %in% seal_rows, ] > 0) / n_seal * 100

# Functional group lookup
fg_lookup_clean <- fg_mifish_full
names(fg_lookup_clean) <- gsub("_", " ", names(fg_lookup_clean))

# Build MiFish table with safe matching
mifish_table <- data.frame(
  Taxon         = colnames(otu_mf),
  stringsAsFactors = FALSE
) %>%
  left_join(
    data.frame(
      Taxon            = tt_mf$LCA_taxon_clean,
      LCA_rank         = tt_mf$LCA_rank,
      stringsAsFactors = FALSE
    ) %>% distinct(Taxon, .keep_all = TRUE),
    by = "Taxon"
  ) %>%
  mutate(
    Functional_group = fg_mifish_full[Taxon],
    Functional_group = ifelse(is.na(Functional_group), "Other", Functional_group),
    FOO_cormorant    = round(as.numeric(foo_corm[Taxon]), 1),
    FOO_seal         = round(as.numeric(foo_seal[Taxon]), 1)
  ) %>%
  arrange(Functional_group, desc((FOO_cormorant + FOO_seal) / 2))

cat("MiFish taxa table rows:", nrow(mifish_table), "\n")
print(mifish_table)

write.csv(mifish_table,
          file.path(home_dir, "taxa_list_mifish.csv"),
          row.names = FALSE)
cat("Saved: taxa_list_mifish.csv\n")

# ----------------------------------------------------------
# COI INVERTEBRATE TAXA TABLE
# ----------------------------------------------------------

tt_coi <- as.data.frame(tax_table(ps_coi_01))
tt_coi$LCA_taxon_clean <- gsub("_", " ", tt_coi$LCA_taxon)
tt_coi$is_invertebrate <- tt_coi$Class != "Actinopteri"

sd_coi  <- as.data.frame(sample_data(ps_coi_01))
otu_coi <- as.data.frame(otu_table(ps_coi_01))
if (taxa_are_rows(ps_coi_01)) otu_coi <- as.data.frame(t(otu_coi))
colnames(otu_coi) <- tt_coi$LCA_taxon_clean

# Aggregate duplicates
otu_coi_agg <- t(rowsum(t(otu_coi), tt_coi$LCA_taxon_clean))

# Subset to invertebrates
invert_taxa <- tt_coi$LCA_taxon_clean[tt_coi$is_invertebrate]
invert_taxa <- unique(invert_taxa)
otu_inv <- otu_coi_agg[, colnames(otu_coi_agg) %in% invert_taxa, drop = FALSE]

# Per-predator FOO
corm_rows_coi <- rownames(sd_coi)[sd_coi$species == "cormorant"]
seal_rows_coi  <- rownames(sd_coi)[sd_coi$species == "harbor_seal"]
n_corm_coi <- length(corm_rows_coi)
n_seal_coi <- length(seal_rows_coi)

foo_corm_inv <- colSums(otu_inv[rownames(otu_inv) %in% corm_rows_coi, ] > 0) / n_corm_coi * 100
foo_seal_inv <- colSums(otu_inv[rownames(otu_inv) %in% seal_rows_coi, ] > 0) / n_seal_coi * 100

# Class lookup
class_lookup_inv <- setNames(
  tt_coi$Class[match(colnames(otu_inv), tt_coi$LCA_taxon_clean)],
  colnames(otu_inv)
)

# Build table
coi_table <- data.frame(
  Taxon         = colnames(otu_inv),
  Class         = class_lookup_inv[colnames(otu_inv)],
  LCA_rank      = tt_coi$LCA_rank[match(colnames(otu_inv), tt_coi$LCA_taxon_clean)],
  FOO_cormorant = round(as.numeric(foo_corm_inv), 1),
  FOO_seal      = round(as.numeric(foo_seal_inv), 1)
) %>%
  arrange(Class, desc((FOO_cormorant + FOO_seal) / 2))

cat("\nCOI invertebrate taxa table rows:", nrow(coi_table), "\n")
print(coi_table)

# Save as CSV
write.csv(coi_table,
          file.path(home_dir, "taxa_list_coi_invertebrates.csv"),
          row.names = FALSE)
cat("Saved: taxa_list_coi_invertebrates.csv\n")


# ---- Original full COI taxa list, chunk 12 ----
# ============================================================
# FULL COI TAXA LIST — all 172 taxa
# ============================================================

tt_coi_full <- as.data.frame(tax_table(ps_coi_01))
tt_coi_full$LCA_taxon_clean <- gsub("_", " ", tt_coi_full$LCA_taxon)

sd_coi  <- as.data.frame(sample_data(ps_coi_01))
otu_coi <- as.data.frame(otu_table(ps_coi_01))
if (taxa_are_rows(ps_coi_01)) otu_coi <- as.data.frame(t(otu_coi))
colnames(otu_coi) <- tt_coi_full$LCA_taxon_clean

# Aggregate duplicates
otu_coi_agg <- t(rowsum(t(otu_coi), tt_coi_full$LCA_taxon_clean))

# Per predator FOO
corm_rows_coi <- rownames(sd_coi)[sd_coi$species == "cormorant"]
seal_rows_coi  <- rownames(sd_coi)[sd_coi$species == "harbor_seal"]
n_corm_coi <- length(corm_rows_coi)
n_seal_coi <- length(seal_rows_coi)

foo_corm_full <- colSums(otu_coi_agg[rownames(otu_coi_agg) %in% corm_rows_coi, ] > 0) / n_corm_coi * 100
foo_seal_full <- colSums(otu_coi_agg[rownames(otu_coi_agg) %in% seal_rows_coi, ] > 0) / n_seal_coi * 100

# Class and rank lookup after aggregation
tax_agg_full <- tt_coi_full %>%
  group_by(LCA_taxon_clean) %>%
  summarise(
    Class    = first(Class),
    LCA_rank = first(LCA_rank),
    .groups  = "drop"
  )

# Build table
coi_full_table <- data.frame(
  Taxon         = colnames(otu_coi_agg),
  FOO_cormorant = round(as.numeric(foo_corm_full), 1),
  FOO_seal      = round(as.numeric(foo_seal_full), 1)
) %>%
  left_join(tax_agg_full, by = c("Taxon" = "LCA_taxon_clean")) %>%
  mutate(
    prey_type = ifelse(Class == "Actinopteri", "Primary prey (fish)", "Secondary prey (invertebrate)")
  ) %>%
  arrange(prey_type, Class, desc((FOO_cormorant + FOO_seal) / 2))

cat("Total taxa after aggregation:", nrow(coi_full_table), "\n")
cat("\nClass breakdown:\n")
print(table(coi_full_table$Class))

write.csv(coi_full_table,
          file.path(home_dir, "taxa_list_coi_full.csv"),
          row.names = FALSE)
cat("Saved: taxa_list_coi_full.csv\n")

# ---- Provenance and deliberately separated QA -------------------------------
# Note original output tables are in OUT/tables, not source home directory.
expected_figs <- file.path(figs_dir, 'diet_FOO_four_panels.png')
expected_tabs <- file.path(home_dir,c('taxa_list_mifish.csv',
 'taxa_list_coi_invertebrates.csv','taxa_list_coi_full.csv'))
if(!all(file.exists(expected_figs))) stop('Missing final figure(s): ',
 paste(basename(expected_figs[!file.exists(expected_figs)]),collapse=', '))
if(!all(file.exists(expected_tabs))) stop('Missing original taxa list(s): ',
 paste(basename(expected_tabs[!file.exists(expected_tabs)]),collapse=', '))
# Reproduction diagnostics only, never manuscript figures/tables.
output_manifest <- data.frame(
  file = c(basename(expected_figs), basename(expected_tabs)),
  kind = c(rep('figure', length(expected_figs)),
           rep('table', length(expected_tabs))),
  stringsAsFactors = FALSE
)
write.csv(output_manifest,
          file = file.path(VAL, '05_OUTPUT_MANIFEST.csv'),
          row.names = FALSE)
writeLines(c(
 'SOURCE: Tarjei prevalence.distribution.coi.mifish.Rmd (R chunks 4, 6-12).',
 'MiFish figure FOO: per-ASV FOO then max within shortened taxon names (not union FOO).',
 'COI figure FOO: sum ASVs into full LCA names then compute per-label FOO.',
 'MiFish diet: second original same-filename ggsave is final figure.',
 'Original first MiFish diet figure saved under validation only.',
 'Original MiFish plot used stale pre-normalisation df_all; this version rebuilds from normalised predator means and groups display taxa into the original labelled categories.',
 'Original COI taxon composition uses top10 by mean RRA (not all 172 displayed).',
 'MiFish final stacked plot corrects stale data flow; only ten original displayed taxa remain and all others are pooled as Other.',
 'Final figure combines A MiFish FOO, B COI FOO, C MiFish diet and D COI diet. Standalone previews are validation only.',
 'Different source object ps_01 could change exact output; visually verify against original.',
 'No extra manuscript outputs added; all QA in validation.'),
 file.path(VAL,'05_METHOD_NOTES.txt'))
capture.output(sessionInfo(),file=file.path(VAL,'05_SESSION_INFO.txt'))
cat('\n04 COMPLETE\nFigures: ',figs_dir,'\nTables: ',home_dir,
 '\nValidation: ',VAL,'\n',sep='')
