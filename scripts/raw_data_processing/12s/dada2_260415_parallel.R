#!/usr/bin/env Rscript

#In R
#install.packages("/cluster/projects/nn9745k/Rpackages/BiocManager_1.30.10.tar.gz", repos = NULL, type="source")
#BiocManager::install("dada2", version = "3.10")

# load libraries
if (!require("stringr")) {
   install.packages("stringr", dependencies = TRUE)
   library(stringr)
   }
if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
   library(ggplot2)
   }
if (!require("dplyr")) {
   install.packages("dplyr", dependencies = TRUE)
   library(dplyr)
   }
if (!require("tibble")) {
   install.packages("tibble", dependencies = TRUE)
   library(tibble)
   }
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

if (!require("Biostrings", quietly = TRUE)){
  BiocManager::install("Biostrings", lib="/cluster/projects/nn9745k/Rpackages_4_2_1")
  library(Biostrings)
  }
if (!require("parallel")) {
   install.packages("parallel", dependencies = TRUE)
   library(parallel)
   }

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.15")

BiocManager::install(version = '3.15')
library("devtools")
if (!require("phyloseq")) {
   source("https://bioconductor.org/biocLite.R"); biocLite("phyloseq", dependencies = TRUE)
   library(phyloseq); packageVersion("phyloseq")
   }

if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
   library(ggplot2); packageVersion("ggplot2")
   }

if (!require("vegan")) {
   install.packages("vegan", dependencies = TRUE)
   library(vegan)
   }
if (!require("optparse")) {
   install.packages("optparse", dependencies = TRUE)
   library(optparse)
   }
library("dada2")

#### Set up flags
option_list = list(

    make_option(c("-t", "--truncate_length_fwd"),
        type="integer",
        default=0,
        help="Truncate length of forward reads [default %default]",
        metavar="number"),

    make_option(c("-l", "--truncate_length_rev"),
        type="integer",
        default=0,
        help="Truncate length of reverse reads [default %default]",
        metavar="number"),

    make_option(c("-q", "--truncate_quality"),
        type="integer",
        default=10,
        help="Truncates based on quality score [default %default]",
        metavar="character"),

    make_option(c("-f", "--fastq_dir"),
        type="character",
        default=" ",
        help="path to input",
        metavar="character"),

    make_option(c("-o", "--output_dir"),
        type="character",
        default=" ",
        help="path to output",
        metavar="character"),

    make_option(c("-s", "--database"),
        type="character",
        default=" ",
        help="use either silva or pr2 databases",
        metavar="character"),

    make_option(c("-n", "--max_n"),
        type="integer",
        default=0,
        help="maximum of ambiguous bases (Ns) [default %default]",
        metavar="number"),

    make_option(c("-m", "--min_overlap"),
        type="integer",
        default=10,
        help="minimum overlap of bases merging overlap [default %default]",
        metavar="number"),

    make_option(c("-E", "--max_EE_f"),
        type="integer",
        default=3,
        help="Maximum expected error R1/forward read [default %default]",
        metavar="number"),

    make_option(c("-R", "--max_EE_r"),
        type="integer",
        default=3,
        help="Maximum expected error R2/reverse read [default %default]",
        metavar="number"),
    
    make_option(c("-c", "--cores"),
        type="integer",
        default=10,
        help="number_of_cores [default %default]",
        metavar="number")

    )

opt = parse_args(OptionParser(option_list=option_list))

###############################################################################

### Number of cores to use

n_cores = opt$cores # the number of cores from slurm script 

#### Define the name of directories to use and create them if necessary
fastq_dir <- opt$fastq_dir  # fastq directory with the samples we are using
#fastq_dir <- "/cluster/projects/nn9745k/02_results/52_eivind202505/AdaptersRemoved"

output_dir <- opt$output_dir
#output_dir <- "/cluster/projects/nn9745k/02_results/52_eivind202505"
ifelse(!dir.exists(file.path(output_dir)), dir.create(file.path(output_dir)), FALSE)

database <- opt$database  # directory with database

filtered_dir <- file.path(output_dir, "Fastq_filtered")  # fastq filtered
ifelse(!dir.exists(file.path(filtered_dir)), dir.create(file.path(filtered_dir)), FALSE)

qual_dir <- file.path(output_dir, "Qual_pdf")  # qual pdf
ifelse(!dir.exists(file.path(qual_dir)), dir.create(file.path(qual_dir)), FALSE)

dada2_dir <- file.path(output_dir, "Dada2")  # dada2 results
ifelse(!dir.exists(file.path(dada2_dir)), dir.create(file.path(dada2_dir)), FALSE)

figs_dir <- file.path(output_dir, "Figs")
ifelse(!dir.exists(file.path(figs_dir)), dir.create(file.path(figs_dir)), FALSE)


#### Examine fastq files
# get a list of all fastq files in the fastq" directory and separate R1 and R2
fns <- sort(list.files(fastq_dir, full.names = TRUE))
fns <- fns[str_detect(basename(fns), ".fastq")]
fns_R1 <- fns[str_detect(basename(fns), "R1")]
fns_R2 <- fns[str_detect(basename(fns), "R2")]

sample.names = sapply(strsplit(basename(fns_R1), "_"), `[`, 1)

df <- data.frame()

for (i in 1:length(fns_R1)) {

  # use the dada2 function fastq.geometry
  geom <- fastq.geometry(fns_R1[i])

  # extract the information on number of sequences and file name
  df_one_row <- data.frame(n_seq = geom[1], file_name = basename(fns_R1[i]))

  # add one line to data frame
  df <- bind_rows(df, df_one_row)
}

#### Make some output tables and figures
write.table(df, file = file.path(figs_dir, "n_seq.txt"), sep='\t', row.names = FALSE, na='',
      quote=FALSE)



pdf(file.path(figs_dir, "quantity_profiles.pdf"))
  ggplot(df, aes(x = n_seq)) + geom_histogram(alpha = 0.5, position = "identity",
                                            binwidth = 1000)
  hist(df$n_seq)
dev.off()

pdf(file.path(figs_dir, "quality_profiles.pdf"))
  plotQualityProfile(fns_R1[25:34])
  plotQualityProfile(fns_R2[25:34])
dev.off()

write(paste0("The number of raw files is ",
        length(fns_R1) + length(fns_R2),
        " \n"),
        stdout())


#### set up filtering and trimming
filt_R1 <- file.path(filtered_dir, basename(fns_R1))
filt_R2 <- file.path(filtered_dir, basename(fns_R2))

filt.list = vector(mode = "list", length = length(fns_R1))
for (i in 1:length(fns_R1)){
  filt.list[[i]] = c(fns_R1[[i]], filt_R1[[i]], fns_R2[[i]], filt_R2[[i]])
}

#### Filter and Trim
parallel.filtering = function(input.list){
  fnFs.f = input.list[1]
  filtFs.f = input.list[2]
  fnRs.f = input.list[3]
  filtRs.f = input.list[4]
  filterAndTrim(fnFs.f, filtFs.f, fnRs.f, filtRs.f,
      truncLen = c(opt$truncate_length_fwd, opt$truncate_length_rev),
      #Quality profiles are good enough to keep full sequence lengths
      maxN=opt$max_n, maxEE=c(opt$max_EE_f,opt$max_EE_r), truncQ=opt$truncate_quality, rm.phix=TRUE,
      compress=TRUE, multithread=FALSE, matchIDs = TRUE)
}

out = mclapply(filt.list, parallel.filtering, mc.preschedule = FALSE, mc.set.seed = TRUE,
      mc.silent = FALSE, mc.cores = n_cores, mc.cleanup = TRUE)

# check in out if files are empty
# correct code
check_files = data.frame(cbind(unlist(out)[c(TRUE, FALSE)], unlist(out)[c(FALSE, TRUE)]))
if (nrow(check_files[check_files[,2] == 0,])) {
    write(paste0("Filtering did not result in any sequences from some files.\n",
        "The number of filtered sequence files is ",
        nrow(check_files[check_files[,2] == 0,]),
        " \n",
        "while the number of raw files is ",
        nrow(check_files),
        " \n"),
        stdout())
    }

filt <- sort(list.files(filtered_dir, full.names = TRUE))
filt <- filt[str_detect(basename(filt), ".fastq")]
filt_R1 <- filt[str_detect(basename(filt), "R1")]
filt_R2 <- filt[str_detect(basename(filt), "R2")]

filt_R1_seq = makeSequenceTable(as.list(filt_R1))

pdf(file.path(figs_dir, "length_profile_R1.pdf"))
  plot(table(nchar(getSequences(filt_R1_seq)))) #simple plot of length distribution
dev.off()

filt_R2_seq = makeSequenceTable(as.list(filt_R2))

pdf(file.path(figs_dir, "length_profile_R2.pdf"))
  plot(table(nchar(getSequences(filt_R2_seq)))) #simple plot of length distribution
dev.off()

save.image(file.path(dada2_dir, "dada2.Rdata"))

#### DADA2 error models
errF = mcparallel(learnErrors(filt_R1, multithread=TRUE, mc.affinity = n_cores, nbases = 1e8))$pid
errR = mcparallel(learnErrors(filt_R2, multithread=TRUE, mc.affinity = n_cores, nbases = 1e8))$pid
errF = mccollect(errF)[[1]]
errR = mccollect(errR)[[1]]

# plot error models
pdf(file.path(figs_dir, "error_plots.pdf"))
  plotErrors(errF, nominalQ=TRUE)
  plotErrors(errR, nominalQ=TRUE)
dev.off()

save.image(file.path(dada2_dir, "dada2.Rdata"))


#### Dereplication
derepFs = mclapply(filt_R1, derepFastq, mc.preschedule = FALSE, mc.set.seed = TRUE,
      mc.silent = FALSE, mc.cores = n_cores, mc.cleanup = TRUE)
derepRs = mclapply(filt_R2, derepFastq, mc.preschedule = FALSE, mc.set.seed = TRUE,
      mc.silent = FALSE, mc.cores = n_cores, mc.cleanup = TRUE)

names(derepFs) = basename(filt_R1)
names(derepRs) = basename(filt_R2)

sample.names = sapply(strsplit(basename(filt_R1), "_"), `[`, 1)

save.image(file.path(dada2_dir, "dada2.Rdata"))

#### Inference
parallel.inference.f = function(input){
  dada(input, err=errF, multithread=TRUE)
}
parallel.inference.r = function(input){
  dada(input, err=errR, multithread=TRUE)
}

dadaFs = mclapply(derepFs, parallel.inference.f, mc.preschedule = FALSE, mc.set.seed = TRUE,
  mc.silent = FALSE, mc.cores = n_cores, mc.cleanup = TRUE)
dadaRs = mclapply(derepRs, parallel.inference.r, mc.preschedule = FALSE, mc.set.seed = TRUE,
  mc.silent = FALSE, mc.cores = n_cores, mc.cleanup = TRUE)
dadaFs[[1]]

save.image(file.path(dada2_dir, "dada2.Rdata"))


#### Merging
merge.list = vector(mode = "list", length = length(dadaFs))
for (i in 1:length(dadaFs)){
  merge.list[[i]] = c(dadaFs[[i]], derepFs[[i]], dadaRs[[i]], derepRs[[i]])
}

parallel.merging = function(input.list){
  dadaFs.f = input.list[1]
  derepFs.f = input.list[2]
  dadaRs.f = input.list[3]
  derepRs.f = input.list[4]
  mergePairs(dadaFs.f, derepFs.f, dadaRs.f, derepRs.f, verbose=TRUE, minOverlap = opt$min_overlap)
}

mergers = mclapply(merge.list, parallel.merging, mc.preschedule = FALSE, mc.set.seed = TRUE,
  mc.silent = FALSE, mc.cores = n_cores, mc.cleanup = TRUE)


seqtab = makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))



#### Remove chimeras
seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
rownames(seqtab.nochim) = sample.names

#### Write ASV and taxonomy tables

saveRDS(seqtab.nochim, file.path(output_dir, "seqtab_nochim.rds"))

uniquesToFasta(seqtab.nochim, fout = file.path(output_dir, "rep-seqs.fna"), ids=colnames(seqtab.nochim))

#### Get length information on sequences
length_table = table(nchar(getSequences(seqtab.nochim)))
write.table(length_table, file.path(figs_dir, "length_table.tsv"))

pdf(file.path(figs_dir, "length_profile.pdf"))
  plot(table(nchar(getSequences(seqtab.nochim)))) #simple plot of length distribution
dev.off()




#### Transforming and saving the ASVs sequences
# seqtab.nochim_trans <- as.data.frame(t(seqtab.nochim)) %>% rownames_to_column(var = "sequence")
#     %>% rowid_to_column(var = "OTUNumber") %>% mutate(OTUNumber = sprintf("otu%04d",OTUNumber))
#     %>% mutate(sequence = str_replace_all(sequence, "# (-|\\.)", ""))
#
# seq_out <- Biostrings::DNAStringSet(seqtab.nochim_trans$sequence)
# names(seq_out) <- seqtab.nochim_trans$OTUNumber
#
# Biostrings::writeXStringSet(seq_out, file.path(dada2_dir, "ASV_no_taxo.fasta"),
#                             compress = TRUE, width = 20000)

## make training data from mothur output
# dada2:::makeTaxonomyFasta_Silva(file.path("/cluster/projects/nn9745k/03_databases/silva138/silva.nr_v138.align"),
#      file.path("/cluster/projects/nn9745k/03_databases/silva138/silva.nr_v138.tax"),
#      "/cluster/projects/nn9745k/03_databases/silva138/silva_nr_v138_train_set.fa.gz")



### Assign taxonomy

if (database == " ") {
  print("no database selected")

} else if (database == "silva") {

silva_file = "/cluster/projects/nn9745k/03_databases/silva138/silva_nr_v138_train_set.fa.gz"

taxa_silva = mcparallel(assignTaxonomy(seqtab.nochim,
    silva_file,
    minBoot=80,
    multithread=TRUE, tryRC = TRUE))$pid
taxa_silva = mccollect(taxa_silva)[[1]]

head(taxa_silva)
#taxa_silva[which.max(colSums(seqtab.nochim)),]

saveRDS(taxa_silva, file.path(output_dir, "silva_tax_final.rds"))

} else if (database == "pr2") {

pr2_file <- "/cluster/projects/nn8075k/03_databases/pr2/5.0.0/pr2_version_5.0.0_SSU_dada2.fasta"
PR2_tax_levels <- c("Kingdom", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family",
                    "Genus", "Species")

taxa_pr2 = mcparallel(assignTaxonomy(seqtab.nochim, refFasta = pr2_file, minBoot=80, taxLevels = PR2_tax_levels,
    multithread=TRUE, tryRC = TRUE))$pid
taxa_pr2 = mccollect(taxa_pr2)[[1]]

head(taxa_pr2)
#taxa_pr2[which.max(colSums(seqtab.nochim)),]

saveRDS(taxa_pr2, file.path(output_dir, "pr2_tax_final.rds"))
} else {
    print("this database is not implemented")
}

save.image(file.path(dada2_dir, "dada2.Rdata"))

#### Get summary information on sequences
getN = function(x) sum(getUniques(x))
track = cbind(unlist(out)[c(TRUE, FALSE)], unlist(out)[c(FALSE, TRUE)],
    sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names

write.table(track, file.path(figs_dir, "track_table.tsv"))

save.image(file.path(dada2_dir, "dada2.Rdata"))


write(paste0("Analysis run to the end."), stdout())
