# PORTABLE STANDALONE RECONSTRUCTION. Original calculations retained.
# 02_Rarefaction_Threshold_Test.R
# Mesopredator_check_2 | FINAL reconstruction, module 01
# Original sources: Tarjei Rostum, updated.masterchunk.sample.quality.Rmd,
# ACE.alfa.beta.readthresholds.Rmd, percentage.test.Rmd.
# Reconstructed source: 05_read_threshold_diagnostics.R and
# 06_ACE_diversity_reproduction.R, plus validated manuscript threshold audit.
#
# From the PROJECT ROOT: source("00_final_reconstruction_scripts/02_Rarefaction_Threshold_Test.R")
# Outputs: only source figures and source threshold table in final_outputs/;
# reproducibility data, simulations and comparison checks in ../../docs/validation/.
# No legacy reconstruction scripts or stored derived objects are sourced.

if (!requireNamespace("phyloseq", quietly = TRUE) ||
    !requireNamespace("vegan", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("tidyverse", quietly = TRUE)) {
  stop("Install required packages: phyloseq, vegan, ggplot2, tidyverse")
}
suppressPackageStartupMessages({
  library(phyloseq)
  library(vegan)
  library(tidyverse)
})

# Resolve project root from either root or scripts folder.
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
input_candidates <- file.path(BASE, 'rawfiles', 'ps_decon_mifish.rds')
input_found <- input_candidates[file.exists(input_candidates)]
if (!length(input_found)) stop("Missing original ps_decon_mifish.rds; checked:\n",
                               paste(input_candidates, collapse = "\n"))
INPUT <- input_found[1L]
if (length(input_found)>1L) {
  hashes <- unname(tools::md5sum(input_found))
  if (length(unique(hashes))>1L) stop("Multiple nonidentical original ps_decon files; choose and audit input first")
}
OUT <- file.path(BASE, "processed_data", "02_Rarefaction_Threshold_Test")
FIG <- file.path(OUT, "figures")
TAB <- file.path(OUT, "tables")
VAL <- file.path(BASE, "validation", "02_Rarefaction_Threshold_Test")
for (d in c(FIG,TAB,VAL)) dir.create(d, recursive=TRUE, showWarnings=FALSE)

ps_decon <- readRDS(INPUT)
stopifnot(inherits(ps_decon,"phyloseq"))
os19 <- c("d-4-29","d-4-30","d-4-31","d-4-32")
stopifnot(all(os19 %in% sample_names(ps_decon)))
ps_decon <- prune_samples(!(sample_names(ps_decon) %in% os19), ps_decon)
metadata <- base::data.frame(sample_data(ps_decon), check.names=FALSE)
stopifnot(all(c("Sample_or_Control","species") %in% names(metadata)))
true_ids <- rownames(metadata)[!is.na(metadata$Sample_or_Control) &
                               as.character(metadata$Sample_or_Control)=="True_Sample"]
ps_decon <- prune_samples(true_ids, ps_decon)
ps_decon <- prune_taxa(taxa_sums(ps_decon)>0,ps_decon)
stopifnot(nsamples(ps_decon)==161L, ntaxa(ps_decon)==167L,
          sum(sample_sums(ps_decon))==580002)
cat("01 input: ",nsamples(ps_decon)," samples; ",ntaxa(ps_decon),
    " ASVs; ",sum(sample_sums(ps_decon))," reads\n",sep="")

# FIGURE 1A/B: exact two-panel plotting instructions from
# updated.masterchunk.sample.quality.Rmd (not a new visual redesign).
tab_rare <- as(otu_table(ps_decon),"matrix")
if (taxa_are_rows(ps_decon)) tab_rare <- t(tab_rare)
md_rare <- base::data.frame(sample_data(ps_decon),check.names=FALSE)
md_rare <- md_rare[rownames(tab_rare),,drop=FALSE]
sp <- gsub("top_cormorant","cormorant",as.character(md_rare$species),fixed=TRUE)
stopifnot(all(sp %in% c("cormorant","harbor_seal")))
colors_rare <- ifelse(sp=="cormorant","#d6604d","#2166ac")
n100 <- sum(rowSums(tab_rare)>=100)
n500 <- sum(rowSums(tab_rare)>=500)
stopifnot(n100==103L,n500==50L)
png(file.path(FIG,"rarefaction_curves.png"),width=12,height=6,units="in",res=300,bg="white")
tryCatch({
  par(mfrow=c(1,2),mar=c(5,4,2,1))
  rarecurve(tab_rare,step=10,lwd=1.5,ylab="ASVs",label=FALSE,
            col=colors_rare,xlab="Number of reads",xlim=c(0,2000),main="")
  abline(v=100,lty=2,lwd=1.5,col="black")
  abline(v=500,lty=3,lwd=1.5,col="black")
  text(110,11.5,"100 reads\n(n = 103)",adj=0,cex=0.75)
  text(510,9.5,"500 reads\n(n = 50)",adj=0,cex=0.75)
  mtext("a",side=3,adj=0,font=2,cex=1.1)
  rarecurve(tab_rare,step=20,lwd=1.5,ylab="",label=FALSE,
            col=colors_rare,xlab="Number of reads",main="")
  legend("bottomright",legend=c("Phalacrocorax carbo","Phoca vitulina"),
         col=c("#d6604d","#2166ac"),lwd=2,text.font=3,bty="n",cex=0.9)
  mtext("b",side=3,adj=0,font=2,cex=1.1)
}, finally=dev.off())

# PERCENT-ABUNDANCE THRESHOLD: global ASV retention by max percentage
# (this is NOT the per-sample zeroing variant explored in percentage.test.Rmd).
ps_percent <- transform_sample_counts(ps_decon,function(x) x/sum(x)*100)
thresholds <- c(0,0.1,0.5,1,2,5)
threshold_summary <- do.call(rbind,lapply(thresholds,function(t){
  z <- if(t==0) ps_percent else filter_taxa(ps_percent,function(x) max(x)>t,prune=TRUE)
  tx <- base::data.frame(tax_table(z),check.names=FALSE)
  stopifnot(all(c("LCA_taxon","LCA_rank") %in% names(tx)))
  lca <- as.character(tx$LCA_taxon)
  ranks <- tolower(trimws(as.character(tx$LCA_rank)))
  data.frame(Threshold=if(t==0) "None" else paste0(t,"%"),
    Samples=nsamples(z),ASVs=ntaxa(z),
    Unique_LCA_labels=length(unique(lca[!is.na(lca)])),
    Identified_LCA_taxa=length(unique(lca[!is.na(lca) &
      tolower(trimws(lca))!="no identification"])),
    Species_level_ASVs=sum(ranks=="species",na.rm=TRUE))
}))
stopifnot(identical(as.integer(threshold_summary$ASVs),c(167L,160L,145L,136L,130L,114L)),
  identical(as.integer(threshold_summary$Unique_LCA_labels),c(50L,50L,47L,46L,43L,40L)),
  identical(as.integer(threshold_summary$Species_level_ASVs),c(67L,67L,62L,58L,54L,45L)))
write.csv(threshold_summary,file.path(TAB,"threshold_summary.csv"),row.names=FALSE)
cat("\nVerified manuscript global threshold table:\n")
print(threshold_summary,row.names=FALSE)

# Read-depth diagnostics: original function below, copied from
# 06_ACE_diversity_reproduction.R (itself copied from Tarjei's ACE Rmd).
# Reference settings are 500 simulations, not a speed-optimized approximation.
depth_cutoff_diagnostics <- function(ps, cutoffs =
                                       c(100, 500, 1000, 2000, 3000, 4000, 5000,
                                         10000),
                                     run_simulations = FALSE,
                                     truth_min_depth = NULL,
                                     n_iter = 100,
                                     seed = 123,
                                     distance_method = "bray",
                                     distance_correlation_method = "spearman") {
  
  #Setting the seed inside the function to have reproducible data
  set.seed(seed)
  
  #Helper function to convert the phyloseq OTU-table into a matrix format, and
  #checking if it is formated correctly
  get_comm_matrix <- function(ps_object) {
    
    #Converts the OTU-table to a matrix
    mat <- as(otu_table(ps_object), "matrix")
    
    #Checks if the format is correct, with sample on rows
    if(taxa_are_rows(ps_object)) {
      
      #Transpose the matrix if it is in the wrong format
      mat <- t(mat)
    }
    
    #Returning the matrix in the correct format
    return(mat)
  }
  
  
  #Helper function to claculate the Gini coeficcient for abundance unevenness 
  gini_coefficient <- function(x) {
    
    #Removes all zero or negative values, so it only calulates unevenness among
    #categories that are actually present. E.g. species with abundance 0
    #is ignored
    x <- x[x > 0]
    
    #If after filtering there is only 0 or 1 positive values left, return NA,
    #since it is impossible to calculate the Gini coeficcient
    if(length(x) <= 1) return(NA_real_)
    
    #Sort the values from smallest to largest, since Gini assumes the
    #values are ordered
    x <- sort(x)
    
    #Finds the number of remaining positive values
    n <- length(x)
    
    #Calculates the gini coefficient. The seq_len(n) will make a vector of numbers
    #from 1 to n, where the smallest abundance gets multiplied by 1, the second
    #smallest abundance gets multiplied by 2, and so on, giving larger abundance
    #bigger weight in the calculation. The sum is the weighted abundance sum,
    #which is then multiplied by 2, to scale the weighted sum so that the result
    #end up on the same Gini coefficient scale.
    #The sum(x) is the total of the positive abundance values, and is scaled by
    #number of positive abundance values. This makes the calculation relative,
    #instead of depending on raw abundance size.
    #The last part is the correction term, which adjust the results so that
    #perfectly even abundances gives a Gini coefficient of 0. Essentially it
    #corrects for the ratio between the nominatior and denominator by removing
    #the ratio and scale the cofficiant to 0-1 
    gini <- (2 * sum(seq_len(n) * x) / (n * sum(x))) - ((n + 1) / n)
    return(gini)
  }
  
  
  
  #Calculate several alpha diversity metrics for each sample, in the
  #abundance matrix
  alpha_metrics <- function(mat) {
    
    #Counts how many taxa is present in each sample (ignoring 0 values). This
    #does not count abundance, but instances of taxa
    observed <- rowSums(mat > 0)
    
    #Caclulates the Shannon index. The Shannon index accounts for both richness
    #by how many taxa are present, in addition to eveness, which is how evenly
    #distributed the abundance is among taxa. A sample with many taxa and similar
    #abundance of the taxa get a high Shannon index
    shannon <- diversity(mat, index = "shannon")
    
    #Calculating the Gini-Simpson index. The value outputted represent the
    #probability that two randomly chosen individuals belong to different taxa.
    #Higher values, mean higher diversity, since it is less likely to randomly
    #choose two individuls of the same taxa
    simpson <- diversity(mat, index = "simpson")
    
    #The inverse Simpson index. It gives a value of "effective number of taxa". A
    #high inverse simpson value means that the sample behaves like it has more
    #equally abundant taxa. It is sensitive to dominant taxa, so if one taxa is
    #very dominant it will decrease the value
    invsimpson <- diversity(mat, index = "invsimpson")
    
    #Calculates the Pielou eveness, by dividing the Shannon diversity by the
    #maximum possible Shannon diversity for that number of observed taxa. (The
    #maximum possible Shannon diversity occurs when all observed taxa are
    #equally abundant). The value is usually between 0-1, where a value close
    #to 1 is very evenly distributed, and around 0 the sample is dominated by one
    #or few taxa
    pielou <- shannon / log(observed)
    
    #Sets the Pielou eveness to NA for samples with one or zero observed taxa,
    #because eveness is not meaningful if there is only one taxon present, since
    #there is nothing to compare it against
    pielou[observed <= 1] <- NA
    
    #Calculating the Abundance-based Coverage Estimator (ACE), which tries to
    #estimate the true richness by accounting for rare taxa. The ACE asks how many
    #taxa might be present, including taxa that may have been undersampled.
    #The apply() applies the function on the matrix across the rows (the 1)
    ace <- apply(mat, 1, function(x) {
      
      #Calculates the ACE value and extract only the ACE value, and hides
      #warnings since estimateR can produce warnings for samples with too few
      #taxa, too few rare taxa, or abundance patterns that make ACE difficult
      #to estimate. It returns one ACE value per sample 
      suppressWarnings(estimateR(x)["S.ACE"])
    })
    
    #Applies the ginni_coefficient function above to the samples by row (the 1  )
    gini <- apply(mat, 1, gini_coefficient)
    
    #creates a data-frame with all the information from the alpha diversity
    #metrics. These are all by sample, so the data-frame makes a sample column
    #with the sample names from the rows, and columns for all the metrics:
    #observed row, with how many taxa is present in the sample, the ACE value
    #for the sample, the Shannon index for the sample, the Gini_Simpson index for
    #the sample, the Inverse Simpson value for the sample, the Pielou eveness for
    #the sample, and the Gini abundance index for the sample. Lastly it removes
    #row names from the data-frame
    data.frame(
      Sample = rownames(mat),
      Observed = observed,
      ACE = as.numeric(ace),
      Shannon = shannon,
      Gini_Simpson = simpson,
      InvSimpson = invsimpson,
      Pielou_evenness = pielou,
      Gini_abundance = gini,
      row.names = NULL
    )
  }
  
  #Making a dataframe containing the sample names and the sequencing depth for
  #each sample 
  depth_df <- data.frame(
    Sample = sample_names(ps),
    Read_depth = sample_sums(ps)
  )
  
  
  #Finds the total taxa in the phyloseq object. This is used as a referance point
  #to calculate how many ASV's have been removed with each cutoff threshold
  total_asvs <- ntaxa(ps)
  
  #Making a data-frame by looping over every value of cutoffs to apply the
  #function which calculates the samples kept, samples removed, percent kept,
  #percent removed, ASV's kept, ASV's removed, percent of ASV's kept, percent of
  #ASV's removed, the mean ACE, and the median ACE for each cutoff threshold
  cutoff_df <- lapply(cutoffs, function(cutoff) {
    
    #Filtering out the samples which does not have enough sequencing depth, by
    #the looped threshold
    ps_cut <- prune_samples(sample_sums(ps) >= cutoff, ps)
    
    #Count the samples kept after filtering
    samples_kept <- nsamples(ps_cut)
    
    #Count how many samples was removed with the filtering
    samples_removed <- nsamples(ps) - samples_kept
    
    #If-test to check if any samples remains after the filtering, in case the
    #filtering is too strict, it will skip this step
    if(samples_kept > 0) {
      
      #Count how many ASV's are kept (ignoring ASV's with 0 abundance). This
      #calculates the prevalance of ASV's, not the abundance
      asvs_kept <- sum(taxa_sums(ps_cut) > 0)
      
      #Calculates how many ASV's are removed
      asvs_removed <- total_asvs - asvs_kept
      
      #Converts the phyloseq object to a matrix, with the correct format
      mat_cut <- get_comm_matrix(ps_cut)
      
      #Calculates the alpha metrics with the alpha metrics function above 
      #on the matrix
      alpha_cut <- alpha_metrics(mat_cut)
      
      #Calculates the mean ACE value by extracting the ACE from the alpha metrics
      #and using the mean() function on the ACE values
      mean_ACE <- mean(alpha_cut$ACE, na.rm = TRUE)
      
      #Calculates the median ACE value by extracting the ACE from the alpha metrics
      #and using the median() function on the ACE values
      median_ACE <- median(alpha_cut$ACE, na.rm = TRUE)
      
      #If all samples have been removed, the values will instead be default set
      #to either 0 or NA, and use the total_asvs as the ASV's removed
    } else {
      asvs_kept <- 0
      asvs_removed <- total_asvs
      mean_ACE <- NA
      median_ACE <- NA
    }
    
    #Making a data-frame with all the values from the calculations for one of the
    #looped cutoff thresholds. It only contains one row and the columns are: the
    #cutoff threshold, how many samples are kept, how many samples are removed,
    #the percent samples kept, the percent samples removed, how many ASV's are
    #kept, how many ASV's are removed, the percent of ASV's kept, the percent of
    #ASVS's removed, the mean ACE value for the samples, and the median ACE value
    #for the samples
    data.frame(
      Cutoff = cutoff,
      Samples_kept = samples_kept,
      Samples_removed = samples_removed,
      Percent_kept = round(100 * samples_kept / nsamples(ps), 1),
      Percent_removed = round(100 * samples_removed / nsamples(ps), 1),
      ASVs_kept = asvs_kept,
      ASVs_removed = asvs_removed,
      Percent_ASVs_kept = round(100 * asvs_kept / total_asvs, 1),
      Percent_ASVs_removed = round(100 * asvs_removed / total_asvs, 1),
      Mean_ACE = mean_ACE,
      Median_ACE = median_ACE)
    
    #When the small data-frame with a single row and all the metrics on the
    #columns is created, it gets piped with the |> which applies the data-frame
    #on the left hand side and applies the bind_rows() to the data-frame, which
    #combines the single row from this loop to all the loops before. Essentially
    #building the finished data-frame row by row for each cutoff threshold
  }) |>
    bind_rows()
  

  
  #Plotting the retention curve for the cutoff threshold values. The x-axis is
  #the different cutoff thresholds and the y-axis is how many samples remain
  p_samples <- ggplot(cutoff_df, aes(x = Cutoff, y = Samples_kept)) +
    geom_line() +
    geom_point(size = 3) +
    theme_bw() +
    labs(
      title = "Samples retained at different read-depth cutoffs",
      x = "Read-depth cutoff",
      y = "Number of samples retained"
    )
  
  #Plotting the retention curve for the cutoff threshold values with ASV's. The
  #x-axis is the different cutoff thresholds and the y-axis is how many ASV's
  #remain after cutoff
  p_asvs <- ggplot(cutoff_df, aes(x = Cutoff, y = ASVs_kept)) +
    geom_line() +
    geom_point(size = 3) +
    theme_bw() +
    labs(
      title = "ASVs retained at different read-depth cutoffs",
      x = "Read-depth cutoff",
      y = "Number of ASVs retained"
    )
  
  #Plotting the median ACE values for the different cutoff threshold values. The
  #x-axis is the different cutoff thresholds and the y-axis is the median ACE
  #values after cutoff thresholds
  p_ace <- ggplot(cutoff_df, aes(x = Cutoff, y = Median_ACE)) +
    geom_line() +
    geom_point(size = 3) +
    theme_bw() +
    labs(
      title = "Median ACE estimate at different read-depth cutoffs",
      x = "Read-depth cutoff",
      y = "Median ACE"
    )
  
  #Histogram of samples at different sequencing depths.The x-axis is the
  #sequencing depth and the y-axis is how many samples is at that sequencing depth
  p_depth_hist <- ggplot(depth_df, aes(x = Read_depth)) +
    geom_histogram(bins = 30, color = "black", fill = "skyblue") +
    theme_bw() +
    labs(
      title = "Read-depth distribution",
      x = "Reads per sample",
      y = "Number of samples"
    )
  

  #Placeholder values if run_simulations is set to false, or if fewer than 3
  #samples have read depth necessary after filtration
  alpha_sim_raw <- NULL
  alpha_sim_summary <- NULL
  beta_sim_raw <- NULL
  beta_sim_summary <- NULL
  
  #If-test to check if simulations is wanted
  if(run_simulations) {
    
    #Check to see if a minimum sequencing depth value has been provided
    if(is.null(truth_min_depth)) {
      
      #If a minimum sequencing depth value has not been provided. Stop the function
      stop("If run_simulations = TRUE, you must provide truth_min_depth.")
    }
    
    #Filter the phyloseq object by the minimum sequencing depth value, to only
    #keep samples which is now labeled as the true richness 
    ps_truth <- prune_samples(sample_sums(ps) >= truth_min_depth, ps)
    
    #If there are fewer than 3 samples remaining after filtering the samples
    if(nsamples(ps_truth) < 3) {
      
      #Give a warning that there is fewer than 3 samples remaining after filtering
      #for the threshold value, and therefore cannot do any simulations
      warning("Fewer than 3 samples have read depth >= truth_min_depth. Skipping simulations.")
    } else {
      
      #Converts the pholoseq object to a matrix format with the correct formating
      #of having samples on rows
      mat_truth <- get_comm_matrix(ps_truth)
      
      #Removes ASV's that are absent in the truth dataset. This is done by summing
      #the columns for an ASV and check if the value is larger than 0. If not,
      #the ASV is not present in any of the samples. drop = FALSE is there to
      #ensure that the matrix does not get converted to a vector
      mat_truth <- mat_truth[, colSums(mat_truth) > 0, drop = FALSE]
      
      #Calculates the alpha metrics with the above function and pipe the data-frame
      #to apply the pivot_longer() function. The pivot_longer converts the "wide"
      #format data-frame, where all the metrics have their own column, to a "long"
      #format data-frame, where all the metrics are in a single column. The Sample
      #column is excluded from this with the "-" to keep it as identifiers. The
      #names of the metrics get its own column, and the values another column
      truth_alpha <- alpha_metrics(mat_truth) |>
        pivot_longer(
          cols = -Sample,
          names_to = "Metric",
          values_to = "Truth_value"
        )
      
      #Make a subset of cutoff values which is smaller than the mimimum sequencing
      #depth of the samples, since you cannot rarefy samples to increase
      #sequencing depth, only decreasing sequencing depth
      valid_cutoffs <- cutoffs[
        
        #Only keeps cutoff thresholds which is larger than 0 AND (&) is smaller
        #than or equal to the minimum sequencing depth of the samples
        cutoffs > 0 &
          
          #This calculates the read depth for each of the samples, and finds the
          #smallest of them all. Then it checks with the cutoff threshold values
          #and only keeps threshold values which is less than or equal to the
          #minimum of the truth samples. This is important because all samples
          #needs to be able to be reduced to the chosen cutoff
          cutoffs <= min(rowSums(mat_truth))
      ]
      
     
      
      #Creating a variable to store a data-frame of the raw threshold data for all
      #iterations for all cutoff thresholds. Each of the threshold cutoff values
      #are looped over for the amount of iterations specified, to make a data-frame
      #containing all iterations for a threshold value, which is then combined into
      #a large data-frame that contains all iterations for all threshold values.
      #Each of the cutoff thresholds represents the read depth which the true
      #samples gets rarefied down to, to check how much alpha diversity is lost
      alpha_sim_raw <- lapply(valid_cutoffs, function(cutoff) {
        
        #Reiterates the rarefaction to the length of the number of iterations.
        #The reason to iterate the rarefaction, is because it is random and will
        #therefore get different ASV's each time. So the iterations is to average
        #the effect of the rarefaction
        lapply(seq_len(n_iter), function(i) {
          
          #Rarefy each sample in the matrix down to the looped cutoff threshold
          rare_mat <- rrarefy(mat_truth, sample = cutoff)
          
          #Calculating the alpha diversity on the rarefield matrix with the
          #alpha metrics function above. Then pipe the data-frame to the
          #The pivot_longer converts the "wide" format data-frame, where all the
          #metrics have their own column, to a "long" format data-frame, where all
          #the metrics are in a single column.
          #The Sample column is excluded from this with the "-" to keep it as
          #identifiers. The names of the metrics get its own column, and the
          #values another column
          alpha_i <- alpha_metrics(rare_mat) |>
            pivot_longer(
              cols = -Sample,
              names_to = "Metric",
              values_to = "Rarefied_value"
              
              #Pipes the data-frame to another function
            ) |>
            
            #Joins the alpha_i data-frame with the truth_alpha data-frame. This
            #is to have both values from before and after rarefying the data in
            #one data-frame. The joining is by Sample and Metric, so the function
            #joins the data-frames where both the sample name and the metric from
            #truth_alpha matches with the columns in the alpha_i. This makes an
            #additional column with the rarefied alpha metrics. Then pipes the
            #data-frame into mutate() to add additional columns
            left_join(truth_alpha, by = c("Sample", "Metric")) |>
            
            #Creates columns with mutate() for what the cutoff value is, how many
            #iterations was used to achieve the results, the difference between
            #the original alpha diversity values and the rarefied alpha diversity
            #values, in addition to the relative difference of change to the
            #true value. Giving a percent change
            mutate(
              Cutoff = cutoff,
              Iteration = i,
              Difference = Rarefied_value - Truth_value,
              Relative_difference = Difference / Truth_value)
          
          #Pipes the data-frame created from an iteration to the other iterations
        }) |>
          
          #Combines all the smaler data-frames from iterations into a large
          #data-frame containing all iterations
          bind_rows()
        
        #Piping the completed iterations data-frame for one cutoff threshold to
        #be combined with the other cutoff thresholds
      }) |>
        
        #Combining all cutoff threshold data-frames into a large data-frame
        #containing all iterations for all cutoff thresholds
        bind_rows()
      
      #Summarizing the raw alpha diversity data-frame by grouping each cutoff
      #threshold value and their metrics together. This is done by piping the
      #raw data-frame to the group_by() function
      alpha_sim_summary <- alpha_sim_raw |>
        
        #Grouping each of the cutoff threshold values to their metrics, to
        #independently look at each of the threshold at a time. Then piping
        #the grouped object to be summarized 
        group_by(Cutoff, Metric) |>
        
        #Summarize each of the metrics for each cutoff threshold value.
        #The metrics are: the mean original alpha diversity for the samples, the
        #mean rarefied alpha diversity for the samples, the mean difference
        #between the original and rarefied alpha diversity, the relative difference
        #between the original and rarefied alpha diversity, and the median
        #absolute (only positive) relative difference between the original and
        #rarefied alpha diversity. Then removing the grouping
        summarise(
          Mean_truth = mean(Truth_value, na.rm = TRUE),
          Mean_rarefied = mean(Rarefied_value, na.rm = TRUE),
          Mean_difference = mean(Difference, na.rm = TRUE),
          Mean_relative_difference = mean(Relative_difference, na.rm = TRUE),
          Median_absolute_relative_difference =
            median(abs(Relative_difference), na.rm = TRUE),
          .groups = "drop")
      
      
      
      #Converts the sample community matrix to relative abundance, so that each
      #sample (row) can be summed to 1. This is done by dividing all abundance
      #data (counts) by the total abundance in the sample, so each ASV becomes
      #relative to the sample instead of being raw abundance data
      truth_rel <- decostand(mat_truth, method = "total")
      
      #Calculates a distance matrix between samples, by using the relative
      #abundance data with the specified distance method in the outer function
      truth_dist <- vegdist(truth_rel, method = distance_method)
      
      #Converts the distance matrix into a vector, to compare distance matrices
      #before and after rarefaction
      truth_dist_vec <- as.vector(truth_dist)
      
      #Runs a NMDS plot on the distance matrix. The k is for how many dimensions
      #the plot should use, the trymax is how many times the plot is allowed to
      #restart to find a better solution, the autotransform is to not
      #automatically transform the data, and trace is to suppress messages
      truth_nmds <- metaMDS(
        truth_dist,
        k = 2,
        trymax = 100,
        autotransform = FALSE,
        trace = FALSE)
      
      
      #Creating a variable to store a data-frame of the raw threshold data for all
      #iterations for all cutoff thresholds. Each of the threshold cutoff values
      #are looped over for the amount of iterations specified, to make a data-frame
      #containing all iterations for a threshold value, which is then combined into
      #a large data-frame that contains all iterations for all threshold values.
      #Each of the cutoff thresholds represents the read depth which the true
      #samples gets rarefied down to, to check how much alpha diversity is lost
      beta_sim_raw <- lapply(valid_cutoffs, function(cutoff) {
        
        #Reiterates the rarefaction to the length of the number of iterations.
        #The reason to iterate the rarefaction, is because it is random and will
        #therefore get different ASV's each time. So the iterations is to average
        #the effect of the rarefaction
        lapply(seq_len(n_iter), function(i) {
          
          #Rarefy each sample in the matrix down to the looped cutoff threshold
          rare_mat <- rrarefy(mat_truth, sample = cutoff)
          
          #Converts the rarefied sample community matrix to relative abundance,
          #so that each sample (row) can be summed to 1. This is done by dividing
          #all abundance data (counts) by the total abundance in the sample,
          #so each ASV becomes relative to the sample instead of being raw
          #abundance data
          rare_rel <- decostand(rare_mat, method = "total")
          
          #Calculates a distance matrix between the rarefied samples, by using
          #the relative abundance data with the specified distance method in
          #the outer function
          rare_dist <- vegdist(rare_rel, method = distance_method)
          
          
          #Converts the rarefied distance matrix into a vector, to compare
          #distance matrices before and after rarefaction
          rare_dist_vec <- as.vector(rare_dist)
          
          #Correlating the vectorized distance matrix in this iteration for the
          #looped threshold with the original non-rarefied vectorized distance
          #matrix, to see how much the correlation change when the samples
          #are rarefied
          distance_correlation <- cor(truth_dist_vec, rare_dist_vec,
                                      method = distance_correlation_method)
          
          
          #Runs a NMDS plot on the distance matrix for the rarefied samples.
          #The k is for how many dimensions the plot should use, the trymax is
          #how many times the plot is allowed to restart to find a better
          #solution, the autotransform is to not automatically transform the
          #data, and trace is to suppress messages
          rare_nmds <- metaMDS(
            rare_dist,
            k = 2,
            trymax = 100,
            autotransform = FALSE,
            trace = FALSE)
          
          #Compares the two NMDS plots using Procrustes analysis, since NMDS
          #plots can be rotated, flipped, or scaled differently, while still
          #showing the same relationships. Essentially measuring how well
          #they match. The symmetric makes the test treat the two ordinations
          #more symetrically
          proc <- procrustes(truth_nmds, rare_nmds, symmetric = TRUE)
          
          #Stores the results of one cutoff threshold and one iteration in a
          #data-frame with columns for: the cutoff threshold that is looped, which
          #iteration the loop is in, the Procrustes sum of squares which is a
          #measure of mismatch between the two NMDS plots (original and rarefied),
          #Procrustes correlation, the stress of the original NMDS plot, and the
          #stress of the rarefied NMDS plot
          data.frame(
            Cutoff = cutoff,
            Iteration = i,
            Procrustes_ss = proc$ss,
            Procrustes_correlation = sqrt(1 - proc$ss),
            Distance_matrix_correlation = distance_correlation,
            Truth_stress = truth_nmds$stress,
            Rarefied_stress = rare_nmds$stress)
          
          #Pipes the data-frame created from an iteration to the other iterations
        }) |>
          
          #Combines all the smaler data-frames from iterations into a large
          #data-frame containing all iterations
          bind_rows()
        
        #Piping the completed iterations data-frame for one cutoff threshold to
        #be combined with the other cutoff thresholds
      }) |>
        
        #Combining all cutoff threshold data-frames into a large data-frame
        #containing all iterations for all cutoff thresholds
        bind_rows()
      
      
      #Summarizing the raw beta diversity data-frame by grouping each cutoff
      #threshold value and their metrics together. This is done by piping the
      #raw data-frame to the group_by() function
      beta_sim_summary <- beta_sim_raw |>
        
        #Grouping each of the cutoff threshold values to independently look at
        #each of the threshold at a time. Then piping the grouped object
        #to be summarized
        group_by(Cutoff) |>
        
        #Summarize the metrics for each of the cutoff thresholds. The metrics are:
        #the mean Procrustes correlation, the median Procrustes correlation, the
        #mean Procrustes sum of squares, the distance matrix correlation for the
        #lower quantile, the mean distance matrix correlation, the median distance
        #correlation, the distance matrix correlation for the upper quantile,
        #the original NMDS stress value, the mean rarefied NMDS plot stress,
        #and the mean stress difference.
        #Then ungroups the grouping
        summarise(
          Mean_Procrustes_correlation = mean(Procrustes_correlation, na.rm = TRUE),
          Median_Procrustes_correlation = median(Procrustes_correlation,
                                                 na.rm = TRUE),
          Mean_Procrustes_ss = mean(Procrustes_ss, na.rm = TRUE),
          
          Lower_distance_matrix_correlation =
            quantile(Distance_matrix_correlation, 0.25, na.rm = TRUE),
          Mean_distance_matrix_correlation =
            mean(Distance_matrix_correlation, na.rm = TRUE),
          Median_distance_matrix_correlation =
            median(Distance_matrix_correlation, na.rm = TRUE),
          Upper_distance_matrix_correlation =
            quantile(Distance_matrix_correlation, 0.75, na.rm = TRUE),
          
          Truth_stress = first(Truth_stress),
          Mean_rarefied_stress = mean(Rarefied_stress, na.rm = TRUE),
          Mean_stress_difference = mean(Rarefied_stress - Truth_stress,
                                        na.rm = TRUE),
          .groups = "drop")
    }
  }
  
  
  #Placeholder variables if the tests are not done
  p_alpha_sim <- NULL
  p_beta_sim <- NULL
  p_beta_dist_sim <- NULL
  
  #Checking if the alpha_sim_summary has been created in the eariler simulations
  #or if it is NULL, since no simulations has been done
  if(!is.null(alpha_sim_summary)) {
    
    #Making a plot of the alpha diversity simulations with the x-axis being the
    #cutoff threshold value, and the y-axis being the absolute median relative
    #difference between the original samples alpha diversity and the rarefied
    #samples alpha diversity
    p_alpha_sim <- ggplot(alpha_sim_summary,
                          aes(x = Cutoff, y = Median_absolute_relative_difference)  ) +
      geom_line() +
      geom_point(size = 3) +
      
      #Makes all other metrics get it's own plot, where the y-axis is free for
      #that specific plot
      facet_wrap(~Metric, scales = "free_y") +
      theme_bw() +
      labs(
        title = "Effect of subsampling depth on alpha diversity estimates",
        x = "Rarefaction depth",
        y = "Median absolute relative difference from high-depth reference")
  }
  
  
  #Checking if the beta_sim_summary has been created in the eariler simulations
  #or if it is NULL, since no simulations has been done
  if(!is.null(beta_sim_summary)) {
    
    #Making a plot of the beta diversity distance matrix correlation with the
    #x-axis being the cutoff threshold value, and the y-axis being the median
    #beta diversity distance matrix correlation
    p_beta_dist_sim <- ggplot(beta_sim_summary,
                              aes(x = Cutoff, 
                                  y = Median_distance_matrix_correlation)) +
      geom_line() +
      geom_point(size = 3) +
      theme_bw() +
      labs(
        title = "Preservation of Bray-Curtis distances across rarefaction depths",
        x = "Rarefaction depth",
        y = paste("Median distance-matrix correlation",
                  paste0("(", distance_correlation_method, ")")))
    
    #Making a plot of the beta diversity simulations with the x-axis being the
    #cutoff threshold value, and the y-axis being the median Procrustes correlation
    p_beta_sim <- ggplot(beta_sim_summary,
                         aes(x = Cutoff, y = Median_Procrustes_correlation)) +
      geom_line() +
      geom_point(size = 3) +
      theme_bw() +
      labs(
        title = "Conservation of beta-diversity pattern across rarefaction depths",
        x = "Rarefaction depth",
        y = "Median Procrustes correlation")
  }
  

  #Returning the information from the plot
  return(list(
    summary = list(
      cutoff_summary = cutoff_df,
      alpha_sim_summary = alpha_sim_summary,
      beta_sim_summary = beta_sim_summary),
    plots = list(
      samples_retained = p_samples,
      asvs_retained = p_asvs,
      ace_cutoff = p_ace,
      read_depth_histogram = p_depth_hist,
      alpha_simulation = p_alpha_sim,
      beta_simulation = p_beta_sim,
      beta_dist_simulation = p_beta_dist_sim),
    raw_data = list(
      alpha_sim_raw = alpha_sim_raw,
      beta_sim_raw = beta_sim_raw)
  ))
}

# Original simulation parameters, as explicitly written in ACE source.
a <- depth_cutoff_diagnostics(
  ps=ps_decon, cutoffs=c(100,500,1000,2000,3000,4000,5000),
  run_simulations=TRUE, truth_min_depth=5000, n_iter=500,
  seed=123, distance_method="bray",distance_correlation_method="spearman")
stopifnot(is.list(a),is.list(a$summary),
          a$summary$cutoff_summary$Samples_kept[a$summary$cutoff_summary$Cutoff==100]==103,
          a$summary$cutoff_summary$Samples_kept[a$summary$cutoff_summary$Cutoff==500]==50)
# Raw simulation records are audit information, not additional manuscript tables.
saveRDS(a,file.path(VAL,"ACE_original_diagnostics.rds"))
for (name in c("cutoff_summary","alpha_sim_summary","beta_sim_summary")) {
  x <- a$summary[[name]]
  if(!is.null(x)) write.csv(as.data.frame(x),file.path(VAL,paste0("ACE_",name,".csv")),row.names=FALSE)
}

# The final 3-panel figure in ACE.alfa.beta.readthresholds.Rmd references
# neutral_col, threshold_col, theme_methods and patchwork without defining
# them in that source file. Specify them here; geometry/data/annotations
# below otherwise reproduce the explicit plotting block in the original.
if (!requireNamespace("patchwork",quietly=TRUE)) stop("Install patchwork for original combined diagnostic figure")
neutral_col <- "#404040"
threshold_col <- "#d6604d"
theme_methods <- function() ggplot2::theme_classic(base_size=11)
library(patchwork)
# ── Panel A: Samples retained ─────────────────────────────────────────────────
p_samples <- ggplot(a$summary$cutoff_summary,
                    aes(x = Cutoff, y = Samples_kept)) +
  geom_line(color = neutral_col, linewidth = 0.8) +
  geom_point(color = neutral_col, size = 3) +
  geom_point(data = subset(a$summary$cutoff_summary, Cutoff == 100),
             aes(x = Cutoff, y = Samples_kept),
             color = threshold_col, size = 4) +
  geom_vline(xintercept = 100, linetype = "dashed",
             color = threshold_col, linewidth = 0.6) +
  geom_hline(yintercept = 50, linetype = "dotted",
             color = "grey60", linewidth = 0.6) +
  annotate("label", x = 2500, y = 44,
           label = "n = 50 (500 reads)",
           size = 3, color = "grey30", hjust = 0.5,
           fontface = "bold", label.size = 0, fill = "white") +
  annotate("label", x = 2500, y = 103,
           label = "n = 103 (100 reads)",
           size = 3, color = threshold_col, hjust = 0.5,
           fontface = "bold", label.size = 0, fill = "white") +
  scale_x_continuous(breaks = a$summary$cutoff_summary$Cutoff) +
  scale_y_continuous(limits = c(10, 120)) +
  labs(x = "Read-depth threshold", y = "Samples retained",
       title = "A  Sample retention",
       subtitle = "Dashed line = chosen threshold (100 reads)") +
  theme_methods()

# ── Panel B: Beta diversity ───────────────────────────────────────────────────
p_beta <- ggplot(a$summary$beta_sim_summary,
                 aes(x = Cutoff, y = Mean_Procrustes_correlation)) +
  geom_line(color = neutral_col, linewidth = 0.8) +
  geom_point(color = neutral_col, size = 3) +
  geom_point(data = subset(a$summary$beta_sim_summary, Cutoff == 100),
             aes(x = Cutoff, y = Mean_Procrustes_correlation),
             color = threshold_col, size = 4) +
  geom_vline(xintercept = 100, linetype = "dashed",
             color = threshold_col, linewidth = 0.6) +
  geom_hline(yintercept = 0.95, linetype = "dotted",
             color = "grey60", linewidth = 0.6) +
  annotate("label", x = 2500, y = 0.95,
           label = "r = 0.95",
           size = 3, color = "grey30", hjust = 0.5,
           fontface = "bold", label.size = 0, fill = "white") +
  scale_x_continuous(breaks = a$summary$beta_sim_summary$Cutoff) +
  scale_y_continuous(limits = c(0.4, 1.0)) +
  labs(x = "Read-depth threshold", y = "Procrustes correlation",
       title = "B  Beta diversity preservation",
       subtitle = "Mean Procrustes correlation across 500 iterations") +
  theme_methods()

# ── Panel C: Alpha diversity ──────────────────────────────────────────────────
alpha_observed <- subset(a$summary$alpha_sim_summary, Metric == "Observed")

p_alpha <- ggplot(alpha_observed, aes(x = Cutoff)) +
  geom_hline(yintercept = unique(alpha_observed$Mean_truth),
             color = "grey60", linewidth = 0.8, linetype = "dashed") +
  geom_line(aes(y = Mean_rarefied), color = neutral_col, linewidth = 0.8) +
  geom_point(aes(y = Mean_rarefied), color = neutral_col, size = 3) +
  geom_point(data = subset(alpha_observed, Cutoff == 100),
             aes(x = Cutoff, y = Mean_rarefied),
             color = threshold_col, size = 4) +
  geom_vline(xintercept = 100, linetype = "dashed",
             color = threshold_col, linewidth = 0.6) +
  annotate("label", x = 2500, y = unique(alpha_observed$Mean_truth),
           label = "True richness (reference)",
           size = 3, color = "grey40", hjust = 0.5,
           fontface = "bold", label.size = 0, fill = "white") +
  annotate("label", x = 2500, y = 2.93,
           label = "Rarefied richness",
           size = 3, color = neutral_col, hjust = 0.5,
           fontface = "bold", label.size = 0, fill = "white") +
  scale_x_continuous(breaks = alpha_observed$Cutoff) +
  scale_y_continuous(limits = c(2.3, 6.5)) +
  labs(x = "Read-depth threshold", y = "Mean observed ASV richness",
       title = "C  Alpha diversity recovery",
       subtitle = "Dashed grey = reference richness at 5,000 reads") +
  theme_methods()

# ── Combine and save ──────────────────────────────────────────────────────────
p_methods <- p_samples / p_beta / p_alpha


ggplot2::ggsave(file.path(FIG,"read_threshold_diagnostics.png"),
                plot=p_methods,width=7,height=10,dpi=300,bg="white")

write.csv(data.frame(metric=c("input","input_md5","samples","asvs","reads",
 "reads_ge_100","reads_ge_500","simulation_iterations","global_threshold_rule",
 "LCA_label_definition","diagnostic_palette_note"),
 value=c(INPUT,unname(tools::md5sum(INPUT)),nsamples(ps_decon),ntaxa(ps_decon),
 sum(sample_sums(ps_decon)),n100,n500,500,"max(sample_percentage)>cutoff",
 "includes no identification","ACE source omits definitions for diagnostic plot theme and colours; explicit replacements in this script")),
 file.path(VAL,"02_PROVENANCE.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(VAL,"02_SESSION_INFO.txt"))
cat("\n01 complete. Original-design figures:",FIG,
    "\nThreshold table:",file.path(TAB,"threshold_summary.csv"),
    "\nValidation:",VAL,"\n")
