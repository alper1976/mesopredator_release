#!/

if (!require("ggplot2")) {
   install.packages("ggplot2", dependencies = TRUE)
   library(ggplot2); packageVersion("ggplot2")
   }
if (!require("xlsx")) {
   install.packages("xlsx", dependencies = TRUE)
   library(xlsx)
   }
if (!require("stringr")) {
   install.packages("stringr", dependencies = TRUE)
   library(stringr)
   }
if (!require("R.utils")) {
   install.packages("R.utils", dependencies = TRUE)
   library(R.utils)
   }
if (!require("reshape2")) {
   install.packages("reshape2", dependencies = TRUE)
   library(reshape2)
   }
if (!require("dplyr")){
  install.packages("dplyr")
  library(dplyr)
  }
if (!require("psych")) {
   install.packages("psych", dependencies = TRUE)
   library(psych)
   }
 if (!require("optparse")) {
   install.packages("optparse", dependencies = TRUE)
   library(optparse)
   }

option_list = list(
   make_option(c("-i", "--Input_path"), type="character", default=NULL, 
              help="path to fasta.otu.txt", metavar="character"),
   make_option(c("-o", "--Output_path"), type="character", default=NULL, 
              help="path to output", metavar="character")
  
)
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$Input_path)){
  print_help(opt_parser)
  stop("Input argument(s) missing.n", call=FALSE)
}

#Reformats the annotation to work with LCA annotation.

RootPath <- file.path(opt$Input_path) # Change to the path, where the result from DADA2 pipeline or similar pipelines is stored.
FigsPath <- file.path(opt$Output_path) # Path where output will be stored.


print("RootPath")
print(RootPath)

#Using the other file for this
#blast_result_file = (list.files(path = RootPath, pattern ="output_blast_results$",full.names=TRUE))

#blast_results_db = read.table(file=blast_result_file,sep="\t") # Results from BLAST analysis above.
#colnames(blast_results_db) <- c( "Query ID", "Subject title", "Subject accession","Subject Taxonomy ID", "Identity percentage", "Query Coverage Per Subject", "Evalue", "Bitscore" ) # NCBI annotation EchoPipe buildup

blast_result_tax_file = (list.files(path = RootPath, pattern ="taxonomy_output_covfiltered.tabular$",full.names=TRUE))
blast_results_tax_db = read.table(file=blast_result_tax_file,sep="\t", fill=TRUE) # Results from BLAST analysis above.
colnames(blast_results_tax_db) <- c( "Query_ID", "Subject_title", "Subject_accession","Subject_Taxonomy_ID", "Identity_percentage", "Query_Coverage_Per_Subject", "Evalue", "Bitscore", "Source", "Taxonomy" ) # NCBI annotation EchoPipe buildup
blast_results_tax_db = na.omit(blast_results_tax_db)

#create a dataframe with only taxa redundant
#tax_db = data.frame(stringr::str_replace_all(blast_results_tax_db$Taxonomy, ";", " / ")) # Replaces ";" with "/" in the taxonomic ranking to become compatible with downstream analysis, redundant
#colnames(tax_db) = c("Subject")
#taxonomy = sapply(strsplit(tax_db$Subject, split="\\|"), function(x) {
   # The taxonomy part is the third element in the split list
 #  taxonomy = x[3]
  
   # Remove potential trailing characters (like numbers and the last pipe) after the taxonomy
  # taxonomy_clean = sub("\\|.*$", "", taxonomy)
  
   # Return the clean taxonomy
   #return(taxonomy_clean)
#})

#blast_results_db$"Taxonomy" = taxonomy # creates a column with taxonomy
#blast_results_db$Subject = (stringr::str_replace_all(blast_results_db$Subject, ";", "/")) # Replaces ";" with "/" in the taxonomic ranking to become compatible with downstream analysis
#blast_results_db$Subject = (stringr::str_replace_all(blast_results_db$Subject, "[[|]]", "/")) # Replaces ";" with "/" in the taxonomic ranking to become compatible with downstream analysis


#tax_acc = blast_results_tax_db %>% select(Taxonomy, Subject_accession)
#tax_acc = data.frame((blast_results_tax_db$Taxonomy),(blast_results_tax_db$Subject accession))
#colnames(tax_acc)= c("Taxonomy")
#tax_accno = str_extract(tax_acc$Taxonomy, "(?<=gb\\/)[^/]+")

#tax_accno = data.frame(cbind(tax_acc, tax_accno))
#colnames(tax_accno) = c("Subject", "Subject accession") # Stores accession number and subject ID in same column.

#blast_results_db = left_join(blast_results_db, tax_accno, by = "Subject") # Joins blast_results by accession and subject ID.
#blast_results_db$"Subject Taxonomy ID" = blast_results_db$"Subject" # Data formating
#blast_results_db$"Source" = "EchoPipe"


#blast_results_tax_db_final = data.frame(cbind(blast_results_tax_db$"Query_ID",
#                                          blast_results_tax_db$"Subject_title",
#                                          blast_results_tax_db$"Subject_accession",
#                                          blast_results_tax_db$"Subject Taxonomy ID",
#                                          blast_results_tax_db$"Identity percentage",
#                                          blast_results_tax_db$"Coverage",
#                                          blast_results_tax_db$"Evalue",
#                                          blast_results_tax_db$"Bitscore",
#                                          blast_results_tax_db$"Source",
#                                          blast_results_tax_db$"Taxonomy"
#                                          ))

colnames(blast_results_tax_db) = c("#Query ID",
                                     "#Subject",
                                     "#Subject accession",
                                     "#Subject Taxonomy ID",
                                     "#Identity percentage",
                                     "#Coverage",
                                     "#Evalue",
                                     "#Bitscore",
                                     "#Source",
                                     "#Taxonomy"
                                          )



write.table(blast_results_tax_db,file.path(FigsPath, file="blast_result_file.lca.tabular"), row.names = FALSE, quote = FALSE, sep = "\t") # Stores the BLAST results in a table accepted by the LCA analysis.



