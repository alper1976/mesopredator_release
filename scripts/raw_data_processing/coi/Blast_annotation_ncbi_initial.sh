#!/bin/bash
# Job name: NCBI_annotation_mic
#SBATCH --job-name=NCBI_annotation_misc
#SBATCH --account=nn9744k
#SBATCH --time=10:30:00
#SBATCH --mem=4GB
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

### 
### seq, annotation COI, LCA + blast
# sequencing run
# 2025_09_02 sequencing run
# COI leray misc

Help()
{
   #Display Help
   echo "Unique Molecular Identifier based error correction for eDNA metabarcoding."
   echo " "
   echo "The script expects files to be in a fastq.gz format."
   echo "The script expects read1 and read2 files being seperated as R1 and R2 files respectively."
   echo " "
   echo "Arguments: [-h|-i|-o|-s|-r|-n]"
   echo " "
   echo "Options:"
   echo "-h    Display this help."
   echo "-i    Path to directory where rawdata file is stored (fastq.gz), char."
   echo "-o    Path to directory to store output data, char"
   echo "-s    Path to script"
   echo "-r    Path to the folder of reference database"
   echo "-n    Name of reference database in the path"
}


##Setup job environment
module purge   # clear any inherited modules
#set -o errexit # exit on errors


while getopts "hi:o:s:r:n:b:c:z:l:v:" option; do
   case $option in
      h) # display Help
         Help
         exit 0
         ;;
      i)
         Input_path="$OPTARG"
         ;;
      o)
         Output_path="$OPTARG"
         ;;
      s)
         Script_path="$OPTARG"
         ;;
      r)
         Reference_database_path="$OPTARG"
         ;;
      n) 
         Reference_database_name="$OPTARG"
         ;;
      b)
         Blast_id_cutoff="$OPTARG"
         ;;
      c)
         Blast_cov_cutoff="$OPTARG"
         ;;
      z)
         LCA_bitscore="$OPTARG"
         ;;
      l)
         LCA_id_cutoff="$OPTARG"
         ;;
      v)
         LCA_cov_cutoff="$OPTARG"
         ;;


   esac
done

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # set curret directory as working directory if non is selected.
#SCANDIFISH_DIR=/cluster/projects/nn9745k/03_databases/fish/ScandiFish_12s_v1.4/ # set reference database directory of ScandiFish if non is selected. # on saga

if [ -z "$Input_path" ]; then Input_path="/cluster/projects/nn9744k/02_results/45_Olli_Stensrud/COI/NCBI_all"; fi #where the rep-seqs.fna output file is
if [ -z "$Output_path" ]; then Output_path="/cluster/projects/nn9744k/02_results/45_Olli_Stensrud/COI/NCBI_all"; fi #Where to store the output, standard the same as $Input_path
if [ -z "$Script_path" ]; then Script_path="/cluster/projects/nn9744k/scripts/45_Olli_Stensrud/NCBI_annotation"; fi #if the script directory is not selected, script path is as standard for saga
if [ -z "$Reference_database_path" ]; then Reference_database_path="/cluster/projects/nn9745k/03_databases/COI/Invasive_species/241218"; fi #if the reference database is not selected, ScandiFish v1.4 will be used
if [ -z "$Reference_database_name" ]; then Reference_database_name="Invasive_Scandi_Coi_2024-12-18_db"; fi #if the reference database is not selected, ScandiFish v1.4 will be used
if [ -z "$Blast_id_cutoff" ]; then Blast_id_cutoff=99; fi #Blast species annotation identity cutoff #99
if [ -z "$Blast_cov_cutoff" ]; then Blast_cov_cutoff=99; fi #Blast species annotation coverage cutoff #99
if [ -z "$LCA_bitscore" ]; then LCA_bitscore=8; fi #LCA Bitscore top percentage threshold
if [ -z "$LCA_id_cutoff" ]; then LCA_id_cutoff=98; fi #LCA Bitscore top percentage threshold #98
if [ -z "$LCA_cov_cutoff" ]; then LCA_cov_cutoff=98; fi #LCA Bitscore top percentage threshold #98


#download reference database
#cd $Reference_database_path

#wget -nc https://raw.githubusercontent.com/EivindStensrud/ScandiFish/main/ScandiFish_12s_v1.4/ScandiFish_12s_v1.4_nf.fasta -P $Reference_database_path #Downloads database
#wget -nc https://raw.githubusercontent.com/EivindStensrud/ScandiFish/main/ScandiFish_12s_v1.4/ScandiFish_12s_v1.4_nf.fasta.md5 -P $Reference_database_path # Downloads md5sum

#md5sum -c $Reference_database_path/ScandiFish_12s_v1.4_nf.fasta.md5 > $Output_path/checkmd5.md5 # Creates files which can be checked if md5sum is correct. Needs to updated for new database version
# Downloads LCA script from GitHub.

printf "Initiation the annotation using BLAST with $Reference_database_path/$Reference_database_name \n"
echo " "
module purge

module load BLAST+/2.13.0-gompi-2022a

## Assigning taxonomy by BLASTN

blastn -query $Input_path/rep-seqs.fna -db /cluster/shared/databases/blast/2025-04-10/nt -task megablast -num_threads 16 -max_hsps 1 -perc_identity $LCA_id_cutoff -qcov_hsp_perc $LCA_cov_cutoff -out $Output_path/output_blast_results -outfmt "6 qseqid stitle sacc staxid pident qcovs evalue bitscore" -max_target_seqs 20 #works

