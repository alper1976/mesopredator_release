#!/bin/bash

# Job name: dada2_lst_Vg
#SBATCH --job-name=dada2_lst_Vg
#SBATCH --account=nn9745k
#SBATCH --time=10:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1

##Setup job environment
module --force purge   # clear any inherited modules
#set -o errexit # exit on errors

module load StdEnv           # Establish the standard environment path
#module load foss/2022a       # Load the toolchain
module load R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1
set -o errexit # exit on errors

SAMPEL_FOLDER="/cluster/projects/nn9744k/02_results/51_Stensrud_students/tnt/coi/AdaptersRemoved"
OUTPUT_FOLDER="/cluster/projects/nn9744k/02_results/51_Stensrud_students/tnt/coi/Figs"

Rscript /cluster/projects/nn9744k/scripts/51_Stensrud_students/tnt/coi/dada2_260414_parallel.R -f $SAMPEL_FOLDER -o $OUTPUT_FOLDER -m 12 -q 10 
