#!/bin/bash

# Job name: Eiler Group dada2
#SBATCH --job-name=Eiler_Group_dada2
#SBATCH --account=nn9745k
#SBATCH --time=5:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=12
#SBATCH --ntasks=1

#change job name and account before use.

##Setup job environment
module purge   # clear any inherited modules
set -o errexit # exit on errors

module load R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1

RESDIR=/cluster/projects/nn9744k/02_results/52_Stensrud_students/tnt

SAMPLE_FOLDER=$RESDIR/AdaptersRemoved
OUTPUT_FOLDER=$RESDIR/Figs

Rscript /cluster/projects/nn9744k/scripts/52_Stensrud_students/tnt/dada2_260415_parallel.R -f $SAMPLE_FOLDER -o $OUTPUT_FOLDER -m 10 -q 10 -E 1 -R 1 -c ${SLURM_CPUS_PER_TASK}


