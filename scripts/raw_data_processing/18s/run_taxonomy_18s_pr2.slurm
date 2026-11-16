#!/bin/bash
#SBATCH --job-name=tax_pr2
#SBATCH --account=nn9745k
#SBATCH --time=5:00:00      # Give it 10 hours just in case
#SBATCH --mem=96GB           # Massive memory for PR2 assignment
#SBATCH --cpus-per-task=16   # 16 cores for multithreading

module --force purge
module load StdEnv
module load foss/2022a
module load R-bundle-Bioconductor/3.15-foss-2022a-R-4.2.1

Rscript taxonomy_18s_pr2.R
