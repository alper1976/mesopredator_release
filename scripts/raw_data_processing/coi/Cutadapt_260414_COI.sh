#!/bin/bash
# Job name: Eiler Group Cutadapt
#SBATCH --job-name=Cutadapt_Eiler
#SBATCH --account=nn9745k
#SBATCH --time=3:30:00
#SBATCH --mem=12GB
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=10


######################################
### 


#change job name and account before use.

set -o errexit # exit on errors

### cutadapt
module load cutadapt/4.4-GCCcore-12.2.0


## Change Rawdir and Resdir location to where your samples are, and where you want them.

RAWDIR=/cluster/projects/nn9744k/01_raw_data/51_Stensrud_students/260409_SH01004.Project_Stensrud-amplicon5-09-04-26/tnt #location of rawdata files
RESDIR=/cluster/projects/nn9744k/02_results/51_Stensrud_students/tnt/coi #location of where you want the data to be stored.

# transfer files to scratch

mkdir -p $RESDIR/AdaptersRemoved
mkdir -p $RESDIR/Figs


##############################################
### run cutadapt to remove primer sequences
##############################################

cd $RAWDIR

#Change primers to your experiment, the default is the COI-Leray-XT primer set.
# Primer sequences for R1 and R2 (replace with your sequences, change overlap value as well)
# Note that cutadapt prefer working with N's instead of inosine bases, earlier version of Cutadapts do not tolerate Inosine
minimum_length=200
ADAPTER_SEQ_R1="GGWACWRGWTGRACWNTNTAYCCYCC;anywhere;o=18"
ADAPTER_SEQ_R2="TANACYTCNGGRTGNCCRAARAAYCA;anywhere;o=18"

#This code assumes this folder structure, and read structures.

#RAWDIR
#|
#|-- Sample_1_folder
#|          |-- Sample_1_1.fq.gz        
#|          |-- Sample_1_2.fq.gz
#|          |-- MD5.txt
#|
#|-- Sample_2_folder
#           |-- Sample_2_1.fq.gz
#           |-- Sample_2_2.fq.gz
#           |-- MD5.txt
#              

# if this is the folder structure, -mindepth 0 -maxdepth 1
#RAWDIR
#|
#|
#|-- Sample_1_1.fq.gz        
#|-- Sample_1_2.fq.gz
#|
#|
#|
#|-- Sample_2_1.fq.gz
#|-- Sample_2_2.fq.gz
#|
#              



# Find R1 files and infer R2 files, then process with cutadapt
## See comment above
find "$RAWDIR" -mindepth 0 -maxdepth 1 -type f -name "*R1_001.fastq.gz" | while read -r r1_file; do
    # Infer the R2 filename
    r2_file="${r1_file/R1_001.fastq.gz/R2_001.fastq.gz}"

    # Extract the filenames without the directory path
    r1_filename=$(basename "$r1_file")
    r2_filename=$(basename "$r2_file")

    # Define the output paths
    output_r1="$RESDIR/AdaptersRemoved/${r1_filename/.fastq.gz/_trimmed.fastq.gz}"
    output_r2="$RESDIR/AdaptersRemoved/${r2_filename/.fastq.gz/_trimmed.fastq.gz}"

    # Check if R2 file exists
    if [[ -f "$r2_file" ]]; then
        # Run cutadapt for paired-end reads
        echo "Processing $r1_file and $r2_file -> $output_r1, $output_r2"
        cutadapt -g "$ADAPTER_SEQ_R1" -G "$ADAPTER_SEQ_R2" -o "$output_r1" -p "$output_r2" --nextseq-trim=20 --minimum-length $minimum_length --discard-untrimmed "$r1_file" "$r2_file"
    else
        echo "Error: Matching R2 file not found for $r1_file"
    fi
done
