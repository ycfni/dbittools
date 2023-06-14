#!/bin/bash
#SBATCH --partition=week
#SBATCH --job-name=stpipeline1_50
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marcello.distasio@yale.edu

# Set up enviroment such that st_pipeline_run.py is in $PATH

# NOTE: Load miniconda module if you need the job to have access to a particular conda environment
# (uncomment if needed)
module load miniconda
conda init
conda activate $(jq -r '.conda_environment_name' ../config.json)


# NEED TO RUN OpticNerveHead_DBiT/Analysis/Python/fastq_process.py on *_[Read|R|]*2.fastq first
# FASTQ reads
RV=$(jq -r '.fastq_RV' ../config.json)
FW=$(jq -r '.fastq_FW_processed' ../config.json)

# References for mapping, annotation and nonRNA-filtering
MAP=`pwd`/../../ReferenceFiles/$(jq -r '.reference_files_folder_name' ../config.json)/
ANN=`pwd`/../../ReferenceFiles/genes.gtf

# Barcodes settings
ID=`pwd`/../$(jq -r '.barcodes_file' ../config.json)

# Output folder and experiment name
PATH_TO_OUTPUT=`pwd`/../../data/$(jq -r '.experiment_name' ../config.json)/$(jq -r '.output_relative_path' ../config.json)
OUTPUT=$PATH_TO_OUTPUT/st_pipeline_out/
mkdir -p $OUTPUT

PATH_TO_TEMP=`pwd`/../../data/$(jq -r '.experiment_name' ../config.json)/tmp/
TMP=$PATH_TO_TEMP/st_pipeline_tmp/
mkdir -p $TMP

# EXP is the experiment name Do not add / or \ to the experiment name
EXP=$(jq -r '.experiment_name' ../config.json)

# Running the pipeline
st_pipeline_run.py \
  --output-folder $OUTPUT \
  --ids $ID \
  --ref-map $MAP \
  --ref-annotation $ANN \
  --expName $EXP \
  --htseq-no-ambiguous \
  --verbose \
  --log-file $OUTPUT/${EXP}_log.txt \
  --temp-folder $TMP \
  --no-clean-up \
  --umi-start-position 16 \
  --umi-end-position 26 \
  --min-length-qual-trimming 10 \
  $FW $RV
