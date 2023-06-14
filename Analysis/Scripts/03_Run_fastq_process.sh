#!/bin/bash
#SBATCH --partition=day
#SBATCH --job-name=fastq_process
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marcello.distasio@yale.edu

# Set up enviroment such that st_pipeline_run.py is in $PATH

# NOTE: Load miniconda module
module load miniconda
conda init
conda activate $(jq -r '.conda_environment_name' ../config.json)

# NEED TO RUN OpticNerveHead_DBiT/Analysis/Python/fastq_process.py on *_[Read|R|]*2.fastq first
# FASTQ reads

FW_orig=$(jq -r '.fastq_FW' ../config.json)
python ../Python/fastq_process.py -i $FW_orig
