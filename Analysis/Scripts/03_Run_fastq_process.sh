#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=fastq_process
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$(jq '.email' config.json)

# Set up enviroment such that st_pipeline_run.py is in $PATH

# NOTE: Load miniconda module if you need the job to have access to a particular conda environment
# (uncomment if needed)
module load $(jq '.conda_environment_name' config.json)
conda init

conda activate $(jq '.conda_environment_name' config.json)

# NEED TO RUN OpticNerveHead_DBiT/Analysis/Python/fastq_process.py on *_[Read|R|]*2.fastq first
# FASTQ reads

FW_orig=`pwd`/../../data/$(jq '.fastq_FW' config.json)
python ../Python/fastq_process.py -i $FW_orig
