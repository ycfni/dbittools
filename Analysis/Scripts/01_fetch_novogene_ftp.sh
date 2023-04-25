#!/bin/bash
#SBATCH --job-name=NovogeneFetch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --mail-user=$(jq '.email' config.json)
#SBATCH --mail-type=ALL
#SBATCH -o slurm-%j.out

HOST=$(jq '.host' credentials.json)
USER=$(jq '.host_user' credentials.json)
PASSWORD=$(jq '.host_password' credentials.json)

mkdir -p $(jq '.data_folder' config.json)/$(jq '.experiment_name' config.json)/raw/raw_molecular_data

wget -r --user=$USER --password=$PASSWORD $HOST -P $(jq '.data_folder' config.json)/$(jq '.experiment_name' config.json)/raw/raw_molecular_data/

