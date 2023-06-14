#!/bin/bash
#SBATCH --job-name=NovogeneFetch
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --mail-user=marcello.distasio@yale.edu
#SBATCH --mail-type=ALL
#SBATCH -o slurm-%j.out


HOST="ftp://usftp21.novogene.com"
USER="X202SC23043586-Z01-F001"
PASSWORD="13exnz0d"

DATADIR=/home/mmd47/labshare/BrainEpendyma_DBiT/data/raw/raw_molecular_data/
mkdir -p $DATADIR

wget -r --user=$USER --password=$PASSWORD $HOST -P $DATADIR

