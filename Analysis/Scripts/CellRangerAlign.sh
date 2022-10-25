#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=mmd_align_201374
#SBATCH --cpus-per-task=12
#SBATCH --mem=100000
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=marcello.distasio@yale.edu

module load cellranger/6.0.1
cellranger count --id=mmd_align_AMD201374_HHT \
                 --transcriptome=$HOME/data/Retina_scRNAseq/ReferenceGenome/GRCh38-3.0.0_premrna \
                 --fastqs=$HOME/data/Retina_scRNAseq/AMD201374/ \
                 --sample=AMD201374_HHT \
                 --expect-cells=5000 \
                 --include-introns
