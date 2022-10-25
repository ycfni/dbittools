#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=RetrieveData
#SBATCH --cpus-per-task=12
#SBATCH --mem=100000
#SBATCH --time=72:00:00
#SBATCH --mail-type=ALL
#SBATCH -o tmp/slurm-%j.out


################################################################################
# Marcello DiStasio
# December, 2021
#
# Usage:
#
# . GetDataFromArchive.sh 
#
# OR
#
# sbatch GetDataFromArchive.sh 
#
################################################################################



### ----->>>> EDIT TO POINT TO THE ARCHIVED DATA FOR YOUR PROJECT, IF IT EXISTS <<<<---------

# Path to archived data for this project (e.g. src="/SAY/standard/ycfni01-CC1410-MEDPAT/Retina_snRNAseq_01-data__2022-01-18_13-46-47_mmd47.tar.gz")

src="/SAY/standard/ycfni01-CC1410-MEDPAT/ARCHIVE_FILENAME.tar.gz"

#----------------------------------------------------------------------------------------------------




# Where to backup to.
dest=$(pwd)/../..

set -e
if test -d ${dest}/data; then
    echo "${dest}/data exists already. Aborting retrieval"
    kill -INT $$
fi


# Log file
TIMESTAMP=`date +%Y-%m-%d_%H-%M-%S`
mkdir -p log
LOG_FILE="log/ARCHIVE_RETRIEVAL_${TIMESTAMP}.log"


# Print start status message.
echo "Retrieving archive $src to $dest"

# Log the command used to archive
TIMESTAMP=`date +%Y-%m-%d_%H-%M-%S`
command="tar xvzf $src --directory $dest"
echo "$TIMESTAMP:: $command" >> $LOG_FILE




# Execute the command
eval $command



echo "$TIMESTAMP:: Done. " >> $LOG_FILE


# DONE
echo "Extracted archive file to: $dest/data"
echo "Log file at $LOG_FILE"

