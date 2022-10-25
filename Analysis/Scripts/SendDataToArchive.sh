#!/bin/bash

################################################################################
# Marcello DiStasio
# December, 2021
#
# Usage:
#
# . SendDataToArchive.sh /path/to/folder/
#
################################################################################


# Path to backup
src=$(readlink -f $(pwd)/../../data)

# Where to backup to.
dest="/SAY/standard/ycfni01-CC1410-MEDPAT"

# Create archive filename. This uses the TWO tail directories to set the filename (e.g. .../MyProject/data --> MyProject-data)
IFS="/" read -ra PARTS <<< "$src"
PARTSlen=${#PARTS[@]}
outname=${PARTS[$PARTSlen-2]}-${PARTS[$PARTSlen-1]}

TIMESTAMP=`date +%Y-%m-%d_%H-%M-%S`
#outname=${src%\/}
#outname=${outname//\//-}
archive_file="${outname}__${TIMESTAMP}_${USER}.tar.gz"

# Log file
LOG_FILE="${archive_file}_ARCHIVAL.log"

# Print start status message.
echo "Backing up $src to $dest"

# Log the command used to archive
TIMESTAMP=`date +%Y-%m-%d_%H-%M-%S`
command="tar cvzf "${dest}/${archive_file}" $src"
echo "$TIMESTAMP:: $command" >> $LOG_FILE




# Execute the command
eval $command



echo "$TIMESTAMP:: Done. " >> $LOG_FILE


# DONE
echo "Created archive file: $dest/$archive_file."
echo "Log file at $LOG_FILE"

