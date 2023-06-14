# ------------------------------------------------------------------------------------------------------------------------
# fastq_process.py
# Marcello DiStasio, modified from https://github.com/edicliuyang/DBiT-seq_FFPE/tree/master/Rawdata_processing
# Jan 2023

# Run as:
# python fastq_process.py -i path/to/input/file.fastq.gz
#
# ------------------------------------------------------------------------------------------------------------------------

from Bio import SeqIO
import gzip
import shutil
import os, sys
import argparse
 
# Initialize parser
parser = argparse.ArgumentParser()
# Adding optional argument
parser.add_argument("-i", "--input", help = "output file path")
args = parser.parse_args()

if args.input:
    input_filename = args.input
else:
    print("Input filename needed")
    sys.exit()

basedir = os.path.dirname(input_filename)
fname = os.path.basename(input_filename)

output_filename = os.path.join(basedir,''.join(fname.split('.')[:-2]) + '_processed.fastq')


# Process the file
output_handle = open(output_filename, "w")
with gzip.open(input_filename, "rt") as handle:
    for record in SeqIO.parse(handle, "fastq"):
      cut_record = record[32:40] + record[70:78] + record[22:32]  # BC2 + BC1 + UMI
      SeqIO.write(cut_record, output_handle, "fastq")
output_handle.close()


with open(output_filename, 'rb') as f_in:
    with gzip.open(output_filename+'.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

