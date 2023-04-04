## OBJECTIVE
Produce a counts tsv output file from the st_pipeline (Spatial Transcriptomics Pipeline - https://github.com/SpatialTranscriptomicsResearch/st_pipeline) tool for downstream analysis of fastq files. 

## WHAT YOU'LL NEED
- Access to the folder at /gpfs/gibbs/pi/b_hafler/OpticNerveHead_DBiT, where we will use the scripts in /gpfs/gibbs/pi/b_hafler/OpticNerveHead_DBiT/Analysis/Scripts for the steps later described.
- As required by the above, to quote from the st_pipeline documentation (at the aforementioned GitHub link):
	- FASTQ files (Read 1 containing the spatial information and the UMI and read 2 containing the genomic sequence)
	- A genome index generated with STAR
An annotation file in GTF or GFF3 format (optional when using a transcriptome)
	- The file containing the barcodes and array coordinates (look at the folder "ids" to use it as a reference). Basically this file contains 3 columns (BARCODE, X and Y), so if you provide this file with barcodes identinfying cells (for example), the ST pipeline can be used for single cell data. This file is also optional if the data is not barcoded (for example RNA-Seq data). 

## STEPS
All steps are submittable as SBATCH jobs.
- Prepare your environment. Load conda environment ENVIRONMENT_NAME. In a terminal window, run:
	- conda activate ENVIRONMENT_NAME
- Get your data. Fetch from novogene (wget download). Use link/command provided by novogene data download email.
- Prepare your data - fastq process on the second read
	- python ../Python/fastq_process.py $FW_orig in script
- Run st_pipeline
- Run convertoname.sh to get final expression matrix from st_pipeline output

All you should have to edit are the config.json and credentials.json files. As a reference point, though, refer to the individual scripts, among the four scripts, to see per the comments where the file or folder paths may require replacement by you as a matter of a new run or a new user point of access.

In essence, you will be running four .sh scripts in sequence.
The four are to: get the data, make the STAR referencei (anticipate at least 20 minutes of runtime on the cluster once this job is submitted, potentially up to a few hours), run the fastq process to strip the barcodes for each read and prepend them, run st_pipeline to get the counts tsv file, and (not necessarily required by the newest versions of st_pipeline) run converttoname to account for gene names.

### FIRST, EDIT PATHS IN THE CONFIG FILE ON BEHALF OF EACH SCRIPT

To edit files, you can open them directly in the terminal using vim or emacs, or you can edit them externally in a preferred text editor. The Ruddle dashboard, for example, allows for files to be downloaded and re-uploaded if you want to download a config file to edit locally before updating the repository or directly putting the new version back on Ruddle.

In the crentials.json and config.json files, any field you wish to change can have the text directly replaced as with any text file, being sure to keep the quotation marks and commas surrounding each entry and to only alter the quoted text on the right-hand side of each colon. The format of the JSON file is to have the field name on the left, e.g. host_user being the variable name for your username, whereas what you would edit would be the username itself on the right-hand side of host_user. The end result is something like {... "host_user": "EDIT_YOUR_USERNAME_HERE", ...}

1. CREDENTIALS.JSON

Both credentials.json and config.json will need an email to reach you at to alert you of your job submission status once you run each of the four scripts as separate jobs on the cluster. Change "email" in each of these two JSON files to the email you'd like to receive job status updates.

Upon receiving an email from Novogene about your data being ready, which should be the case by now if you're here trying to use st_pipeline on some data, you should populate the credentials.JSON file. Replace the right-hand side of the "host" field with the link that Novogene provides and similarly populate "host_user" and "host_password" accordingly. Do not upload your an edited crendentials configuration file to Github with such fields populated, as a matter of standard security practice; the credentials file should be in the git ignore accordingly. The remainder of the changes you make to the other config.json file may be uploaded freely to this Github repository as the edited config.json.

Only the credentials.json file is required to run the first 01_fetch_novogene_ftp.sh script that will download our data. Once that data is downloaded, we move on to populate the config.json file with some fields being dependent upon this first script having already been run.

In summary, your credentials.json file should not be in the Github repository by nature of its content, but an example of the full file you can create on your own would be the following. Note that only the right-hand side's values are examples, whereas the left-hand side of each colon contains requisitely consistent keys you should not change the names of and should copy into your own file verbatim:

```
{
"email": "my_email@anydomain.example"
"host" : "ftp://usftp21.novogene.com",
"host_user" : "FILL_EMAILED_USER_ID",
"host_password" : "FILL_EMAILED_PASSWORD"
}
```

2. CONFIG.JSON

02_BuildSTARreference.sh will put together the gene reference files that st_pipeline requires. To the right side of each of these fields, we need to replace "reference_files_folder_name" and "genome_folder_path" if not using the prepopulated (per the original repository configuration file) defaults of hg38_STAR and /gpfs/ycga/datasets/genomes/Homo_sapiens/NCBI/GRCh38/ respectively.

Aside from describing the genome index files, fastq files, and other requirements & inner workings of st_pipeline, the source repository at https://github.com/jfnavarro/st_pipeline also describes needing to have STAR installed, for example. We handle installations through an Anaconda/Miniconda environment that we load up for certain scripts to run so that the environment may consistently offer certain versions of certain installations in order to make st_pipeline behavior consistent. The environment should already be available on Ruddle and not of your immediate concern to regenerate. This environment that we load by default is "stpipeline" in the "conda_environment_name" field. You should not need to change this environment, nor should you have any difficulties with loading this environment while running the 01 through 04 scripts sourced from the same top directory. 

Populate fastq_FW and fastq_RV with the forward and reverse fastq file paths you got from the Novogene download once you've run the 01_fetch_novogene_ftp.sh script to job completion. The script's Novogene download will have created a folder named after the link Novogene gave you and that contains a 01.RawData directory. That 01.RawData directory should contain a directory named after your sample, and that inner directory should contain your two fastq files. These files should end in .fq.gz per file and have either a 1 or a 2 before the file extension to indicate RV or FW respectively. We will use the fastq FW file for the 03_Run_fastq_process.sh job submission.

04_Run_stpipeline.sh, which finally runs st_pipeline, needs an output directory path for where the processed data is stored. The raw directory containing the Novogene data download, would have been in DATA_NAME/raw holding our fastq files. The st_pipeline output directory would correspondingly be the "processed" folder under your DATA_NAME folder as opposed to the "raw" folder in that DATA_NAME folder, e.g. would be "output_relative_path": "A22-1834_slide3/processed/" if "A22-1834_slide3" was the name of your data for processing. This fourth 04_Run_stpipeline.sh script will also need the spatial barcodes file, e.g. "barcodes_file": "Library/spatial_barcodes_index_1_50.txt", and an experiment name, e.g. "experiment_name": "ONH_A22-1834_S3".

### NOW RUN THE SCRIPTS IN SEQUENCE

Edit the credentials.json file with the Novogene email's indicated information for data retrieval, as described in the previous section. Open a Terminal window either in a Ruddle command-line interactive session or remote desktop, or ssh into Ruddle from a local Terminal/command-line window, and run (copy, then hit Enter) the following lines in sequence.

```
sbatch 01_fetch_novogene_ftp.sh 
```

Now edit the config.json file with updates to the data paths now that the data is downloaded from this first job's completion. Once the data paths and other config.json information as described in the prior section all get updated, run each othe following lines in order, paying attention to not start the next job (i.e. to enter/run the next line) until each submitted job has finished successfully. You can check the status of your submitted jobs (confirm that it's been submitted and is still running vs. crashed/finished, though you should receive the same alerts by email) by entering squeue --me in the terminal as described under the Common Slurm Commands section of this page: https://docs.ycrc.yale.edu/clusters-at-yale/job-scheduling/. So, in order and waiting for the completion of the job each line creates before you move on to running the next line, run each of the following lines individually in your terminal:
```
sbatch 02_BuildSTARreference.sh
sbatch 03_Run_fastq_process.sh
sbatch 04_Run_stpipeline.sh
```

Your final output should be a .tsv counts file, which we can use downstream for statistical and image analysis using read_dbit (a separate process/repository).
