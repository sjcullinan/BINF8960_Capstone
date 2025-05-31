#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=capstone_raw_qc
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=2gb
#SBATCH --output=/scratch/sjc78466/capstone/job_logs/log.%j

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sjc78466@uga.edu

##########################################################################
###   LOADING ALL SOFTWARE MODULES NEEDED FOR TRIMMING AND QC CHECK    ###
##########################################################################

module load FastQC/0.11.9-Java-11
module load MultiQC/1.14-foss-2022a

##############################
###   DEFINING VARIABLES   ###
##############################

# Variable for absolute path again so don't need to keep typing it out
workdir="/scratch/sjc78466/capstone"

##################################################################
###   MAKING OUTPUT DIRECTORY FOR FASTQC AND MULTIQC RESULTS   ###
##################################################################

mkdir $workdir/results/raw_fastqc
mkdir $workdir/results/raw_multiqc

#########################################################
###   RUNNING FASTQC ON RAW READS TO CHECK QUALITY   ###
#########################################################

fastqc -o $workdir/results/raw_fastqc $workdir/data/raw_fastq/*.fastq.gz

####################################################################
###   RUNNING MULTIQC TO COMPILE QC RESULTS INTO SINGLE REPORT   ###
####################################################################

multiqc -o $workdir/results/raw_multiqc $workdir/results/raw_fastqc

# Moving MultiQC HTML file from cluster to local computer using FileZilla to visualize
