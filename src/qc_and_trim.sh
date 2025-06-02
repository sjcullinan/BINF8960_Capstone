#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=capstone_trim_and_qc
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=2gb
#SBATCH --output=/scratch/sjc78466/capstone/job_logs/log.%j

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sjc78466@uga.edu

##########################################################################
###   LOADING ALL SOFTWARE MODULES NEEDED FOR TRIMMING AND QC CHECK    ###
##########################################################################

module load FastQC/0.11.9-Java-11
module load MultiQC/1.14-foss-2022a
module load Trimmomatic/0.39-Java-13

##############################
###   DEFINING VARIABLES   ###
##############################

# Variable for absolute path again so don't need to keep typing it out
workdir="/scratch/sjc78466/capstone"

# Variable for running Trimmomatic program
TRIMMOMATIC="java -jar /apps/eb/Trimmomatic/0.39-Java-13/trimmomatic-0.39.jar"
#############################################################
###   MAKING DIRECTORIES NEEDED FOR JOB TO RUN PROPERLY   ###
#############################################################

mkdir $workdir/data/trimmed_fastq
mkdir $workdir/results/trimmed_fastqc
mkdir $workdir/results/trimmed_multiqc

######################################################
### TRIM OFF ADAPTERS AND LOW QUALITY BASE CALLS   ###
######################################################

for fwd in $workdir/data/raw_fastq/*_1.fastq.gz
do
	sample=$(basename $fwd _1.fastq.gz)

	$TRIMMOMATIC PE $workdir/data/raw_fastq/${sample}_1.fastq.gz $workdir/data/raw_fastq/${sample}_2.fastq.gz  \
	$workdir/data/trimmed_fastq/${sample}_1.paired.fastq.gz $workdir/data/trimmed_fastq/${sample}_1.unpaired.fastq.gz \
	$workdir/data/trimmed_fastq/${sample}_2.paired.fastq.gz $workdir/data/trimmed_fastq/${sample}_2.unpaired.fastq.gz \
	ILLUMINACLIP:$workdir/data/raw_fastq/NexteraPE-PE.fa:2:30:10:5:True SLIDINGWINDOW:4:20
done

############################################################
###   RUNNING FASTQC ON TRIMMED READS TO CHECK QUALITY   ###
############################################################

fastqc -o $workdir/results/trimmed_fastqc $workdir/data/trimmed_fastq/*.paired.fastq.gz

####################################################################
###   RUNNING MULTIQC TO COMPILE QC RESULTS INTO SINGLE REPORT   ###
####################################################################

multiqc -o $workdir/results/trimmed_multiqc $workdir/results/trimmed_fastqc

# Moving MultiQC HTML file from cluster to local computer using FileZilla to visualize
