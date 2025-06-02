#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=capstone_align_var_call
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=2gb
#SBATCH --output=/scratch/sjc78466/capstone/job_logs/log.%j

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sjc78466@uga.edu

#######################################################################################
###   LOADING ALL SOFTWARE MODULES NEEDED FOR READ ALIGNMENT AND VARIANT CALLING    ###
#######################################################################################

module load BWA/0.7.18-GCCcore-13.3.0 # Load BWA
module load SAMtools/1.18-GCC-12.3.0 # load SAMtools
module load BCFtools/1.18-GCC-12.3.0 # load BCFtools

##############################
###   DEFINING VARIABLES   ###
##############################

# Variable for absolute path again so don't need to keep typing it out
workdir="/scratch/sjc78466/capstone"

# Setting refseq E. coli genome URL to a variable so easier to work with
genome_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz"

############################################
###   SETTING UP NECESSARY DIRECTORIES   ###
############################################

mkdir $workdir/data/genomes
mkdir $workdir/results/sam
mkdir $workdir/results/bam
mkdir $workdir/results/vcf
mkdir $workdir/results/bcf

##############################################################################################

# Download refseq E. coli genome from NCBI
wget -O $workdir/data/genomes/ecoli_rel606.fna.gz $genome_url
gunzip $workdir/data/genomes/ecoli_rel606.fna.gz # decompressing the genome file

# Index the E. coli genome
bwa index $workdir/data/genomes/ecoli_rel606.fna #calling BWA and telling it to index the unzipped E. coli genome

# Loop over reads and do alignment and variant calling
for fwd in $workdir/data/trimmed_fastq/*_1.paired.fastq.gz
do
	sample=$(basename $fwd _1.paired.fastq.gz)

	# Run alignment
	echo "Aligning $sample" #echo here so we can keep tabs on progress of the program
	rev="/scratch/sjc78466/capstone/data/trimmed_fastq/${sample}_2.paired.fastq.gz" #setting variable for the reverse read
	bwa mem $workdir/data/genomes/ecoli_rel606.fna $fwd $rev > $workdir/results/sam/$sample.sam
	#break

	# Convert to BAM and sort
	samtools view -S -b $workdir/results/sam/$sample.sam > $workdir/results/bam/$sample.bam
	samtools sort -o $workdir/results/bam/$sample.sorted.bam $workdir/results/bam/$sample.bam

	# Do variant calling
	echo "Calling variants in $sample"
	bcftools mpileup -O b -o $workdir/results/bcf/$sample.bcf -f $workdir/data/genomes/ecoli_rel606.fna $workdir/results/bam/$sample.sorted.bam
	bcftools call --ploidy 1 -m -v -o $workdir/results/vcf/$sample.vcf $workdir/results/bcf/$sample.bcf

done
