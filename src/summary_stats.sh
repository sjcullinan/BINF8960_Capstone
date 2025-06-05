#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=capstone_align_var_call
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=2gb
#SBATCH --output=/scratch/sjc78466/capstone/job_logs/log.%j

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sjc78466@uga.edu

# Loading required software
module load SAMtools/1.18-GCC-12.3.0

# Setting variable for absolute path of working directory
workdir="/scratch/sjc78466/capstone"

##########################################################################################################################################
###   COLLECTING SUMMARY STATISTICS DATA ABOUT RESULTS FROM TRIMMING, ALIGNING SAMPLE READS TO REFERENCE GENOME, AND VARIANT CALLING   ###
##########################################################################################################################################

# READS IN RAW DATA FOR EACH SAMPLE DETERMINED WITH MULTIQC HTML REPORT ON RAW READS BY REFERENCING THE SEQUENCE COUNTS SECTION, DIRECTLY BELOW THE FASTQC HEADER THAT FOLLOWS THE GENERAL STATS HEADER
# SRR2584863_1: 1,120,742 unique reads & 432,517 duplicate reads
# SRR2584863_2: 1,180,111 unique reads & 373,148 duplicate reads
# SRR2584866_1: 1,755,424 unique reads & 1,012,974 duplicate reads
# SRR2584866_2: 1,686,058 unique reads & 1,082,340 duplicate reads
# SRR2589044_1: 865,621 unique reads & 241,469 duplicate reads
# SRR2589044_2: 896,687 unique reads & 210,412 duplicate reads



# READS IN TRIMMED PAIRED READ DATA FOR EACH SAMPLE DETERMINED WITH MULTIQC HTML REPORT ON TRIMMED PAIRED READS BY REFERENCING THE SEQUENCE COUNTS SECTION, DIRECTLY BELOW THE FASTQC HEADER THAT FOLLOWS THE GENERAL STATS HEADER
# SRR2584863_1.paired: 1,124,156 unique reads & 361,804 duplicate reads
# SRR2584863_2.paired: 1,180,986 unique reads & 304,974 duplicate reads
# SRR2584866_1.paired: 1,887,684 unique reads & 821,988 duplicate reads
# SRR2584866_2.paired: 1,821,620 unique reads & 888,052 duplicate reads
# SRR2589044_1.paired: 865,494 unique reads & 201,225 duplicate reads
# SRR2589044_2.paired: 894,309 unique reads & 172,410 duplicate reads

# COUNTING NUMBER OF READS THAT ALIGNED TO REFERENCE GENOME FOR EACH SAMPLE
#	Indexing BAM files
for bamfile in $workdir/results/bam/*.sorted.bam
do
	samtools index $bamfile
done

#	samtools view to count reads that aligned to reference genome
samtools view -F 0x4 -c $workdir/results/bam/SRR2584863.sorted.bam CP000819.1 -o $workdir/results/mapped_reads_4863.txt
samtools view -F 0x4 -c $workdir/results/bam/SRR2584866.sorted.bam CP000819.1 -o $workdir/results/mapped_reads_4866.txt
samtools view -F 0x4 -c $workdir/results/bam/SRR2589044.sorted.bam CP000819.1 -o $workdir/results/mapped_reads_9044.txt

# COUNTING NUMBER OF VARIANT SITES IN EACH SAMPLE for sample in $workdir/results/vcf/*.vcf
for sample in $workdir/results/vcf/*.vcf
do
	grep -H -c CP000819.1 $sample >> results/var_sites.txt
done

# Using sed to kinda manually subtract 1 from all values to remove the extra count from the *.vcf file header to obtain true counts of the number of variants per sample
#
#
#	Setting variables b/c trying to find and replace numbers (command was not working with numerical values directly, but that may be because I missed something else on my part!)

var_count_4863=31
true_var_count_4863=30
var_count_4866=819
true_var_count_4866=818
var_count_9044=15
true_var_count_9044=14

#	Finding and replacing variant site counts with the accurate counts

sed -i -e "s/$var_count_4863/$true_var_count_4863/g" $workdir/results/var_sites.txt
sed -i -e "s/$var_count_4866/$true_var_count_4866/g" $workdir/results/var_sites.txt
set -i -e "s/$var_count_9044/$true_var_count_9044/g" $workdir/results/var_sites.txt

########################################################################
###   OUTPUTTING SUMMARY STATISTICS TO A TIDY-FORMATTED *.CSV FILE   ###
########################################################################

# Setting samples object to hold _______
#samples="_________________"

for read in $samples
do

done



