#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=capstone_summary_stats
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=2gb
#SBATCH --output=/scratch/sjc78466/capstone/job_logs/log.%j

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sjc78466@uga.edu

##################
###   SET-UP   ###
##################

# Loading required software
module load SAMtools/1.18-GCC-12.3.0

# Setting variable for absolute path of working directory
workdir="/scratch/sjc78466/capstone"




##########################################################################################################################################
###   COLLECTING SUMMARY STATISTICS DATA ABOUT RESULTS FROM TRIMMING, ALIGNING SAMPLE READS TO REFERENCE GENOME, AND VARIANT CALLING   ###
##########################################################################################################################################

# READS IN RAW DATA FOR EACH SAMPLE DETERMINED WITH MULTIQC_FASTQ.TXT REPORT

cut -f 5 $workdir/results/raw_multiqc/multiqc_data/multiqc_fastqc.txt >> $workdir/results/raw_read_counts.txt
sed -n -i -e '2p;4p;6p' $workdir/results/raw_read_counts.txt




# READS IN TRIMMED PAIRED READ DATA FOR EACH SAMPLE DETERMINED WITH MULTIQC_FASTQC.TXT REPORT

cut -f 5 $workdir/results/trimmed_multiqc/multiqc_data/multiqc_fastqc.txt >> $workdir/results/trimmed_read_counts.txt
sed -n -i -e '2p;4p;6p' $workdir/results/trimmed_read_counts.txt




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

# Compile mapped read count data into single *.txt file in established sample order (SRR2584863 then SRR2584866 then SRR2589044)
cat $workdir/results/mapped_reads_4863.txt $workdir/results/mapped_reads_4866.txt $workdir/results/mapped_reads_9044.txt > $workdir/results/all_mapped_reads.txt

# COUNTING NUMBER OF VARIANT SITES IN EACH SAMPLE for sample in $workdir/results/vcf/*.vcf
for sample in $workdir/results/vcf/*.vcf
do
	grep -c CP000819.1 $sample >> $workdir/results/var_sites.txt #if had added -H flag here, could have added the file name
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
sed -i -e "s/$var_count_9044/$true_var_count_9044/g" $workdir/results/var_sites.txt





########################################################################
###   OUTPUTTING SUMMARY STATISTICS TO A TIDY-FORMATTED *.CSV FILE   ###
########################################################################

for sample in SRR2584863 SRR2584866 SRR2589044
do
	echo $sample >> $workdir/results/samples.txt
done

paste -d "," $workdir/results/samples.txt $workdir/results/raw_read_counts.txt $workdir/results/trimmed_read_counts.txt $workdir/results/all_mapped_reads.txt $workdir/results/var_sites.txt > $workdir/results/summary_stats.csv

# Manually adding headers in: Sample, Raw_Read_Counts, Trimmed_Read_Counts, Mapped_Read_Counts, Variant_Site_Counts
echo "Sample, Raw_Read_Counts, Trimmed_Read_Counts, Mapped_Read_Counts, Variant_Site_Counts" > $workdir/results/headers.csv
cat $workdir/results/summary_stats.csv >> $workdir/results/headers.csv
mv $workdir/results/headers.csv $workdir/results/summary_stats.csv
