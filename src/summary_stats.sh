#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=capstone_align_var_call
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=2gb
#SBATCH --output=/scratch/sjc78466/capstone/job_logs/log.%j

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sjc78466@uga.edu

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


# COUNTING NUMBER OF VARIANT SITES IN EACH SAMPLE


########################################################################
###   OUTPUTTING SUMMARY STATISTICS TO A TIDY-FORMATTED *.CSV FILE   ###
########################################################################




