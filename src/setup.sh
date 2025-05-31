#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=capstone_setup
#SBATCH --ntasks=1
#SBATCH --time=0:30:00
#SBATCH --mem=2gb
#SBATCH --output=/scratch/sjc78466/capstone/log.%j

#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sjc78466@uga.edu

# Copy raw E. coli short read sequence data from BINF 8960 instructor_data directory and make read-only

# Variable for workdir on cluster b/c want to have absolute path in script out of an abundance of caution
workdir="/scratch/sjc78466/capstone"

# Copying raw reads for three E. coli isolates & adapters used for sequencing run to capstone/work directory
cp -r /work/binf8960/instructor_data/raw_fastq $workdir/data

# Changing file permissions for raw data to be read-only so cannot be modified
chmod -w $workdir/data/raw_fastq/*.fastq*
