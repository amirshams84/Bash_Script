#!/bin/bash
#
#SBATCH --job-name=S01star
#SBATCH --mail-type=ALL

# sbatch --cpus-per-task=16 --mem=40g --time=24:00:00 RSEM_S01_STAR.txt
# sinteractive --cpus-per-task=16 --mem=40g
DATA_DIR="/data/SCAF/00-BLI_William/00_RawSequencingData/190710_M01303_0151_000000000-CGN4F/Data/Intensities/BaseCalls"

#this is amir solution
for fastq in ${DATA_DIR}/*R1*.fastq*
do
	fastq_name=$(basename $fastq)
	echo $fastq_name
	sample_name=${fastq_name%R1*}
	echo $sample_name
	forward_fastq=$fastq_name
	echo $forward_fastq
	reverse_fastq=${fastq_name/R1/R2}
	echo $reverse_fastq
done