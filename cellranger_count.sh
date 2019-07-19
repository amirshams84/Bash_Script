#!/bin/bash
# sbatch --module cellranger/3.0.2 -g 60 -t 32 --time 24:00:00 SCAF797_Ctrl-1.sh

DATA_DIR="/data/dadkhahe/single_cell_sample_data/fastq/Seq1_HK7V3BGXB/HK7V3BGXB"
TRANSCRIPTOME=/fdb/cellranger/refdata-cellranger-mm10-3.0.0
OUTPUT_FOLDER=/data/dadkhahe/single_cell_sample_data/count
SAMPLE_SHEET=/data/dadkhahe/single_cell_sample_data/metadata/input_samplesheet.csv
module load cellranger/3.0.2

for sample_name in $(cut -d"," -f2 $SAMPLE_SHEET | grep -v 'Sample')
do
	cellranger count \
	--id=$sample_name \
	--fastqs=$DATA_DIR/ \
	--transcriptome=$TRANSCRIPTOME \
	--localcores=32 \
	--localmem=60 \
	--sample=$sample_name

	mv ./$sample_name $OUTPUT_FOLDER
done



