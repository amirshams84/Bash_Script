#!/bin/bash
#sbatch --mem=50G --cpus-per-task=32 --partition=norm --time=1-00:00:00 cellranger_mkfastq.sh


DATA_DIR=/data/dadkhahe/single_cell_sample_data
RUN_FOLDER=190615_NB552150_0005_AHLYVNBGXB
SAMPLE_SHEET=/data/dadkhahe/single_cell_sample_data/metadata/input_samplesheet.csv
OUTPUT_DIR=/data/dadkhahe/single_cell_sample_data/fastq
TARGET=$DATA_DIR/$RUN_FOLDER

module load cellranger/3.0.2


cellranger mkfastq \
--run=$TARGET \
--csv=$SAMPLE_SHEET \
--localcores=16 \
--localmem=15 \
--output-dir=$OUTPUT_DIR \
--qc
