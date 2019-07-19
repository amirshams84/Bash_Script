#!/bin/bash
# sbatch --module cellranger/3.0.2 -g 60 -t 32 --time 24:00:00 SCAF797_Ctrl-1.sh

DATA_DIR="/data/dadkhahe/single_cell_sample_data/fastq/HLYVNBGXB"
COUNT_FOLDER=/data/dadkhahe/single_cell_sample_data/count
TARGET=$DATA_DIR/$SAMPLE_SET
TRANSCRIPTOME=/fdb/cellranger/refdata-cellranger-mm10-3.0.0
OUTPUT_FOLDER=/data/dadkhahe/single_cell_sample_data/aggr
SAMPLE_SHEET=/data/dadkhahe/single_cell_sample_data/metadata/input_samplesheet.csv
#building aggr file
#each study has aggr.csv file including two columns 1 sample name 2. the path  to a file in $sample_name/count/outs 
#named molecule_info.h5
# we are going to build this file with the following script
# first we search our study folder for all generated molecule_info.h5 file with "find" command
#then we are going to extract the sample name from the path of molecule_info.h5 file
#this method is safe so we are not confusing mulecule file with different study
AGGR_INPUT_SHEET=/data/dadkhahe/single_cell_sample_data/aggr/aggr.csv
echo "library_id,molecule_h5" > $AGGR_INPUT_SHEET
for molecule_info_path in $(find $COUNT_FOLDER -name "molecule_info.h5")
do
	#echo $molecule_info_path
	SAMPLE_NAME=$(echo $molecule_info_path | sed 's/\/outs.*//')
	#this command remove everything after /outs including /outs
	SAMPLE_NAME=$(echo $SAMPLE_NAME | sed 's/.*SCAF/SCAF/')
	#this command replace everything before SCAF with just SCAF
	#echo $SAMPLE_NAME
	#printhing the sample name 

	#then we are going the concat sample_name with molecule_info.h5 path
	echo "$SAMPLE_NAME,$molecule_info_path" >> $AGGR_INPUT_SHEET

done

module load cellranger/3.0.2

cellranger aggr \
--id=AggregatedDatasets \
--csv=$AGGR_INPUT_SHEET \
--normalize=mapped

mv ./AggregatedDatasets /data/dadkhahe/single_cell_sample_data/aggr/
#/data/dadkhahe/single_cell_sample_data/count/SCAF802_CTLA-2/outs/molecule_info.h5


