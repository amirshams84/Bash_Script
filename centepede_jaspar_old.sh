#! /bin/bash
set -o pipefail
set -e

PROCESSORS=$((SLURM_CPUS_PER_TASK - 2))
export SINGULARITY_BINDPATH="/gpfs,/spin1,/data/$USER:/data,/scratch,/fdb,/lscratch"

CENTEPEDE_IMAGE=/data/RTB/datashare/Amir/Projects/CENTEPEDE/centepede_v1.simg
JASPAR_DB=/data/$USER/CENTEPEDE/JASPAR
JASPAR_SINGULARITY_DB=/data/CENTEPEDE/JASPAR
WORKDIR_SINGULARITY=/data/CENTEPEDE
WORKDIR=/data/$USER/CENTEPEDE
mkdir -p $WORKDIR/MATRIX
mkdir -p $WORKDIR/Centepede_Result
CELL1_BAM_FILE_SINGULARITY=/data/jfarber-20160119-P1/spsingh-20171019-D2/VA2.7+CCR6+CCR2+/hg38/align/rep1/AB2955-L1_R1.trim.PE2SE.nodup.bam
CELL2_BAM_FILE_SINGULARITY=/data/jfarber-20160119-P1/spsingh-20171019-D2/VA2.7+CCR6+CCR2-/hg38/align/rep1/AB2954-L1_R1.trim.PE2SE.nodup.bam
CELL1_TITLE=MAIT_PLUS
CELL2_TITLE=MAIT_MINUS
key_term='cebp';

for each_motif in $JASPAR_DB/*.meme
do
	#####################################
	# Grab the Motif Info
	#####################################
	motif_info=$(grep 'MOTIF' $each_motif | cut -f2-3 -d ' ' | sed -e 's/[^[:alnum:]]/_/g')
	lower_motif=$(echo $motif_info | sed -e 's/\(.*\)/\L\1/')
	echo $lower_motif
	if [[ $lower_motif = *"$key_term"* ]];
		then
		motif_length=$(grep 'letter-probability matrix' $each_motif | cut -f6 -d ' ')
		filename=$(basename $each_motif)
		TARGET_SCANNED_MOTIF_BED=$JASPAR_SINGULARITY_DB/${filename%.meme}_out/fimo.bed.gz
		####################################
		# Making Discreet Cut Matrix
		####################################
		CELL1_CUT_MATRIX=$WORKDIR/MATRIX/${CELL1_TITLE}.${motif_info}.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL1_BAM_FILE_SINGULARITY $TARGET_SCANNED_MOTIF_BED | gzip > $CELL1_CUT_MATRIX

		CELL2_CUT_MATRIX=$WORKDIR/MATRIX/${CELL2_TITLE}.${motif_info}.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL2_BAM_FILE_SINGULARITY $TARGET_SCANNED_MOTIF_BED | gzip > $CELL2_CUT_MATRIX
		####################################
		# RUN CENTEPEDE
		####################################
		CELL1_CUT_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL1_TITLE}.${motif_info}.matrix.gz
		CELL1_CENTEPEDE_DISCREET_TABLE=$WORKDIR_SINGULARITY/Centepede_Result/${CELL1_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL1_CENTEPEDE_PDF=$WORKDIR_SINGULARITY/Centepede_Result/${CELL1_TITLE}.${motif_info}.motif_centepede.pdf
		singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL1_CUT_MATRIX $TARGET_SCANNED_MOTIF_BED $CELL1_CENTEPEDE_DISCREET_TABLE $CELL1_CENTEPEDE_PDF $motif_length

		CELL2_CUT_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL2_TITLE}.${motif_info}.matrix.gz
		CELL2_CENTEPEDE_DISCREET_TABLE=$WORKDIR_SINGULARITY/Centepede_Result/${CELL2_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL2_CENTEPEDE_PDF=$WORKDIR_SINGULARITY/Centepede_Result/${CELL2_TITLE}.${motif_info}.motif_centepede.pdf
		singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL2_CUT_MATRIX $TARGET_SCANNED_MOTIF_BED $CELL2_CENTEPEDE_DISCREET_TABLE $CELL2_CENTEPEDE_PDF $motif_length
		####################################
		# FILTER CENTEPEDE
		####################################
		CELL1_CENTEPEDE_DISCREET_TABLE=$WORKDIR/Centepede_Result/${CELL1_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL1_CENTEPEDE_BOUND_TABLE=$WORKDIR/Centepede_Result/${CELL1_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		BOUND_POSTERIOR_PROBS_LIMIT=0.95
		zcat $CELL1_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL1_CENTEPEDE_BOUND_TABLE
		

		CELL2_CENTEPEDE_DISCREET_TABLE=$WORKDIR/Centepede_Result/${CELL2_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL2_CENTEPEDE_BOUND_TABLE=$WORKDIR/Centepede_Result/${CELL2_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		BOUND_POSTERIOR_PROBS_LIMIT=0.95
		zcat $CELL2_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL2_CENTEPEDE_BOUND_TABLE
		####################################
		# Making Aggregate Cut Matrix 
		####################################
		CELL1_CENTEPEDE_BOUND_TABLE=$WORKDIR_SINGULARITY/Centepede_Result/${CELL1_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		CELL2_CENTEPEDE_BOUND_TABLE=$WORKDIR_SINGULARITY/Centepede_Result/${CELL2_TITLE}.${motif_info}.motif_centepede.bound.bed.gz


		CELL1_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR/MATRIX/${CELL1_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		CELL2_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR/MATRIX/${CELL2_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz

		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL1_BAM_FILE_SINGULARITY $CELL1_CENTEPEDE_BOUND_TABLE | gzip > $CELL1_CENTEPEDE_BOUND_AGGREGATE_MATRIX

		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL2_BAM_FILE_SINGULARITY $CELL2_CENTEPEDE_BOUND_TABLE | gzip > $CELL2_CENTEPEDE_BOUND_AGGREGATE_MATRIX
		####################################
		# PLOT BOUNDED 
		####################################
		CELL1_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL1_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		CELL2_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL2_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		CENTEPEDE_MOTIF_AGGREGATE_PDF=$WORKDIR_SINGULARITY/Centepede_Result/${CELL1_TITLE}.VS.${CELL2_TITLE}.${motif_info}.motif_aggregate_centepede.pdf
		singularity exec $CENTEPEDE_IMAGE plot_aggregate_matrix.R $CELL1_CENTEPEDE_BOUND_AGGREGATE_MATRIX $CELL2_CENTEPEDE_BOUND_AGGREGATE_MATRIX "${CELL1_TITLE}.VS.${CELL2_TITLE} ${motif_info} regions with motifs oriented by strand" $CENTEPEDE_MOTIF_AGGREGATE_PDF

		break
	fi

done


