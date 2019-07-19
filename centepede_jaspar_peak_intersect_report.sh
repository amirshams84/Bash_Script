#! /bin/bash

module load singularity
module load bedtools
PROCESSORS=$((SLURM_CPUS_PER_TASK - 2))
export SINGULARITY_BINDPATH="/gpfs,/spin1,/data/$USER:/data,/scratch,/fdb,/lscratch"

CENTEPEDE_IMAGE=/data/RTB/datashare/Amir/Projects/CENTEPEDE/centepede_v1.simg
JASPAR_DB=/data/$USER/CENTEPEDE/JASPAR
JASPAR_SINGULARITY_DB=/data/CENTEPEDE/JASPAR
WORKDIR_SINGULARITY=/data/CENTEPEDE
WORKDIR=/data/$USER/CENTEPEDE
mkdir -p $WORKDIR/MATRIX
mkdir -p $WORKDIR/Centepede_Result


CELL1_PEAK_FILE=/data/RTB/datashare/Amir/Projects/jfarber-20160119-P1/spsingh-20180828-D1/NAIIVE/hg38/peak/macs2/pooled_rep/NAIIVE.narrowPeak.gz
CELL2_PEAK_FILE=/data/RTB/datashare/Amir/Projects/jfarber-20160119-P1/spsingh-20180828-D1/VA2.7+CCR6+CCR2-/hg38/peak/macs2/pooled_rep/VA2.7+CCR6+CCR2-.narrowPeak.gz
CELL3_PEAK_FILE=/data/RTB/datashare/Amir/Projects/jfarber-20160119-P1/spsingh-20180828-D1/VA2.7+CCR6+CCR2+/hg38/peak/macs2/pooled_rep/VA2.7+CCR6+CCR2+.narrowPeak.gz
CELL4_PEAK_FILE=/data/RTB/datashare/Amir/Projects/jfarber-20160119-P1/spsingh-20180828-D1/NAIIVE/hg38/peak/macs2/pooled_rep/NAIIVE.narrowPeak.gz
CELL5_PEAK_FILE=/data/RTB/datashare/Amir/Projects/jfarber-20160119-P1/spsingh-20180828-D1/CCR6-CCR2-/hg38/peak/macs2/pooled_rep/CCR6-CCR2-.narrowPeak.gz
CELL6_PEAK_FILE=/data/RTB/datashare/Amir/Projects/jfarber-20160119-P1/spsingh-20180828-D1/VA2.7-CCR6+/hg38/peak/macs2/pooled_rep/VA2.7-CCR6+.narrowPeak.gz



CELL1_BAM_FILE_SINGULARITY=$WORKDIR_SINGULARITY/BAMS/NAIIVE.processed.bam
CELL2_BAM_FILE_SINGULARITY=$WORKDIR_SINGULARITY/BAMS/VA2.7+CCR6+CCR2-.processed.bam
CELL3_BAM_FILE_SINGULARITY=$WORKDIR_SINGULARITY/BAMS/VA2.7+CCR6+CCR2+.processed.bam
CELL4_BAM_FILE_SINGULARITY=$WORKDIR_SINGULARITY/BAMS/NAIIVE.processed.bam
CELL5_BAM_FILE_SINGULARITY=$WORKDIR_SINGULARITY/BAMS/CCR6-CCR2-.processed.bam
CELL6_BAM_FILE_SINGULARITY=$WORKDIR_SINGULARITY/BAMS/VA2.7-CCR6+.processed.bam

CELL1_TITLE=Jurkat
CELL2_TITLE=Mait_Minus
CELL3_TITLE=Mait_Plus
CELL4_TITLE=Naive
CELL5_TITLE=Tconv_Minus
CELL6_TITLE=Tconv_Plus


key_term='';

BOUND_POSTERIOR_PROBS_LIMIT=0.95

BOUND_SITES_REPORT=/data/$USER/CENTEPEDE/bound_sites_report.txt
NEWLINE=$'\n'
TAB=$'\t'
for each_motif in $JASPAR_DB/*.meme
do
	#####################################
	# Grab the Motif Info
	#####################################
	motif_title=$(grep 'MOTIF' $each_motif | cut -f2-3 -d ' ' | sed -e 's/ /_/')
	motif_info=$(grep 'MOTIF' $each_motif| cut -f2-3 -d ' ' | sed -e 's/[^[:alnum:]]/_/g')
	echo $motif_title
	lower_motif=$(echo $motif_info | sed -e 's/\(.*\)/\L\1/')
	if [[ $lower_motif = *"$key_term"* ]];
		then
		motif_length=$(grep 'letter-probability matrix' $each_motif | cut -f6 -d ' ')
		filename=$(basename $each_motif)
		TARGET_SCANNED_MOTIF_BED_SINGULARITY=$JASPAR_SINGULARITY_DB/${filename%.meme}_out/fimo.bed.gz
		TARGET_SCANNED_MOTIF_BED=$JASPAR_DB/${filename%.meme}_out/fimo.bed.gz

		bound_sites="$motif_title"
		bound_sites+=$TAB
		####################################
		# INTERSECT WITH ATAC-SEQ PEAKS
		####################################
		peak_number=$(zcat $TARGET_SCANNED_MOTIF_BED | wc -l)
		if [[ "$peak_number" < 1  ]];
			then
				bound_sites+="Empty records"
				echo "$bound_sites" >> $BOUND_SITES_REPORT
				continue
		fi
		bound_sites+=$(zcat $TARGET_SCANNED_MOTIF_BED | wc -l)
		bound_sites+=$TAB

		bound_sites+=$(zcat $CELL1_PEAK_FILE | wc -l)
		bound_sites+=$TAB

		bound_sites+=$(zcat $CELL2_PEAK_FILE | wc -l)
		bound_sites+=$TAB

		bound_sites+=$(zcat $CELL3_PEAK_FILE | wc -l)
		bound_sites+=$TAB

		bound_sites+=$(zcat $CELL4_PEAK_FILE | wc -l)
		bound_sites+=$TAB

		bound_sites+=$(zcat $CELL5_PEAK_FILE | wc -l)
		bound_sites+=$TAB

		bound_sites+=$(zcat $CELL6_PEAK_FILE | wc -l)
		bound_sites+=$TAB
		{

		echo "Intersecting with peaks"
		CELL1_SCANNED_MOTIF_BED_INTERSECT_PEAKS=$WORKDIR/Centepede_Result/${CELL1_TITLE}.${motif_info}.intersect_peak.bed.gz
		CELL1_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY=$WORKDIR_SINGULARITY/Centepede_Result/${CELL1_TITLE}.${motif_info}.intersect_peak.bed.gz
		intersectBed -a $TARGET_SCANNED_MOTIF_BED -b $CELL1_PEAK_FILE -wa | gzip > $CELL1_SCANNED_MOTIF_BED_INTERSECT_PEAKS
		
		bound_sites+=$(zcat $CELL1_SCANNED_MOTIF_BED_INTERSECT_PEAKS | wc -l)
		bound_sites+=$TAB

		CELL2_SCANNED_MOTIF_BED_INTERSECT_PEAKS=$WORKDIR/Centepede_Result/${CELL2_TITLE}.${motif_info}.intersect_peak.bed.gz
		CELL2_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY=$WORKDIR_SINGULARITY/Centepede_Result/${CELL2_TITLE}.${motif_info}.intersect_peak.bed.gz
		intersectBed -a $TARGET_SCANNED_MOTIF_BED -b $CELL2_PEAK_FILE -wa | gzip > $CELL2_SCANNED_MOTIF_BED_INTERSECT_PEAKS
		
		bound_sites+=$(zcat $CELL2_SCANNED_MOTIF_BED_INTERSECT_PEAKS | wc -l)
		bound_sites+=$TAB

		CELL3_SCANNED_MOTIF_BED_INTERSECT_PEAKS=$WORKDIR/Centepede_Result/${CELL3_TITLE}.${motif_info}.intersect_peak.bed.gz
		CELL3_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY=$WORKDIR_SINGULARITY/Centepede_Result/${CELL3_TITLE}.${motif_info}.intersect_peak.bed.gz
		intersectBed -a $TARGET_SCANNED_MOTIF_BED -b $CELL3_PEAK_FILE -wa | gzip > $CELL3_SCANNED_MOTIF_BED_INTERSECT_PEAKS
		
		bound_sites+=$(zcat $CELL3_SCANNED_MOTIF_BED_INTERSECT_PEAKS | wc -l)
		bound_sites+=$TAB

		CELL4_SCANNED_MOTIF_BED_INTERSECT_PEAKS=$WORKDIR/Centepede_Result/${CELL4_TITLE}.${motif_info}.intersect_peak.bed.gz
		CELL4_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY=$WORKDIR_SINGULARITY/Centepede_Result/${CELL4_TITLE}.${motif_info}.intersect_peak.bed.gz
		intersectBed -a $TARGET_SCANNED_MOTIF_BED -b $CELL4_PEAK_FILE -wa | gzip > $CELL4_SCANNED_MOTIF_BED_INTERSECT_PEAKS
		
		bound_sites+=$(zcat $CELL4_SCANNED_MOTIF_BED_INTERSECT_PEAKS | wc -l)
		bound_sites+=$TAB

		CELL5_SCANNED_MOTIF_BED_INTERSECT_PEAKS=$WORKDIR/Centepede_Result/${CELL5_TITLE}.${motif_info}.intersect_peak.bed.gz
		CELL5_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY=$WORKDIR_SINGULARITY/Centepede_Result/${CELL5_TITLE}.${motif_info}.intersect_peak.bed.gz
		intersectBed -a $TARGET_SCANNED_MOTIF_BED -b $CELL5_PEAK_FILE -wa | gzip > $CELL5_SCANNED_MOTIF_BED_INTERSECT_PEAKS
		
		bound_sites+=$(zcat $CELL5_SCANNED_MOTIF_BED_INTERSECT_PEAKS | wc -l)
		bound_sites+=$TAB

		CELL6_SCANNED_MOTIF_BED_INTERSECT_PEAKS=$WORKDIR/Centepede_Result/${CELL6_TITLE}.${motif_info}.intersect_peak.bed.gz
		CELL6_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY=$WORKDIR_SINGULARITY/Centepede_Result/${CELL6_TITLE}.${motif_info}.intersect_peak.bed.gz
		intersectBed -a $TARGET_SCANNED_MOTIF_BED -b $CELL6_PEAK_FILE -wa | gzip > $CELL6_SCANNED_MOTIF_BED_INTERSECT_PEAKS
		
		bound_sites+=$(zcat $CELL6_SCANNED_MOTIF_BED_INTERSECT_PEAKS | wc -l)
		bound_sites+=$TAB

		echo "Done!!!"
		} || {
			bound_sites+="Unable to intersect"
			echo "$bound_sites" >> $BOUND_SITES_REPORT
			continue
		}
		sleep 1
		####################################
		# Making Discreet Cut Matrix
		####################################
		{
		echo "Making discreet cut matrix"
		CELL1_CUT_MATRIX=$WORKDIR/MATRIX/${CELL1_TITLE}.${motif_info}.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL1_BAM_FILE_SINGULARITY $CELL1_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY | gzip > $CELL1_CUT_MATRIX

		CELL2_CUT_MATRIX=$WORKDIR/MATRIX/${CELL2_TITLE}.${motif_info}.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL2_BAM_FILE_SINGULARITY $CELL2_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY | gzip > $CELL2_CUT_MATRIX

		CELL3_CUT_MATRIX=$WORKDIR/MATRIX/${CELL3_TITLE}.${motif_info}.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL3_BAM_FILE_SINGULARITY $CELL3_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY | gzip > $CELL3_CUT_MATRIX

		CELL4_CUT_MATRIX=$WORKDIR/MATRIX/${CELL4_TITLE}.${motif_info}.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL4_BAM_FILE_SINGULARITY $CELL4_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY | gzip > $CELL4_CUT_MATRIX

		CELL5_CUT_MATRIX=$WORKDIR/MATRIX/${CELL5_TITLE}.${motif_info}.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL5_BAM_FILE_SINGULARITY $CELL5_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY | gzip > $CELL5_CUT_MATRIX

		CELL6_CUT_MATRIX=$WORKDIR/MATRIX/${CELL6_TITLE}.${motif_info}.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL6_BAM_FILE_SINGULARITY $CELL6_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY | gzip > $CELL6_CUT_MATRIX
		
		echo "Done!!!"
		} || {
			bound_sites+="Unable to Making discreet cut matrix"
			echo "$bound_sites" >> $BOUND_SITES_REPORT
			continue
		}
		sleep 1
		####################################
		# RUN CENTEPEDE
		####################################
		{
		echo "Running centipede"
		CELL1_CUT_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL1_TITLE}.${motif_info}.matrix.gz
		CELL1_CENTEPEDE_DISCREET_TABLE_SINGULARITY=$WORKDIR_SINGULARITY/Centepede_Result/${CELL1_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL1_CENTEPEDE_DISCREET_TABLE=$WORKDIR/Centepede_Result/${CELL1_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL1_CENTEPEDE_PDF=$WORKDIR_SINGULARITY/Centepede_Result/${CELL1_TITLE}.${motif_info}.motif_centepede.pdf
		singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL1_CUT_MATRIX $CELL1_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY $CELL1_CENTEPEDE_DISCREET_TABLE_SINGULARITY $CELL1_CENTEPEDE_PDF $motif_length

		bound_sites+=$(zcat $CELL1_CENTEPEDE_DISCREET_TABLE| wc -l)
		bound_sites+=$TAB

		CELL2_CUT_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL2_TITLE}.${motif_info}.matrix.gz
		CELL2_CENTEPEDE_DISCREET_TABLE_SINGULARITY=$WORKDIR_SINGULARITY/Centepede_Result/${CELL2_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL2_CENTEPEDE_DISCREET_TABLE=$WORKDIR/Centepede_Result/${CELL2_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL2_CENTEPEDE_PDF=$WORKDIR_SINGULARITY/Centepede_Result/${CELL2_TITLE}.${motif_info}.motif_centepede.pdf
		singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL2_CUT_MATRIX $CELL2_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY $CELL2_CENTEPEDE_DISCREET_TABLE_SINGULARITY $CELL2_CENTEPEDE_PDF $motif_length

		bound_sites+=$(zcat $CELL2_CENTEPEDE_DISCREET_TABLE | wc -l)
		bound_sites+=$TAB
		
		CELL3_CUT_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL3_TITLE}.${motif_info}.matrix.gz
		CELL3_CENTEPEDE_DISCREET_TABLE_SINGULARITY=$WORKDIR_SINGULARITY/Centepede_Result/${CELL3_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL3_CENTEPEDE_DISCREET_TABLE=$WORKDIR/Centepede_Result/${CELL3_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL3_CENTEPEDE_PDF=$WORKDIR_SINGULARITY/Centepede_Result/${CELL3_TITLE}.${motif_info}.motif_centepede.pdf
		singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL3_CUT_MATRIX $CELL3_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY $CELL3_CENTEPEDE_DISCREET_TABLE_SINGULARITY $CELL3_CENTEPEDE_PDF $motif_length

		bound_sites+=$(zcat $CELL3_CENTEPEDE_DISCREET_TABLE | wc -l)
		bound_sites+=$TAB
		
		CELL4_CUT_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL4_TITLE}.${motif_info}.matrix.gz
		CELL4_CENTEPEDE_DISCREET_TABLE_SINGULARITY=$WORKDIR_SINGULARITY/Centepede_Result/${CELL4_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL4_CENTEPEDE_DISCREET_TABLE=$WORKDIR/Centepede_Result/${CELL4_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL4_CENTEPEDE_PDF=$WORKDIR_SINGULARITY/Centepede_Result/${CELL4_TITLE}.${motif_info}.motif_centepede.pdf
		singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL4_CUT_MATRIX $CELL4_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY $CELL4_CENTEPEDE_DISCREET_TABLE_SINGULARITY $CELL4_CENTEPEDE_PDF $motif_length
		
		bound_sites+=$(zcat $CELL4_CENTEPEDE_DISCREET_TABLE | wc -l)
		bound_sites+=$TAB

		CELL5_CUT_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL5_TITLE}.${motif_info}.matrix.gz
		CELL5_CENTEPEDE_DISCREET_TABLE_SINGULARITY=$WORKDIR_SINGULARITY/Centepede_Result/${CELL5_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL5_CENTEPEDE_DISCREET_TABLE=$WORKDIR/Centepede_Result/${CELL5_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL5_CENTEPEDE_PDF=$WORKDIR_SINGULARITY/Centepede_Result/${CELL5_TITLE}.${motif_info}.motif_centepede.pdf
		singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL5_CUT_MATRIX $CELL5_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY $CELL5_CENTEPEDE_DISCREET_TABLE_SINGULARITY $CELL5_CENTEPEDE_PDF $motif_length

		bound_sites+=$(zcat $CELL5_CENTEPEDE_DISCREET_TABLE | wc -l)
		bound_sites+=$TAB
		
		CELL6_CUT_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL6_TITLE}.${motif_info}.matrix.gz
		CELL6_CENTEPEDE_DISCREET_TABLE_SINGULARITY=$WORKDIR_SINGULARITY/Centepede_Result/${CELL6_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL6_CENTEPEDE_DISCREET_TABLE=$WORKDIR/Centepede_Result/${CELL6_TITLE}.${motif_info}.motif_centepede.bed.gz
		CELL6_CENTEPEDE_PDF=$WORKDIR_SINGULARITY/Centepede_Result/${CELL6_TITLE}.${motif_info}.motif_centepede.pdf
		singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL6_CUT_MATRIX $CELL6_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY $CELL6_CENTEPEDE_DISCREET_TABLE_SINGULARITY $CELL6_CENTEPEDE_PDF $motif_length
		
		bound_sites+=$(zcat $CELL6_CENTEPEDE_DISCREET_TABLE | wc -l)
		bound_sites+=$TAB
		
		echo "Done!!!"
		} || {
			bound_sites+="Unable to Running centipede"
			echo "$bound_sites" >> $BOUND_SITES_REPORT
			continue
		}
		sleep 1
		####################################
		# FILTER CENTEPEDE
		####################################
		
		echo "Filtering Centipede result based on Posterior probability"
		
		CELL1_CENTEPEDE_BOUND_TABLE=$WORKDIR/Centepede_Result/${CELL1_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		
		zcat $CELL1_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL1_CENTEPEDE_BOUND_TABLE
		bound_sites+=$(zcat $CELL1_CENTEPEDE_BOUND_TABLE | wc -l)
		bound_sites+=$TAB

		
		CELL2_CENTEPEDE_BOUND_TABLE=$WORKDIR/Centepede_Result/${CELL2_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		
		zcat $CELL2_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL2_CENTEPEDE_BOUND_TABLE
		bound_sites+=$(zcat $CELL2_CENTEPEDE_BOUND_TABLE | wc -l)
		bound_sites+=$TAB

		
		CELL3_CENTEPEDE_BOUND_TABLE=$WORKDIR/Centepede_Result/${CELL3_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		
		zcat $CELL3_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL3_CENTEPEDE_BOUND_TABLE
		bound_sites+=$(zcat $CELL3_CENTEPEDE_BOUND_TABLE | wc -l)
		bound_sites+=$TAB

		
		CELL4_CENTEPEDE_BOUND_TABLE=$WORKDIR/Centepede_Result/${CELL4_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		
		zcat $CELL4_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL4_CENTEPEDE_BOUND_TABLE
		bound_sites+=$(zcat $CELL4_CENTEPEDE_BOUND_TABLE | wc -l)
		bound_sites+=$TAB

		CELL5_CENTEPEDE_BOUND_TABLE=$WORKDIR/Centepede_Result/${CELL5_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		
		zcat $CELL5_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL5_CENTEPEDE_BOUND_TABLE
		bound_sites+=$(zcat $CELL5_CENTEPEDE_BOUND_TABLE | wc -l)
		bound_sites+=$TAB

		
		CELL6_CENTEPEDE_BOUND_TABLE=$WORKDIR/Centepede_Result/${CELL6_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		
		zcat $CELL6_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL6_CENTEPEDE_BOUND_TABLE
		bound_sites+=$(zcat $CELL6_CENTEPEDE_BOUND_TABLE | wc -l)
		
		

		echo "$bound_sites" >> $BOUND_SITES_REPORT
		
		echo "Done!!!"
		sleep 1
		####################################
		# Making Aggregate Cut Matrix 
		####################################
		echo "Making Aggrgate cut matrix"
		CELL1_CENTEPEDE_BOUND_TABLE=$WORKDIR_SINGULARITY/Centepede_Result/${CELL1_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		CELL1_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR/MATRIX/${CELL1_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL1_BAM_FILE_SINGULARITY $CELL1_CENTEPEDE_BOUND_TABLE | gzip > $CELL1_CENTEPEDE_BOUND_AGGREGATE_MATRIX


		CELL2_CENTEPEDE_BOUND_TABLE=$WORKDIR_SINGULARITY/Centepede_Result/${CELL2_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		CELL2_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR/MATRIX/${CELL2_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL2_BAM_FILE_SINGULARITY $CELL2_CENTEPEDE_BOUND_TABLE | gzip > $CELL2_CENTEPEDE_BOUND_AGGREGATE_MATRIX


		CELL3_CENTEPEDE_BOUND_TABLE=$WORKDIR_SINGULARITY/Centepede_Result/${CELL3_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		CELL3_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR/MATRIX/${CELL3_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL3_BAM_FILE_SINGULARITY $CELL3_CENTEPEDE_BOUND_TABLE | gzip > $CELL3_CENTEPEDE_BOUND_AGGREGATE_MATRIX


		CELL4_CENTEPEDE_BOUND_TABLE=$WORKDIR_SINGULARITY/Centepede_Result/${CELL4_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		CELL4_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR/MATRIX/${CELL4_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL4_BAM_FILE_SINGULARITY $CELL4_CENTEPEDE_BOUND_TABLE | gzip > $CELL4_CENTEPEDE_BOUND_AGGREGATE_MATRIX

		CELL5_CENTEPEDE_BOUND_TABLE=$WORKDIR_SINGULARITY/Centepede_Result/${CELL5_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		CELL5_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR/MATRIX/${CELL5_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL5_BAM_FILE_SINGULARITY $CELL5_CENTEPEDE_BOUND_TABLE | gzip > $CELL5_CENTEPEDE_BOUND_AGGREGATE_MATRIX


		CELL6_CENTEPEDE_BOUND_TABLE=$WORKDIR_SINGULARITY/Centepede_Result/${CELL6_TITLE}.${motif_info}.motif_centepede.bound.bed.gz
		CELL6_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR/MATRIX/${CELL6_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL6_BAM_FILE_SINGULARITY $CELL6_CENTEPEDE_BOUND_TABLE | gzip > $CELL6_CENTEPEDE_BOUND_AGGREGATE_MATRIX
		echo "Done!!!"
		sleep 1
		####################################
		# PLOT BOUNDED 
		####################################
		echo "Plotting Bound sites"
		CELL1_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL1_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		CELL2_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL2_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		CELL3_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL3_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		CELL4_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL4_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		CELL5_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL5_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz
		CELL6_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$WORKDIR_SINGULARITY/MATRIX/${CELL6_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz



		CENTEPEDE_MOTIF_AGGREGATE_PDF=$WORKDIR_SINGULARITY/Centepede_Result/${CELL1_TITLE}.${CELL2_TITLE}.${CELL3_TITLE}.${CELL4_TITLE}.${CELL5_TITLE}.${CELL6_TITLE}.${motif_info}.motif_aggregate_centepede.pdf
		singularity exec $CENTEPEDE_IMAGE $WORKDIR_SINGULARITY/plot_aggregate_matrix_new.R $CELL1_CENTEPEDE_BOUND_AGGREGATE_MATRIX $CELL2_CENTEPEDE_BOUND_AGGREGATE_MATRIX $CELL3_CENTEPEDE_BOUND_AGGREGATE_MATRIX $CELL4_CENTEPEDE_BOUND_AGGREGATE_MATRIX $CELL5_CENTEPEDE_BOUND_AGGREGATE_MATRIX $CELL6_CENTEPEDE_BOUND_AGGREGATE_MATRIX "${motif_title} regions with motifs oriented by strand" $CENTEPEDE_MOTIF_AGGREGATE_PDF
		echo "Done!!!"
		
	fi

done

echo "Pipeline Successfully finished at: "$(date)


