#! /bin/bash

module load singularity
module load bedtools
PROCESSORS=$((SLURM_CPUS_PER_TASK - 2))
export SINGULARITY_BINDPATH="/gpfs,/spin1,/data/$USER:/data,/scratch,/fdb,/lscratch"

CENTEPEDE_IMAGE=/data/RTB/datashare/Amir/Projects/CENTEPEDE/centepede_v1.simg

WORKDIR_SINGULARITY=/data/CENTEPEDE
WORKDIR=/data/$USER/CENTEPEDE

SWARM_FILE=$WORKDIR/centipede_swarm.txt

MOTIF_DB=/data/$USER/CENTEPEDE/JASPAR
MOTIF_SINGULARITY_DB=/data/CENTEPEDE/JASPAR
MOTIF_FORMAT=meme
MOTIF_DB_TITLE=JASPAR

mkdir -p $WORKDIR
mkdir -p $WORKDIR/MOTIF_${MOTIF_DB_TITLE}_DATA
mkdir -p $WORKDIR/CENTIPEDE_DATA
mkdir -p $WORKDIR/CENTIPEDE_RESULTS


CELL1_PEAK_FILE=$WORKDIR/PEAKS/Jurkat.optimal.narrowPeak.gz
CELL2_PEAK_FILE=$WORKDIR/PEAKS/Mait_Minus.optimal.narrowPeak.gz
CELL3_PEAK_FILE=$WORKDIR/PEAKS/Mait_Plus.optimal.narrowPeak.gz
CELL4_PEAK_FILE=$WORKDIR/PEAKS/Naive.optimal.narrowPeak.gz
CELL5_PEAK_FILE=$WORKDIR/PEAKS/Tconv_Minus.optimal.narrowPeak.gz
CELL6_PEAK_FILE=$WORKDIR/PEAKS/Tconv_Plus.optimal.narrowPeak.gz



CELL1_BAM_FILE_SINGULARITY=$WORKDIR_SINGULARITY/BAMS/Jurkat.bam
CELL2_BAM_FILE_SINGULARITY=$WORKDIR_SINGULARITY/BAMS/Mait_Minus.bam
CELL3_BAM_FILE_SINGULARITY=$WORKDIR_SINGULARITY/BAMS/Mait_Plus.bam
CELL4_BAM_FILE_SINGULARITY=$WORKDIR_SINGULARITY/BAMS/Naive.bam
CELL5_BAM_FILE_SINGULARITY=$WORKDIR_SINGULARITY/BAMS/Tconv_Minus.bam
CELL6_BAM_FILE_SINGULARITY=$WORKDIR_SINGULARITY/BAMS/Tconv_Plus.bam

CELL1_TITLE=Jurkat
CELL2_TITLE=Mait_Minus
CELL3_TITLE=Mait_Plus
CELL4_TITLE=Naive
CELL5_TITLE=Tconv_Minus
CELL6_TITLE=Tconv_Plus

BOUND_SITES_REPORT=$WORKDIR/bound_sites_report.txt
BOUND_POSTERIOR_PROBS_LIMIT=0.95
######################################################
# CREATE BACKGROUND GENOME NUCLEOTIDES FREQUENCIES
######################################################
GENOME=/fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa
GENOME_TITLE=hg38
GENOME_SINGULARITY_BG=$WORKDIR_SINGULARITY/${GENOME_TITLE}.bg
GENOME_BG=$WORKDIR/${GENOME_TITLE}.bg
MARKOV_REPORT=$WORKDIR/${GENOME_TITLE}_fasta_get_markov_report.txt
if [ ! -f $GENOME_BG ]
	then
	singularity exec $CENTEPEDE_IMAGE fasta-get-markov $GENOME $GENOME_SINGULARITY_BG > $MARKOV_REPORT
fi


MEME_DB=$WORKDIR/MOTIF_${MOTIF_DB_TITLE}_DATA
MEME_SINGULARITY_DB=$WORKDIR_SINGULARITY/MOTIF_${MOTIF_DB_TITLE}_DATA

CENTIPEDE_DATA=$WORKDIR/CENTIPEDE_DATA
CENTIPEDE_SINGULARITY_DATA=$WORKDIR_SINGULARITY/CENTIPEDE_DATA
AWK_COMMAND=awk 'BEGIN {IFS="\t"; OFS="\t";} {print $1,$2,$3,$4,$5,$6}'
TAB=$'\t'
for each_motif in $MOTIF_DB/*.$MOTIF_FORMAT
do
	######################################################
	# pwm to meme motif conversion
	######################################################
	if [ $MOTIF_FORMAT == "pwm" ];
		then
			pwm=$(basename $each_motif)
			singularity exec $CENTEPEDE_IMAGE matrix2meme < <(tail -n+2 $each_motif | cut -f2-) > $MEME_DB/${pwm%.txt}.meme
			MEME=$MEME_DB/${pwm%.txt}.meme
			MEME_SINGULARITY=$MEME_SINGULARITY_DB/${pwm%.txt}.meme
		else
			cp $each_motif $MEME_DB/
			MEME=$MEME_DB/$(basename $each_motif)
			MEME_SINGULARITY=$MEME_SINGULARITY_DB/$(basename $each_motif)
	fi
	######################################################
	# meme motif process
	######################################################
	echo "module load singularity samtools bedtools bedops \
	singularity exec $CENTEPEDE_IMAGE fimo --verbosity 4 --oc ${MEME_SINGULARITY}_out --thresh 1e-4 --bgfile $GENOME_SINGULARITY_BG $MEME_SINGULARITY $GENOME \
	singularity exec $CENTEPEDE_IMAGE gff2bed < ${MEME_SINGULARITY}_out/fimo.gff | $AWK_COMMAND | gzip -c > ${MEME}_out/fimo.bed.gz \
	motif_title=$(grep 'MOTIF' $MEME | cut -f2-3 -d ' ' | sed -e 's/ /_/') \
	echo $motif_title \
	motif_info=$(grep 'MOTIF' $MEME| cut -f2-3 -d ' ' | sed -e 's/[^[:alnum:]]/_/g') \
	motif_length=$(grep 'letter-probability matrix' $MEME | cut -f6 -d ' ') \
	bound_sites='$motif_title' \
	bound_sites+=$TAB \
	$TARGET_SCANNED_MOTIF_BED=${MEME}_out/fimo.bed.gz \
	peak_number=$(zcat $TARGET_SCANNED_MOTIF_BED | wc -l) \
	bound_sites+=$TAB \
	bound_sites+=$(zcat $CELL1_PEAK_FILE | wc -l) \
	bound_sites+=$TAB \
	bound_sites+=$(zcat $CELL2_PEAK_FILE | wc -l) \
	bound_sites+=$TAB \
	bound_sites+=$(zcat $CELL3_PEAK_FILE | wc -l) \
	bound_sites+=$TAB \
	bound_sites+=$(zcat $CELL4_PEAK_FILE | wc -l) \
	bound_sites+=$TAB \
	bound_sites+=$(zcat $CELL5_PEAK_FILE | wc -l) \
	bound_sites+=$TAB \
	bound_sites+=$(zcat $CELL6_PEAK_FILE | wc -l) \
	bound_sites+=$TAB \
	echo 'Intersecting with peaks' \
	CELL1_SCANNED_MOTIF_BED_INTERSECT_PEAKS=$CENTIPEDE_DATA/${CELL1_TITLE}.${motif_info}.intersect_peak.bed.gz \
	CELL1_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY=$CENTIPEDE_SINGULARITY_DATA/${CELL1_TITLE}.${motif_info}.intersect_peak.bed.gz \
	intersectBed -a $TARGET_SCANNED_MOTIF_BED -b $CELL1_PEAK_FILE -wa | gzip > $CELL1_SCANNED_MOTIF_BED_INTERSECT_PEAKS \
	bound_sites+=$(zcat $CELL1_SCANNED_MOTIF_BED_INTERSECT_PEAKS | wc -l) \
	bound_sites+=$TAB \
	CELL2_SCANNED_MOTIF_BED_INTERSECT_PEAKS=$CENTIPEDE_DATA/${CELL2_TITLE}.${motif_info}.intersect_peak.bed.gz \
	CELL2_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY=$CENTIPEDE_SINGULARITY_DATA/${CELL2_TITLE}.${motif_info}.intersect_peak.bed.gz \
	intersectBed -a $TARGET_SCANNED_MOTIF_BED -b $CELL2_PEAK_FILE -wa | gzip > $CELL2_SCANNED_MOTIF_BED_INTERSECT_PEAKS \
	bound_sites+=$(zcat $CELL2_SCANNED_MOTIF_BED_INTERSECT_PEAKS | wc -l) \
	bound_sites+=$TAB \
	CELL3_SCANNED_MOTIF_BED_INTERSECT_PEAKS=$CENTIPEDE_DATA/${CELL3_TITLE}.${motif_info}.intersect_peak.bed.gz \
	CELL3_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY=$CENTIPEDE_SINGULARITY_DATA/${CELL3_TITLE}.${motif_info}.intersect_peak.bed.gz \
	intersectBed -a $TARGET_SCANNED_MOTIF_BED -b $CELL3_PEAK_FILE -wa | gzip > $CELL3_SCANNED_MOTIF_BED_INTERSECT_PEAKS \
	bound_sites+=$(zcat $CELL3_SCANNED_MOTIF_BED_INTERSECT_PEAKS | wc -l) \
	bound_sites+=$TAB \
	CELL4_SCANNED_MOTIF_BED_INTERSECT_PEAKS=$CENTIPEDE_DATA/${CELL4_TITLE}.${motif_info}.intersect_peak.bed.gz \
	CELL4_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY=$CENTIPEDE_SINGULARITY_DATA/${CELL4_TITLE}.${motif_info}.intersect_peak.bed.gz \
	intersectBed -a $TARGET_SCANNED_MOTIF_BED -b $CELL4_PEAK_FILE -wa | gzip > $CELL4_SCANNED_MOTIF_BED_INTERSECT_PEAKS \
	bound_sites+=$(zcat $CELL4_SCANNED_MOTIF_BED_INTERSECT_PEAKS | wc -l) \
	bound_sites+=$TAB \
	CELL5_SCANNED_MOTIF_BED_INTERSECT_PEAKS=$CENTIPEDE_DATA/${CELL5_TITLE}.${motif_info}.intersect_peak.bed.gz \
	CELL5_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY=$CENTIPEDE_SINGULARITY_DATA/${CELL5_TITLE}.${motif_info}.intersect_peak.bed.gz \
	intersectBed -a $TARGET_SCANNED_MOTIF_BED -b $CELL5_PEAK_FILE -wa | gzip > $CELL5_SCANNED_MOTIF_BED_INTERSECT_PEAKS \
	bound_sites+=$(zcat $CELL5_SCANNED_MOTIF_BED_INTERSECT_PEAKS | wc -l) \
	bound_sites+=$TAB \
	CELL6_SCANNED_MOTIF_BED_INTERSECT_PEAKS=$CENTIPEDE_DATA/${CELL6_TITLE}.${motif_info}.intersect_peak.bed.gz \
	CELL6_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY=$CENTIPEDE_SINGULARITY_DATA/${CELL6_TITLE}.${motif_info}.intersect_peak.bed.gz \
	intersectBed -a $TARGET_SCANNED_MOTIF_BED -b $CELL6_PEAK_FILE -wa | gzip > $CELL6_SCANNED_MOTIF_BED_INTERSECT_PEAKS \
	bound_sites+=$(zcat $CELL6_SCANNED_MOTIF_BED_INTERSECT_PEAKS | wc -l) \
	bound_sites+=$TAB \
	CELL1_CUT_MATRIX=$CENTIPEDE_DATA/${CELL1_TITLE}.${motif_info}.matrix.gz \
	singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL1_BAM_FILE_SINGULARITY $CELL1_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY | gzip > $CELL1_CUT_MATRIX \
	CELL2_CUT_MATRIX=$CENTIPEDE_DATA/${CELL2_TITLE}.${motif_info}.matrix.gz \
	singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL2_BAM_FILE_SINGULARITY $CELL2_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY | gzip > $CELL2_CUT_MATRIX \
	CELL3_CUT_MATRIX=$CENTIPEDE_DATA/${CELL3_TITLE}.${motif_info}.matrix.gz \
	singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL3_BAM_FILE_SINGULARITY $CELL3_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY | gzip > $CELL3_CUT_MATRIX \
	CELL4_CUT_MATRIX=$CENTIPEDE_DATA/${CELL4_TITLE}.${motif_info}.matrix.gz \
	singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL4_BAM_FILE_SINGULARITY $CELL4_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY | gzip > $CELL4_CUT_MATRIX \
	CELL5_CUT_MATRIX=$CENTIPEDE_DATA/${CELL5_TITLE}.${motif_info}.matrix.gz \
	singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL5_BAM_FILE_SINGULARITY $CELL5_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY | gzip > $CELL5_CUT_MATRIX \
	CELL6_CUT_MATRIX=$CENTIPEDE_DATA/${CELL6_TITLE}.${motif_info}.matrix.gz \
	singularity exec $CENTEPEDE_IMAGE make_cut_matrix -d -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 150-324 325-400 1)' $CELL6_BAM_FILE_SINGULARITY $CELL6_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY | gzip > $CELL6_CUT_MATRIX \
	CELL1_CUT_MATRIX=$CENTIPEDE_SINGULARITY_DATA/${CELL1_TITLE}.${motif_info}.matrix.gz \
	CELL1_CENTEPEDE_DISCREET_TABLE_SINGULARITY=$CENTIPEDE_SINGULARITY_DATA/${CELL1_TITLE}.${motif_info}.motif_centepede.bed.gz \
	CELL1_CENTEPEDE_DISCREET_TABLE=$CENTIPEDE_DATA/${CELL1_TITLE}.${motif_info}.motif_centepede.bed.gz \
	CELL1_CENTEPEDE_PDF=$CENTIPEDE_SINGULARITY_DATA/${CELL1_TITLE}.${motif_info}.motif_centepede.pdf \
	singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL1_CUT_MATRIX $CELL1_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY $CELL1_CENTEPEDE_DISCREET_TABLE_SINGULARITY $CELL1_CENTEPEDE_PDF $motif_length \
	bound_sites+=$(zcat $CELL1_CENTEPEDE_DISCREET_TABLE| wc -l) \
	bound_sites+=$TAB \
	CELL2_CUT_MATRIX=$CENTIPEDE_SINGULARITY_DATA/${CELL2_TITLE}.${motif_info}.matrix.gz \
	CELL2_CENTEPEDE_DISCREET_TABLE_SINGULARITY=$CENTIPEDE_SINGULARITY_DATA/${CELL2_TITLE}.${motif_info}.motif_centepede.bed.gz \
	CELL2_CENTEPEDE_DISCREET_TABLE=$CENTIPEDE_DATA/${CELL2_TITLE}.${motif_info}.motif_centepede.bed.gz \
	CELL2_CENTEPEDE_PDF=$CENTIPEDE_SINGULARITY_DATA/${CELL2_TITLE}.${motif_info}.motif_centepede.pdf \
	singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL2_CUT_MATRIX $CELL2_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY $CELL2_CENTEPEDE_DISCREET_TABLE_SINGULARITY $CELL2_CENTEPEDE_PDF $motif_length \
	bound_sites+=$(zcat $CELL2_CENTEPEDE_DISCREET_TABLE | wc -l) \
	bound_sites+=$TAB \
	CELL3_CUT_MATRIX=$CENTIPEDE_SINGULARITY_DATA/${CELL3_TITLE}.${motif_info}.matrix.gz \
	CELL3_CENTEPEDE_DISCREET_TABLE_SINGULARITY=$CENTIPEDE_SINGULARITY_DATA/${CELL3_TITLE}.${motif_info}.motif_centepede.bed.gz \
	CELL3_CENTEPEDE_DISCREET_TABLE=$CENTIPEDE_DATA/${CELL3_TITLE}.${motif_info}.motif_centepede.bed.gz \
	CELL3_CENTEPEDE_PDF=$CENTIPEDE_SINGULARITY_DATA/${CELL3_TITLE}.${motif_info}.motif_centepede.pdf \
	singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL3_CUT_MATRIX $CELL3_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY $CELL3_CENTEPEDE_DISCREET_TABLE_SINGULARITY $CELL3_CENTEPEDE_PDF $motif_length \
	bound_sites+=$(zcat $CELL3_CENTEPEDE_DISCREET_TABLE | wc -l) \
	bound_sites+=$TAB \
	CELL4_CUT_MATRIX=$CENTIPEDE_SINGULARITY_DATA/${CELL4_TITLE}.${motif_info}.matrix.gz \
	CELL4_CENTEPEDE_DISCREET_TABLE_SINGULARITY=$CENTIPEDE_SINGULARITY_DATA/${CELL4_TITLE}.${motif_info}.motif_centepede.bed.gz \
	CELL4_CENTEPEDE_DISCREET_TABLE=$CENTIPEDE_DATA/${CELL4_TITLE}.${motif_info}.motif_centepede.bed.gz \
	CELL4_CENTEPEDE_PDF=$CENTIPEDE_SINGULARITY_DATA/${CELL4_TITLE}.${motif_info}.motif_centepede.pdf \
	singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL4_CUT_MATRIX $CELL4_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY $CELL4_CENTEPEDE_DISCREET_TABLE_SINGULARITY $CELL4_CENTEPEDE_PDF $motif_length \
	bound_sites+=$(zcat $CELL4_CENTEPEDE_DISCREET_TABLE | wc -l) \
	bound_sites+=$TAB \
	CELL5_CUT_MATRIX=$CENTIPEDE_SINGULARITY_DATA/${CELL5_TITLE}.${motif_info}.matrix.gz \
	CELL5_CENTEPEDE_DISCREET_TABLE_SINGULARITY=$CENTIPEDE_SINGULARITY_DATA/${CELL5_TITLE}.${motif_info}.motif_centepede.bed.gz \
	CELL5_CENTEPEDE_DISCREET_TABLE=$CENTIPEDE_DATA/${CELL5_TITLE}.${motif_info}.motif_centepede.bed.gz \
	CELL5_CENTEPEDE_PDF=$CENTIPEDE_SINGULARITY_DATA/${CELL5_TITLE}.${motif_info}.motif_centepede.pdf \
	singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL5_CUT_MATRIX $CELL5_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY $CELL5_CENTEPEDE_DISCREET_TABLE_SINGULARITY $CELL5_CENTEPEDE_PDF $motif_length \
	bound_sites+=$(zcat $CELL5_CENTEPEDE_DISCREET_TABLE | wc -l) \
	bound_sites+=$TAB \
	CELL6_CUT_MATRIX=$CENTIPEDE_SINGULARITY_DATA/${CELL6_TITLE}.${motif_info}.matrix.gz \
	CELL6_CENTEPEDE_DISCREET_TABLE_SINGULARITY=$CENTIPEDE_SINGULARITY_DATA/${CELL6_TITLE}.${motif_info}.motif_centepede.bed.gz \
	CELL6_CENTEPEDE_DISCREET_TABLE=$CENTIPEDE_DATA/${CELL6_TITLE}.${motif_info}.motif_centepede.bed.gz \
	CELL6_CENTEPEDE_PDF=$CENTIPEDE_SINGULARITY_DATA/${CELL6_TITLE}.${motif_info}.motif_centepede.pdf \
	singularity exec $CENTEPEDE_IMAGE run_centepede.R $CELL6_CUT_MATRIX $CELL6_SCANNED_MOTIF_BED_INTERSECT_PEAKS_SINGULARITY $CELL6_CENTEPEDE_DISCREET_TABLE_SINGULARITY $CELL6_CENTEPEDE_PDF $motif_length \
	bound_sites+=$(zcat $CELL6_CENTEPEDE_DISCREET_TABLE | wc -l) \
	bound_sites+=$TAB \
	CELL1_CENTEPEDE_BOUND_TABLE=$CENTIPEDE_DATA/${CELL1_TITLE}.${motif_info}.motif_centepede.bound.bed.gz \	
	zcat $CELL1_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL1_CENTEPEDE_BOUND_TABLE \
	bound_sites+=$(zcat $CELL1_CENTEPEDE_BOUND_TABLE | wc -l) \
	bound_sites+=$TAB \
	CELL2_CENTEPEDE_BOUND_TABLE=$CENTIPEDE_DATA/${CELL2_TITLE}.${motif_info}.motif_centepede.bound.bed.gz \
	zcat $CELL2_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL2_CENTEPEDE_BOUND_TABLE \
	bound_sites+=$(zcat $CELL2_CENTEPEDE_BOUND_TABLE | wc -l) \
	bound_sites+=$TAB \
	CELL3_CENTEPEDE_BOUND_TABLE=$CENTIPEDE_DATA/${CELL3_TITLE}.${motif_info}.motif_centepede.bound.bed.gz \
	zcat $CELL3_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL3_CENTEPEDE_BOUND_TABLE \
	bound_sites+=$(zcat $CELL3_CENTEPEDE_BOUND_TABLE | wc -l) \
	bound_sites+=$TAB \
	CELL4_CENTEPEDE_BOUND_TABLE=$CENTIPEDE_DATA/${CELL4_TITLE}.${motif_info}.motif_centepede.bound.bed.gz \
	zcat $CELL4_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL4_CENTEPEDE_BOUND_TABLE \
	bound_sites+=$(zcat $CELL4_CENTEPEDE_BOUND_TABLE | wc -l) \
	bound_sites+=$TAB \
	CELL5_CENTEPEDE_BOUND_TABLE=$CENTIPEDE_DATA/${CELL5_TITLE}.${motif_info}.motif_centepede.bound.bed.gz \
	zcat $CELL5_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL5_CENTEPEDE_BOUND_TABLE \
	bound_sites+=$(zcat $CELL5_CENTEPEDE_BOUND_TABLE | wc -l) \
	bound_sites+=$TAB \
	CELL6_CENTEPEDE_BOUND_TABLE=$CENTIPEDE_DATA/${CELL6_TITLE}.${motif_info}.motif_centepede.bound.bed.gz \
	
	zcat $CELL6_CENTEPEDE_DISCREET_TABLE | awk -F"\t" -v OFS="\t" '$7 >= "'$BOUND_POSTERIOR_PROBS_LIMIT'" { print $1,$2,$3,$4,$5,$6 }' | gzip -c > $CELL6_CENTEPEDE_BOUND_TABLE \
	bound_sites+=$(zcat $CELL6_CENTEPEDE_BOUND_TABLE | wc -l) \

	CELL1_CENTEPEDE_BOUND_TABLE=$CENTIPEDE_SINGULARITY_DATA/${CELL1_TITLE}.${motif_info}.motif_centepede.bound.bed.gz \
	CELL1_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$CENTIPEDE_DATA/${CELL1_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz \
	singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL1_BAM_FILE_SINGULARITY $CELL1_CENTEPEDE_BOUND_TABLE | gzip > $CELL1_CENTEPEDE_BOUND_AGGREGATE_MATRIX \


	CELL2_CENTEPEDE_BOUND_TABLE=$CENTIPEDE_SINGULARITY_DATA/${CELL2_TITLE}.${motif_info}.motif_centepede.bound.bed.gz \
	CELL2_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$CENTIPEDE_DATA/${CELL2_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz \
	singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL2_BAM_FILE_SINGULARITY $CELL2_CENTEPEDE_BOUND_TABLE | gzip > $CELL2_CENTEPEDE_BOUND_AGGREGATE_MATRIX \


	CELL3_CENTEPEDE_BOUND_TABLE=$CENTIPEDE_SINGULARITY_DATA/${CELL3_TITLE}.${motif_info}.motif_centepede.bound.bed.gz \
	CELL3_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$CENTIPEDE_DATA/${CELL3_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz \
	singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL3_BAM_FILE_SINGULARITY $CELL3_CENTEPEDE_BOUND_TABLE | gzip > $CELL3_CENTEPEDE_BOUND_AGGREGATE_MATRIX \


	CELL4_CENTEPEDE_BOUND_TABLE=$CENTIPEDE_SINGULARITY_DATA/${CELL4_TITLE}.${motif_info}.motif_centepede.bound.bed.gz \
	CELL4_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$CENTIPEDE_DATA/${CELL4_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz \
	singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL4_BAM_FILE_SINGULARITY $CELL4_CENTEPEDE_BOUND_TABLE | gzip > $CELL4_CENTEPEDE_BOUND_AGGREGATE_MATRIX \

	CELL5_CENTEPEDE_BOUND_TABLE=$CENTIPEDE_SINGULARITY_DATA/${CELL5_TITLE}.${motif_info}.motif_centepede.bound.bed.gz \
	CELL5_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$CENTIPEDE_DATA/${CELL5_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz \
	singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL5_BAM_FILE_SINGULARITY $CELL5_CENTEPEDE_BOUND_TABLE | gzip > $CELL5_CENTEPEDE_BOUND_AGGREGATE_MATRIX \


	CELL6_CENTEPEDE_BOUND_TABLE=$CENTIPEDE_SINGULARITY_DATA/${CELL6_TITLE}.${motif_info}.motif_centepede.bound.bed.gz \
	CELL6_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$CENTIPEDE_DATA/${CELL6_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz \
	singularity exec $CENTEPEDE_IMAGE make_cut_matrix -a -p 2 -f 3 -F 4 -F 8 -q 30 -b '(36-149 1) (150-324 2) (325-400 5)' $CELL6_BAM_FILE_SINGULARITY $CELL6_CENTEPEDE_BOUND_TABLE | gzip > $CELL6_CENTEPEDE_BOUND_AGGREGATE_MATRIX \

	CELL1_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$CENTIPEDE_SINGULARITY_DATA/${CELL1_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz \
	CELL2_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$CENTIPEDE_SINGULARITY_DATA/${CELL2_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz \
	CELL3_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$CENTIPEDE_SINGULARITY_DATA/${CELL3_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz \
	CELL4_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$CENTIPEDE_SINGULARITY_DATA/${CELL4_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz \
	CELL5_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$CENTIPEDE_SINGULARITY_DATA/${CELL5_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz \
	CELL6_CENTEPEDE_BOUND_AGGREGATE_MATRIX=$CENTIPEDE_SINGULARITY_DATA/${CELL6_TITLE}.${motif_info}.motif_centepede.bound.agg.matrix.gz \

	CENTEPEDE_MOTIF_AGGREGATE_PDF=$CENTIPEDE_SINGULARITY_DATA/${CELL1_TITLE}.${CELL2_TITLE}.${CELL3_TITLE}.${CELL4_TITLE}.${CELL5_TITLE}.${CELL6_TITLE}.${motif_info}.motif_aggregate_centepede.pdf \

	singularity exec $CENTEPEDE_IMAGE $WORKDIR_SINGULARITY/plot_aggregate_matrix_new.R $CELL1_CENTEPEDE_BOUND_AGGREGATE_MATRIX $CELL2_CENTEPEDE_BOUND_AGGREGATE_MATRIX $CELL3_CENTEPEDE_BOUND_AGGREGATE_MATRIX $CELL4_CENTEPEDE_BOUND_AGGREGATE_MATRIX $CELL5_CENTEPEDE_BOUND_AGGREGATE_MATRIX $CELL6_CENTEPEDE_BOUND_AGGREGATE_MATRIX '${motif_title} regions with motifs oriented by strand' $CENTEPEDE_MOTIF_AGGREGATE_PDF \
	echo '$bound_sites' >> $BOUND_SITES_REPORT \
	" >> $SWARM_FILE
done


echo "Done!!! "$(date)


