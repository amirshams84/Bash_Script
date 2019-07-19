#!/bin/bash

DATA_DIR="/data/dadkhahe/00-BLI/01_DemultiplexedFastqs"

module load mixcr


for fastq in ${DATA_DIR}/*R1*.fastq*
do
	fastq_name=$(basename $fastq)
	echo $fastq_name
	sample_name=${fastq_name%L001*}
	echo $sample_name
	forward_fastq=$fastq_name
	echo $forward_fastqs
	reverse_fastq=${fastq_name/R1/R2}
	echo $reverse_fastq

	mixcr analyze shotgun \
	--species hs \
	--starting-material rna \
	--only-productive \
	$DATA_DIR/$forward_fastq $DATA_DIR/$reverse_fastq ${sample_name}_analysis
	
done