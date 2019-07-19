#!/bin/bash	
module load spades
spades.py -o /data/shamsaddinisha/GREEN_DC/AB1850_contigs/ -k 21,33,55,77,99,127 --phred-offset 33 --only-assembler --memory 700 -t $SLURM_CPUS_PER_TASK  -s /data/shamsaddinisha/GREEN_DC/AB1850_filtered.fastq
