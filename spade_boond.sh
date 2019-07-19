#!/bin/bash	
module load spades
spades.py -o ./BOOND_contigs/ -k 21,33,55,77,99,127 --phred-offset 33 --only-assembler --memory 500 -t $SLURM_CPUS_PER_TASK  -s ./BOOND_fastq/A7544_EE_N_filtered.fastq
