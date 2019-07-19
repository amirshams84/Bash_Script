#!/bin/bash	
module load spades
spades.py -o ./hpminus_scaffold/ -k 21,33,55,77,99,127 --phred-offset 33 --only-assembler --memory 500 -t $SLURM_CPUS_PER_TASK  -s hpminus_total.fastq
