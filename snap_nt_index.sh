#!/bin/bash
module load snap/1.0beta.11
snap index /data/shamsaddinisha/nt_acc/nt_acc.fasta /data/shamsaddinisha/nt_acc/nt_acc_snap_index/ -t$SLURM_CPUS_PER_TASK
