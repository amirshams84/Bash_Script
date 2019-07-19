#!/bin/bash
module load usearch
usearch -makeudb_ublast /data/shamsaddinisha/refseq/bacteria/Genomic/fasta/refseq_bacteria.fasta -wordlength 14 -dbstep 14 -output /data/shamsaddinisha/nt_fasta/refseq_bacteria_ublast.udb
