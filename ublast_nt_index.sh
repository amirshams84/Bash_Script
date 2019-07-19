#!/bin/bash
module load usearch
usearch -makeudb_ublast /data/shamsaddinisha/nt_fasta/nt.fasta -wordlength 14 -dbstep 14 -output /data/shamsaddinisha/nt_fasta/nt_ublast.udb
