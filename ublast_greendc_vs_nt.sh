#!/bin/bash
module load usearch
usearch -ublast /data/shamsaddinisha/GREEN_DC/GREENDC_total_merged_relabeled.fasta -db /data/shamsaddinisha/nt_fasta/nt_ublast.udb -strand both -evalue 1e-9 -maxhits 1 -uc /data/shamsaddinisha/GREEN_DC/GREEN_DC_UBLAST.txt
