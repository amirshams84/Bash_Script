#! /bin/bash
set -o pipefail
set -e
module load java/1.8.0_92 || exit 1

echo "Running on $SLURM_CPUS_PER_TASK CPUs"

FASTA_QUERY=/data/shamsaddinisha/MBAC/Trumpeter_Swan_Project/SnipCalling/Reference/Trumpeter_corrected_contigs.fasta
WORK_DIR=/data/shamsaddinisha/MBAC/Trumpeter_Swan_Project/InterproScan/Result
TEMP_DIR=/data/shamsaddinisha/MBAC/Trumpeter_Swan_Project/InterproScan/Temp
TITLE=Trumpeter_corrected_contigs

/usr/local/apps/interproscan/5.22-61.0/interproscan_app/interproscan.sh --seqtype n --goterms --pathways -f svg,tsv,html,gff3 -iprlookup \
--tempdir $TEMP_DIR --disable-precalc --input $FASTA_QUERY --output-dir $WORK_DIR
