#! /bin/bash
set -o pipefail
set -e

PROCESSORS=$((SLURM_CPUS_PER_TASK - 2))

module load kraken

STANDARD_DB=/fdb/kraken/20180220_standard

INPUTDIR=/data/shamsaddinisha/WGS/fastq
OUTPUTDIR=/data/shamsaddinisha/WGS/kraken
DATA_FORMAT=fastq

for sample in $INPUTDIR/*_R1.${DATA_FORMAT}
do
	echo $sample
	sample_name=$(basename $sample)
	echo "Running kraken on $sample_name"
	echo "kraken --db $STANDARD_DB --threads $PROCESSORS --${DATA_FORMAT}-input --out-fmt paired --${DATA_FORMAT}-output --classified-out $OUTPUTDIR/${sample_name%_R1*}.classified \
	--unclassified-out $OUTPUTDIR/${sample_name%_R1*}.unclassified --output $OUTPUTDIR/${sample_name%_R1*}.kraken --paired \
	$sample ${sample%_R1*}_R2.${DATA_FORMAT} > $OUTPUTDIR/${sample_name%_R1*}.kraken_report.txt 2> $OUTPUTDIR/${sample_name%_R1*}.kraken_error_report.txt"

	kraken --db $STANDARD_DB --threads $PROCESSORS --${DATA_FORMAT}-input --out-fmt paired --${DATA_FORMAT}-output --classified-out $OUTPUTDIR/${sample_name%_R1*}.classified \
	--unclassified-out $OUTPUTDIR/${sample_name%_R1*}.unclassified --output $OUTPUTDIR/${sample_name%_R1*}.kraken --paired \
	$sample ${sample%_R1*}_R2.${DATA_FORMAT} > $OUTPUTDIR/${sample_name%_R1*}.kraken_report.txt 2> $OUTPUTDIR/${sample_name%_R1*}.kraken_error_report.txt
	echo "Interpret kraken result of $sample_name"
	echo "kraken-translate --db $STANDARD_DB $OUTPUTDIR/${sample_name%_R1*}.kraken > $OUTPUTDIR/${sample_name%_R1*}.translated.kraken.txt"
	kraken-translate --db $STANDARD_DB $OUTPUTDIR/${sample_name%_R1*}.kraken > $OUTPUTDIR/${sample_name%_R1*}.translated.kraken.txt
done


echo "Pipeline Successfully finished at: "$(date)