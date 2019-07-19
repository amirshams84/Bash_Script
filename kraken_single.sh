#! /bin/bash
set -o pipefail
set -e

PROCESSORS=$((SLURM_CPUS_PER_TASK - 2))

module load kraken
module load usearch

STANDARD_DB=/fdb/kraken/20180220_standard

INPUTDIR=/data/shamsaddinisha/WGS/fastq
OUTPUTDIR=/data/shamsaddinisha/WGS/kraken
DATA_FORMAT=fastq
for sample in $INPUTDIR/*.${DATA_FORMAT}
do
	echo $sample
	sample_name=$(basename $sample)
	echo "Running kraken on $sample_name"
	echo "kraken --db $STANDARD_DB --threads $PROCESSORS --${DATA_FORMAT}-input --out-fmt legacy --classified-out $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.classified \
	--unclassified-out $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.unclassified --output $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken \
	$sample > $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken_report.txt 2> $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken_error_report.txt"

	kraken --db $STANDARD_DB --threads $PROCESSORS --${DATA_FORMAT}-input --out-fmt legacy --classified-out $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.classified.fasta \
	--unclassified-out $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.unclassified.fasta --output $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken \
	$sample > $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken_report.txt 2> $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken_error_report.txt
	echo "Interpret kraken result of $sample_name"
	echo "kraken-translate --db $STANDARD_DB $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken > $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.translated.kraken.txt"
	kraken-translate --db $STANDARD_DB $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken > $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.translated.kraken.txt

	usearch -fastx_getlabels $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.classified.fasta -output $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.classified_labels.txt > $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken_report.txt 2> $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken_error_report.txt
	usearch -fastx_getseqs $sample -labels $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.classified_labels.txt -${DATA_FORMAT}out $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.classified.fastq > $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken_report.txt 2> $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken_error_report.txt
	usearch -fastx_getlabels $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.unclassified.fasta -output $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.unclassified_labels.txt > $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken_report.txt 2> $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken_error_report.txt
	usearch -fastx_getseqs $sample -labels $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.unclassified_labels.txt -${DATA_FORMAT}out $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.unclassified.fastq > $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken_report.txt 2> $OUTPUTDIR/${sample_name%.$DATA_FORMAT}.kraken_error_report.txt
done


echo "Pipeline Successfully finished at: "$(date)