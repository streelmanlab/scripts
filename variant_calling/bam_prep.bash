#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Prepare the bam file for SNP Calling by sorting the BAM file, marking duplicates, adding or replacing read groups, and building a BAM index.

#Usage statement
usage () {
        echo "Usage: bam_prep.bash bamFile picard_file_path [options]
		bamFile: Name of bam file excluding the .bam file extension
		picard_file_path: Path to picard.jar file
		-N Name of Job
		-l memory
		-t hh:mm:ss time needed, job will be killed if exceeded
		-q specifies process queue
		-j controls what gets written to the ouputfile
		-o name of the job's outputfile
		-m controls when email is sent to the submitter
		-M email of submitter"
}

#Command-line options
get_input() {
	bam=$1
	picard=$2
	shift
	shift
	
	name="$bam"
	memory="mem=128gb"
	time="walltime=8:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputFile="$bam"".out"
	emailOpts="abe"
	email="ggruenhagen3@gatech.edu"
	while getopts "N:l:t:q:j:o:m:M:h" opt; do
		case $opt in
		N ) name=$OPTARG ;;
		l ) memory=$OPTARG ;;
		t ) time=$OPTARG ;;
		q ) cluster=$OPTARG ;;
		j ) writingOpts=$OPTARG ;;
		o ) outputFile=$OPTARG ;;
		m ) emailOpts=$OPTARG ;;
		M ) email=$OPTARG ;;
		h ) usage; exit 0;
		esac
	done
}

check_files() {
	if [ -z "$bam"".bam" ] || [ ! -f "$bam"".bam" ]; then
		echo "The first fastq file does not exist. Exiting the program."
		usage
		exit 1
	fi
	if [ -z "$picard" ] || [ ! -f "$picard" ] || [ ${picard: -4} != ".jar" ]; then
		echo "Not a valid path to a picard.jar file. Exiting the program."
		usage
		exit 1
	fi
}

generate_pbs() {
	echo "#PBS -N $name
#PBS -l $memory
#PBS -l nodes=2:ppn=4
#PBS -l $time
#PBS -q $cluster
#PBS -j $writingOpts
#PBS -o $outputFile
#PBS -m $emailOpts
#PBS -M $email

cd \$PBS_O_WORKDIR
module purge
module load java

if [ ! -d ./tmp/ ]; then
	mkdir -p ./tmp/;
fi
java -jar $picard SortSam I=$bam"".bam O=$bam""_sorted.bam SORT_ORDER=coordinate TMP_DIR= ./tmp/

java -jar $picard MarkDuplicates I=$bam""_sorted.bam O=$bam""_dedup.bam METRICS_FILE=$bam""_metrics.txt TMP_DIR= /tmp/

java -jar $picard AddOrReplaceReadGroups I=$bam""_dedup.bam O=$bam""_RG.bam RGID=$bam  RGLB=lib$bam RGPL=Illumina RGPU=unit$bam RGSM=$bam

java -jar $picard BuildBamIndex I=$bam""_RG.bam" > bam_prep.pbs
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	qsub bam_prep.pbs
}


main "$@"
