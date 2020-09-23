#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Create a bam file from the filtered fastq files using STAR (splice aware)

#Usage statement
usage () {
        echo "Usage: star.bash filtered_fastq_1 filtered_fastq_2 reference refDir outputiPrefix [options]
		-N Name of Job
		-l memory
		-t hh:mm:ss time needed, job will be killed if exceeded
		-q specifies process queue
		-j controls what gets written to the ouputfile
		-o name of the job's outputfile
		-m controls when email is sent to the submitter
		-M email of submitter"
}

check_for_help() {
	if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
		usage;
		exit 0;
	fi
}

#Command-line options
get_input() {
	check_for_help "$1"
	
	fastq1=$1
	fastq2=$2
	ref=$3
	gdir=$4
	out=$5
	shift
	shift
	shift
	shift
	shift
	name="$out"
	memory="mem=128gb"
	time="walltime=40:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputFile="$out"".out"
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
	if [ -z "$fastq1" ] || [ ! -f "$fastq1" ]; then
		echo "The first fastq file does not exist. Exiting the program."
		usage
		exit 1
	fi
	
	if [ -z "$fastq2" ] || [ ! -f "$fastq2" ]; then
		echo "The second fastq file does not exist. Exiting the program."
		usage
		exit 1
	fi

	if [ -z "$ref" ] || [ ! -f "$ref" ]; then
                echo "The reference file does not exist. Exiting the program."
                usage
                exit 1
        fi
	
	if [ -z "$gdir" ]; then
		echo "Please specify a STAR reference directory."
		usage
		exit 1	
	fi

	createRefStr=""
	if [ -z "$(ls -A $gdir)" ]; then
		createRefStr="STAR --runMode genomeGenerate --genomeDir $gdir --genomeFastaFiles $ref"
	fi
	
	if [ -z "$out" ]; then
		echo "Please specify an output file."
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
module load gcc/4.9.0
module load STAR/2.5.3a
module load samtools      #   loads samtools package

# Create reference
$createRefStr

# Align Reads
STAR --genomeDir $gdir --readFilesIn $fastq1 --readFilesIn $fastq2 --outFileNamePrefix $out --outSAMtype BAM SortedByCoordinate --clip5pNbases 6
samtools view -bS $out"".sam > $out"".bam       #Converts to BAM" > star.pbs
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	qsub star.pbs
}


main "$@"
