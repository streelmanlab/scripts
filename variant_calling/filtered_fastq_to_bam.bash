#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Create a bam file from the filtered fastq files

#Usage statement
usage () {
        echo "Usage: filtered_fastq_to_bam.bash filtered_fastq_1 filtered_fastq_2 reference outputFile [options]
		-N Name of Job
		-l memnory
		-t hh:mm:ss time needed, job will be killed if exceeded
		-q specifies process queue
		-j controls what gets written to the ouputfile
		-o name of the job's outputfile
		-m controls when email is sent to the submitter
		-M email of submitter"
}

#Command-line options
get_input() {
	fastq1=$1
	fastq2=$2
	ref=$3
	out=$4
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
module load open64/4.5.1
module load bwa/0.7.4     #loads bwa package
module load samtools      #   loads samtools package
bwa mem -M $ref $fastq1 $fastq2 > $out"".sam     #creates the SAM
samtools view -bS $out"".sam > $out       #Converts to BAM" > filtered_fastq_to_bam.pbs
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	qsub filtered_fastq_to_bam.pbs
}


main "$@"
