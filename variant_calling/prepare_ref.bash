#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Prepare the reference file

#Usage statement
usage () {
        echo "Usage: prepare_ref.bash reference gatk_file_path [options]
		-N Name of Job
		-l memnory
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
	
	ref=$1
	gatk=$2
	shift
	shift
	
	name="prepare_ref"
	memory="mem=8gb"
	time="walltime=40:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputFile="prepare_ref.out"
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
	if [ -z "$ref" ] || [ ! -f "$ref" ]; then
                echo "The reference file does not exist. Exiting the program."
                usage
                exit 1
        fi
	
	if [ -z "$gatk" ] || [ ! -f "$gatk" ]; then
                echo "The gatk file does not exist. Please give the full file path to the gatk binary (include the gatk binary in the path). Exiting the program."
                usage
                exit 1
        fi
}

generate_pbs() {
	echo "#PBS -A GT-js585-biocluster
#PBS -N $name
#PBS -l $memory
#PBS -l nodes=2:ppn=4
#PBS -l $time
#PBS -j $writingOpts
#PBS -o $outputFile
#PBS -m $emailOpts
#PBS -M $email

cd \$PBS_O_WORKDIR
$gatk CreateSequenceDictionary -R $ref

samtools faidx $ref

module load open64/4.5.1
module load bwa/0.7.17
bwa index -a bwtsw $ref" > prepare_ref.pbs
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	qsub prepare_ref.pbs
}


main "$@"
