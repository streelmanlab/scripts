#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Prepare the reference file

#Usage statement
usage () {
        echo "Usage: merge.bash <vcf_directory> <bgzip_path> <tabix_path> [options]
		
		vcf_directory: directory containing the vcf files to be merged
		bgzip_path: path to bgzip binary (links to downloads in the README)
		tabix_path: path to tabix binary (links to downloads in the README)
		
		-N Name of Job
		-l memnory
		-t hh:mm:ss time needed, job will be killed if exceeded (default walltime=40:00:00)
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
	
	dir=$1
	bg=$2
	tab=$3
	shift
	shift
	shift
	
	name="merge"
	memory="mem=128gb"
	time="walltime=40:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputFile="merge.out"
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

zip() {
	zipStr=""
	for file in $dir/*.vcf; do
		if [ -f "$file" ]; then
			zipStr+="$bg $file"
		fi
	done
}

tabix() {
	tabixStr=""
	files=""
	for file in $dir/*.vcf; do
		if [ -f "$file" ]; then
			tabixStr+="$tab -p vcf $file"".gz\n"
			files+="$file"".gz "
		fi
	done
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
module load intel/14.0.2
module load perl/5.14.2
module load samtools/0.1.18
module load vcftools/0.1.14.10

$tabixStr

vcf-merge $files" > merge.pbs
}

main() {
	get_input "$@"
	zip
	tabix
	generate_pbs
	
	qsub merge.pbs
}


main "$@"
