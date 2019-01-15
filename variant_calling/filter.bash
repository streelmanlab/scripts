#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Filter a vcf file

#Usage statement
usage () {
        echo "Usage: bash filter.bash <vcf_file> [options]
		
		vcf_file: vcf file to be filtered
		
		-P percent to filter by (default: 50)
		-O name of output vcf file (default: filtered.vcf)
		-N Name of Job
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
	
	vcf=$1
	shift
	
	perct="50"
	out="filtered.vcf"
	name="filter"
	memory="mem=128gb"
	time="walltime=40:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputFile="filter.out"
	emailOpts="abe"
	email="ggruenhagen3@gatech.edu"
	while getopts "P:O:N:l:t:q:j:o:m:M:h" opt; do
		case $opt in
		P ) perct=$OPTARG ;;
		O ) out=$OPTARG ;;
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
	if [ -z "$vcf" ]; then
		echo "The vcf file does not exist: $vcf"
		usage
		exit 1
	fi
}


generate_pbs() {
	echo "#PBS -N $name
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

bgzip $vcf
tabix -p vcf $vcf"".gz

vcftools --gzvcf $vcf"".gz --max-missing 0.""$perct --recode --out filter_$perct
" > filter.pbs
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	qsub filter.pbs
}


main "$@"
