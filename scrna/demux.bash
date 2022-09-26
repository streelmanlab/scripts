#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Use demuxlet to demultiplex individuals

#Usage statement
usage () {
        echo "Usage: bash demux.bash <cell_vcf> <gt_vcf> [options]
		
		cell_vcf: vcf file originating from reads in cells/nuclei
		gt_vcf:   vcf file originating from genotyped individuals
		
		-O basename for output demuxlet files (default: demux_out)
		-N Name of Job
		-t hh:mm:ss time needed, job will be killed if exceeded (default walltime: 1:00:00)
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
	
	cell_vcf=$1
	gt_vcf=$2
	shift
	shift

	out="demux_out"
	name="demux"
	memory="mem=64gb"
	time="walltime=1:00:00"
	writingOpts="oe"
	outputFile="demux_out.out"
	emailOpts="abe"
	email="ggruenhagen3@gatech.edu"
	while getopts "O:N:l:t:j:o:m:M:h" opt; do
		case $opt in
		O ) out=$OPTARG ;;
		N ) name=$OPTARG ;;
		l ) memory=$OPTARG ;;
		t ) time=$OPTARG ;;
		j ) writingOpts=$OPTARG ;;
		o ) outputFile=$OPTARG ;;
		m ) emailOpts=$OPTARG ;;
		M ) email=$OPTARG ;;
		h ) usage; exit 0;
		esac
	done
}

check_files() {
	if [ -z "$cell_vcf" ]; then
		echo "The vcf file does not exist: $cell_vcf"
		usage
		exit 1
	fi
        if [ -z "$gt_vcf" ]; then
                echo "The vcf file does not exist: $gt_vcf"
                usage
                exit 1
        fi
}


generate_pbs() {
	echo "#PBS -A GT-js585
#PBS -N $name
#PBS -l nodes=2:ppn=4
#PBS -l $time
#PBS -j $writingOpts
#PBS -o $outputFile
#PBS -m $emailOpts
#PBS -M $email

cd \$PBS_O_WORKDIR

bgzip $vcf
tabix -p vcf $vcf"".gz

vcftools --gzvcf $vcf"".gz --max-missing 0.""$perct --recode --out filter_$perct
" > demux.pbs
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	qsub demux.pbs
}


main "$@"
