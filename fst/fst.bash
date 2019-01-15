#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Sort a vcf file

#Usage statement
usage () {
        echo "Usage: bash fst.bash <vcf_file> <pop_1> <pop_2> [options]
		
		vcf_file: vcf file to be sorted
		pop_1: a list of file names of individuals from population 1 separated by a new line
		pop_1: a list of file names of individuals from population 2 separated by a new line
		
		-I Fst window size (default: 10000)
		-T Fst window step (default: 10000)
		-O name of output vcf file (default: fst.vcf)
		-N Name of Job
		-t hh:mm:ss time needed, job will be killed if exceeded (default walltime=30:00:00)
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
	pop1=$2
	pop2=$3
	shift
	shift
	shift
	
	winSize="10000"
	winStep="10000"
	out="fst.vcf"
	name="fst"
	memory="mem=128gb"
	time="walltime=30:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputFile="fst.out"
	emailOpts="abe"
	email="ggruenhagen3@gatech.edu"
	while getopts "I:T:O:N:l:t:q:j:o:m:M:h" opt; do
		case $opt in
		I ) winSize=$OPTARG ;;
		T ) winStep=$OPTARG ;;
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
		echo "The directory does not exist: $dir"
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

vcftools --vcf $vcf --fst-window-size $winSize --fst-window-step $winStep --weir-fst-pop $pop1 --weir-fst-pop $pop2 --out $out
" > fst.pbs
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	qsub fst.pbs
}


main "$@"
