#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Filter fastq files

#Usage statement
usage () {
        echo "Usage: filter_fastq.bash fastq_1 fastq_2 [options]
		-N Name of Job
		-l memory (default: mem=2gb)
		-t hh:mm:ss time needed, job will be killed if exceeded (default: walltime=10:00:00)
		-q specifies process queue
		-j controls what gets written to the ouputfile
		-o name of outputfile (filter_f_JOBID.out)
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
	shift
	shift
	
	name=filter_f_"$fastq1"
	memory="4gb"
	time="10:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputFile="filter_f_\$PBS_JOBID.out"
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
	
}

generate_pbs() {
	echo "#!/bin/bash
#SBATCH -A gts-js585
#SBATCH -J $name
#SBATCH --mem=$memory
#SBATCH -N 2 --ntasks-per-node=4
#SBATCH -t $time
#SBATCH -q inferno
#SBATCH -o $outputFile
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=gwg@gatech.edu

cd \$PBS_O_WORKDIR
module load ngsqc_toolkit/2.3.3/
module load anaconda3
conda activate r4
IlluQC.pl -pe $fastq1 $fastq2 N A" > filter_fastq.sbatch
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	sbatch filter_fastq.sbatch
}


main "$@"

# OLD
#
#readNum=$1
#
#cat filter.pbs | sed -e "s/NUM/$readNum/g" > /nv/hp10/ggruenhagen3/scratch/toothseq/fastq/scripts/"$readNum"_tmp.pbs
#qsub /nv/hp10/ggruenhagen3/scratch/toothseq/fastq/scripts/"$readNum"_tmp.pbs
