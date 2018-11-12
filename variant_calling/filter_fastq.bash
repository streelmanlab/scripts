#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Filter fastq files

#Usage statement
usage () {
        echo "Usage: filter_fastq.bash fastq_1 fastq_2 [options]
		-N Name of Job
		-l hh:mm:ss time needed, job will be killed if exceeded
		-q specifies process queue
		-j controls what gets written to the ouputfile
		-o name of outputfile
		-m controls when email is sent to the submitter
		-M email of submitter"
}

#Command-line options
get_input() {
	fastq1=$1
	fastq2=$2
	name=job_"$fastq_id"
	time="walltime=10:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputFile="job1.out"
	emailOpts="abe"
	email="ggruenhagen3@gatech.edu"
	while getopts "N:l:q:j:o:m:M:h" opt; do
		case $opt in
		N ) name=$OPTARG ;;
		l ) time=$OPTARG ;;
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
	echo "#PBS -N $name
#PBS -l nodes=2:ppn=4
#PBS -l $time
#PBS -q $cluster
#PBS -j $writingOpts
#PBS -o $outputFile
#PBS -m $emailOpts
#PBS -M $email
cd \$PBS_O_WORKDIR
module load ngsqc_toolkit/2.3.3/
IlluQC.pl -pe $fastq1 $fastq2 N A" > filter_fastq.pbs
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	qsub filter_fastq.pbs
}


main "$@"

# OLD
#
#readNum=$1
#
#cat filter.pbs | sed -e "s/NUM/$readNum/g" > /nv/hp10/ggruenhagen3/scratch/toothseq/fastq/scripts/"$readNum"_tmp.pbs
#qsub /nv/hp10/ggruenhagen3/scratch/toothseq/fastq/scripts/"$readNum"_tmp.pbs
