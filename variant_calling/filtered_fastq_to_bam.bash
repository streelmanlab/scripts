#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Create a bam file from the filtered fastq files

#Usage statement
usage () {
        echo "Usage: filtered_fastq_to_bam.bash filtered_fastq_1 filtered_fastq_2 [options]
		-N Name of Job
		-l memnory
		-t hh:mm:ss time needed, job will be killed if exceeded
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
	name=job1
	memory="mem=128gb"
	time="walltime=40:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputFile="job1.out"
	emailOpts="abe"
	email="ggruenhagen3@gatech.edu"
	while getopts "N:l:t:q:j:o:m:M" opt; do
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
	echo "#PBS -N $name
#PBS -l $memory
#PBS -l nodes=2:ppn=4
#PBS -l $time
#PBS -q $cluster
#PBS -j $writingOpts
#PBS -o $outputFile
#PBS -m $emailOpts
#PBS -M $email

cd $PBS_O_WORKDIR
module purge
module load open64/4.5.1
module load bwa/0.7.4     #loads bwa package
module load samtools      #   loads samtools package
bwa mem -M M_zebra_UMD2a.fasta NUM_R1.fastq_filtered NUM_R2.fastq_filtered > NUM_bwa.sam     #creates the SAM
samtools view -bS NUM_bwa.sam > NUM.bam       #Converts to BAM" > filtered_fastq_to_bam.pbs
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	qsub filtered_fastq_to_bam.pbs
}


main "$@"

# OLD
#
#readNum=$1
#
#cat filter.pbs | sed -e "s/NUM/$readNum/g" > /nv/hp10/ggruenhagen3/scratch/toothseq/fastq/scripts/"$readNum"_tmp.pbs
#qsub /nv/hp10/ggruenhagen3/scratch/toothseq/fastq/scripts/"$readNum"_tmp.pbs
