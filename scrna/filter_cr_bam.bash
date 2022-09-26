#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Filter the bam file from cellranger to keep on the confidently mapped reads that it uses

#Usage statement
usage () {
        echo "Usage: bash filter_cr_bam.bash <cell_raw_bam> <reference> [options]
		
		cell_raw_bam: bam file of reads from cells/nuclei (unfiltered file called possorted_genome_bam.bam)
		reference: reference genome	
		
		-N Name of Job
		-t hh:mm:ss time needed, job will be killed if exceeded (default walltime: 40:00:00)
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
	
	cell_raw_bam=$1	
	reference=$2
	shift
	shift

	name="filter_cr_bam"
	memory="mem=64gb"
	time="walltime=40:00:00"
	writingOpts="oe"
	outputFile="filter_cr_bam.out"
	emailOpts="abe"
	email="ggruenhagen3@gatech.edu"
	while getopts "N:l:t:j:o:m:M:h" opt; do
		case $opt in
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
	if [ -z "$cell_raw_bam" ]; then
		echo "The vcf file does not exist: $cell_raw_bam"
		usage
		exit 1
	fi
        if [ -z "$reference" ]; then
                echo "The reference file does not exist: $reference"
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

# Subset good quality reads
echo "\nSubset good quality reads (Loial)\n"
samtools view -S -b -q 10 -F 3844 $cell_raw_bam > filtered_qc.bam
samtools view filtered_qc.bam > filtered_qc.sam
grep "xf:i:25" filtered_qc.sam > filtered_qc_cr.sam

# Prep the BAM file for Variant Calling
echo "\nPrep the BAM file for Variant Calling (Loial)\n"
samtools view -bT $reference filtered_qc_cr.sam > filtered_qc_cr.bam
samtools reheader possorted_genome_bam.bam filtered_qc_cr.bam > filtered_qc_cr_header.bam
samtools sort filtered_qc_cr_header.bam > filtered_qc_cr_header_sort.bam
mv filtered_qc_cr_header_sort.bam filtered_final.bam
samtools index filtered_final.bam

" > filter_cr_bam.pbs
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	qsub filter_cr_bam.pbs
}


main "$@"
