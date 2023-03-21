#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Filter the bam file from cellranger to keep on the confidently mapped reads that it uses

#Usage statement
usage () {
        echo "Usage: bash filter_cr_bam.bash <cell_raw_bam> <reference> [options]
		
		cell_raw_bam: bam file of reads from cells/nuclei (unfiltered file called possorted_genome_bam.bam)
		reference: reference genome	
		
		-N Name of Job
		-t hh:mm:ss time needed, job will be killed if exceeded (default walltime: 4:00:00)
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
	time="4:00:00"
	writingOpts="oe"
	outputFile="filter_cr_bam.out"
	emailOpts="BEGIN,END,FAIL"
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
	echo "#!/bin/bash
#SBATCH -A gts-js585
#SBATCH -J $name
#SBATCH -N 2 --ntasks-per-node=4
#SBATCH -t $time
#SBATCH -o $outputFile
#SBATCH --mail-type=$emailOpts
#SBATCH --mail-user=$email

cd \$SLURM_SUBMIT_DIR
source ~/.bashrc
conda activate r4

# Subset good quality reads
echo '\nSubset good quality reads (George)\n'
samtools view -S -b -q 10 -F 3844 $cell_raw_bam > filtered_qc.bam
samtools view filtered_qc.bam > filtered_qc.sam
grep 'xf:i:25' filtered_qc.sam > filtered_qc_cr.sam

# Keep only cells that Cellranger keeps
cp barcodes.txt filter.txt
sed -i -e 's/^/CB:Z:/' filter.txt
cat filtered_qc_cr.sam | LC_ALL=C grep -F -f filter.txt > filtered_qc_cr_cell.sam

# Prep the BAM file for Variant Calling
echo '\nPrep the BAM file for Variant Calling (George)\n'
samtools view -bT $reference filtered_qc_cr_cell.sam > filtered_qc_cr_cell.bam
samtools reheader possorted_genome_bam.bam filtered_qc_cr_cell.bam > filtered_qc_cr_cell_header.bam
samtools sort filtered_qc_cr_cell_header.bam -@ 24 > filtered_qc_cr_cell_header_sort.bam
mv filtered_qc_cr_cell_header_sort.bam filtered_final.bam
samtools index filtered_final.bam

" > filter_cr_bam.sbatch
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	sbatch filter_cr_bam.sbatch
}


main "$@"
