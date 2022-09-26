#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Use souporcell to demultiplex individuals

#Usage statement
usage () {
        echo "Usage: bash soup.bash <cell_bam> <reference> <gatk> <barcodes> <soup_py> [options]
		
		cell_bam: bam file of reads from cells/nuclei (filtered to keep only good reads)
		reference: reference genome file
		gatk: gatk binary
		barcodes: barcodes.tsv.gz from cellranger
		soup_py: souporcell.py file

		-O basename for output souporcell files (default: soup_out)
		-N Name of Job
		-t hh:mm:ss time needed, job will be killed if exceeded (default walltime: 72:00:00)
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
	
	cell_bam=$1	
	gt_vcf=$2
	shift
	shift

	out="soup_out"
	name="soup"
	memory="mem=64gb"
	time="walltime=72:00:00"
	writingOpts="oe"
	outputFile="soup.out"
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
	if [ -z "$cell_bam" ]; then
		echo "The file does not exist: $cell_bam"
		usage
		exit 1
	fi
        if [ -z "$reference" ]; then
                echo "The reference file does not exist: $gt_vcf"
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

# Call Variants
gunzip -c filtered_feature_bc_matrix/barcodes.tsv.gz > barcodes.txt
/storage/home/hcoda1/6/ggruenhagen3/p-js585-0/George/rich_project_pb1/bin//gatk-4.1.8.1/gatk --java-options "-Xmx4g" HaplotypeCaller -R $reference --disable-read-filter MappingQualityAvailableReadFilter -I $cell_bam  -stand-call-conf 30 -O sample.vcf

# Make alt and ref matrices
vartrix --umi --mapq 30 -b $cell_bam -c barcodes.txt --scoring-method coverage --threads 24 --ref-matrix soup_ref.mtx --out-matrix soup_alt.mtx -v sample.vcf --fasta $reference

# Call Cells
echo "\nUsing souporcell to assign cells (Loial)\n"
python $soup_py --alt_matrix soup_alt.mtx --barcodes barcodes.txt --num_clusters 4 --ref_matrix soup_ref.mtx -o $out""_pred.tsv -t 24

# Call Doublets
echo "\nUsing troublet to assign doublets (Loial)\n"
troublet --alts soup_alt.mtx --clusters soup_pred.tsv -r soup_ref.mtx > $out""_dbl.tsv

" > soup.pbs
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	qsub soup.pbs
}


main "$@"
