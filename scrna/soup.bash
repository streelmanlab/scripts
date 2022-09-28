#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Use souporcell to demultiplex individuals

#Usage statement
usage () {
        echo "Usage: bash soup.bash <cell_bam> <reference> <gatk> <barcodes> <soup_py> <num_individuals> [options]
		
		cell_bam: bam file of reads from cells/nuclei (filtered to keep only good reads)
		reference: reference genome file
		gatk: gatk binary
		barcodes: barcodes.tsv.gz from cellranger
		soup_py: souporcell.py file
		num_individuals: number of individuals in the pool
		
		-V vcf from cells/nuclei if you have it (will skip the variant calling step)
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
	reference=$2
	gatk=$3
	barcodes=$4
	soup_py=$5
	num_individuals=$6
	shift
	shift
	shift
	shift
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
		echo "The bam file does not exist: $cell_bam"
		usage
		exit 1
	fi
        if [ -z "$reference" ]; then
                echo "The reference file does not exist: $reference"
                usage
                exit 1
        fi
        if [ -z "$gatk" ]; then
                echo "Cannot find gatk: $gatk"
                usage
                exit 1
        fi
        if [ -z "$barcodes" ]; then
                echo "Cannot find the barcodes file: $barcodes"
                usage
                exit 1
        fi
        if [ -z "$soup_py" ]; then
                echo "Cannot find souporcell.py: $soup_py"
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
module load anaconda3
conda activate souporcell

# Call Variants
gunzip -c $barcodes > barcodes.txt
/storage/home/hcoda1/6/ggruenhagen3/p-js585-0/George/rich_project_pb1/bin//gatk-4.1.8.1/gatk --java-options '-Xmx4g' HaplotypeCaller -R $reference --disable-read-filter MappingQualityAvailableReadFilter -I $cell_bam  -stand-call-conf 30 -O sample.vcf

# Make alt and ref matrices
rm -f soup_alt.mtx
rm -f soup_ref.mtx
vartrix --umi --mapq 30 -b $cell_bam -c barcodes.txt --scoring-method coverage --threads 24 --ref-matrix soup_ref.mtx --out-matrix soup_alt.mtx -v sample.vcf --fasta $reference

# Call Cells
echo '\nUsing souporcell to assign cells (George)\n'
rm -f $out""_pred.tsv
python $soup_py --alt_matrix soup_alt.mtx --barcodes barcodes.txt --num_clusters $num_individuals --ref_matrix soup_ref.mtx -o $out""_pred.tsv -t 24

# Call Doublets
echo '\nUsing troublet to assign doublets (George)\n'
rm -f $out""_dbl.tsv
troublet --alts soup_alt.mtx --clusters $out""_pred.tsv -r soup_ref.mtx > $out""_dbl.tsv

" > soup.pbs
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	qsub soup.pbs
}


main "$@"
