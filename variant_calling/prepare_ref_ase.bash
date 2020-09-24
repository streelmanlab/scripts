#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Prepare a new masked reference at SNP sites from the vcf.
#         Realigning to this masked reference will remove allelic bias.

#Usage statement
usage () {
        echo "Usage: ase_count.bash <reference> <vcf_file_path> <gatk_file_path> <ASEr_file_path> [options]
		Prepare a new masked reference at SNP sites from the vcf. In downstream scripts, realigning to this masked reference will remove allelic bias. Note that the resulting reference will be named masked.fasta.
		
		reference: the path to the reference genome
		vcf_file_path: the path to the vcf file
		gatk_file_path: the path to a gatk executable binary
		ASEr_file_path: the path to the ASEr folder (include ASEr in the file path)
		
		-N Name of Job
		-l memory
		-t hh:mm:ss time needed, job will be killed if exceeded (default: walltime=40:00:00)
		-q specifies process queue
		-j controls what gets written to the ouputfile
		-o name of the job's outputfile (default: prepare_ref_ase.out)
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
	
	ref=$1
	vcf=$2
	gatk=$3
	aser=$4
	shift
	shift
	shift
	
	out="output.table"
	name="prepare_ref_ase"
	memory="mem=16gb"
	time="walltime=40:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputFile="prepare_ref_ase.out"
	emailOpts="abe"
	email="ggruenhagen3@gatech.edu"
	while getopts "O:N:l:t:q:j:o:m:M:h" opt; do
		case $opt in
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
	if [ -z "$aser" ]; then
		echo "Invalid ASEr file path. Please include ASEr in the file path. Exiting the program."
		usage
		exit 1
	fi
	if [ -z "$gatk" ] || [ ! -f "$gatk" ]; then
		echo "Invalid gatk file path. Exiting the program."
		usage
		exit 1
	fi
	if [ -z "$ref" ] || [ ! -f "$ref" ]; then
		echo "The reference file does not exist. Exiting the program."
		usage
		exit 1
	fi
	if [ -z "$vcf" ] || [ ! -f "$vcf" ]; then
		echo "Invalid vcf file path. Exiting the program."
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
module load java/1.8.0_25
module load samtools/0.1.19
module load intel/14.0.2
module load perl/5.14.2
module load vcftools/0.1.14.10
module anaconda3
conda activate seurat

# Purpose: Remove allelic bias where the alt allele has at least 1 mismatch
# Keep only SNPs, mask reference at SNPs, 
vcftools --vcf $vcf --remove-indels --recode
$gatk VariantsToTable -R $ref -V out.recode.vcf -O output.table -F CHROM -F POS -F REF -F ALT -F HET -F HOM-REF -F HOM-VAR -F NCALLED -GF GT
awk '{ if(\$5 > 0) print \$0 }' output.table > output.table.het
python $aser/bin/MaskReferenceFromGATKTable.py $ref output.table.het --outfasta masked.fasta
#awk 'BEGIN {OFS='\t'} { if ( $0 !~ /^#/ ) print $1 , $2 , $3 }' out.recode.vcf > out.recode.bed
#perl $aser/perl_scripts/MaskReferencefromBED.pl out.recode.bed $ref masked.fasta

module purge
# Copied from prepare_ref.bash
module load java/1.8.0_25
$gatk CreateSequenceDictionary -R masked.fasta
module load samtools/0.1.19
samtools faidx masked.fasta
module load open64/4.5.1
module load bwa/0.7.4
bwa index -a bwtsw masked.fasta
# End copied section
" > prepare_ref_ase.pbs
}

main() {
	get_input "$@"
	check_files
	generate_pbs
	
	qsub prepare_ref_ase.pbs
}


main "$@"
