#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Call the variants in the bam files

#Usage statement
usage () {
	echo "Usage: callVariants.bash <reference> <gatk_file_path> [-I <bam_file>... or -D <bam_dir>] [GATK_options] [pbs_options]
	reference: the reference genome
	gatk_file_path: Path to gatk executable
	-I bam file(s) to call (must use -I or -D)
	-D directory of bam files (must use -I or -D)
	
	GATK Options:
	-O name of output (default: out.vcf)
	-S minimum phred-scaled confidence threshold at which variants are called (default: 30)
	-L interval list
	-X intervals to exclude
	-P run with intervals and again without excluding intervals
	
	PBS Options:
	-N Name of Job
	-l memory (default: mem=128gb)
	-t hh:mm:ss time needed, job will be killed if exceeded (default: walltime=72:00:00)
	-q specifies process queue (default: biocluster-6)
	-j controls what gets written to the ouputfile (default: oe)
	-o name of the job's outputfile (default: vcf.out)
	-m controls when email is sent to the submitter (default: abe)
	-M email of submitter (default: ggruenhagen3@gatech.edu)
	-c cores and processes (default: nodes=1:ppn=20)"
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
	gatk=$2
	shift
	shift
	
	gatkOut="out.vcf"
	minPhred="30"
	intList=""
	exList=""
	bothList=""
	name="variant calling"
	memory="mem=128gb"
	cores="nodes=1:ppn=20"
	time="walltime=72:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputFile="vcf.out"
	emailOpts="abe"
	email="ggruenhagen3@gatech.edu"
	while getopts "I:D:O:S:L:XL:P:N:l:t:q:j:o:m:M:c:h" opt; do
		case $opt in
		I ) bams+=("$OPTARG");;
		D ) bamDir=$OPTARG;;
		O ) gatkOut=$OPTARG;;
		S ) minPhred=$OPTARG;;
		N ) name=$OPTARG;;
		L ) intList=$OPTARG;;
		X ) exList=$OTPARG;;
		P ) bothList=$OPTARG;;
		l ) memory=$OPTARG;;
		t ) time=$OPTARG;;
		q ) cluster=$OPTARG ;;
		j ) writingOpts=$OPTARG;;
		o ) outputFile=$OPTARG;;
		m ) emailOpts=$OPTARG;;
		M ) email=$OPTARG;;
		c ) cores=$OPTARG;;
		h ) usage; exit 0;
		esac
	done
}

check_files() {
	if [ ${#bams[@]} -eq 0 ] && [ -z "$bamDir" ]; then
		echo "Either at least one bam file or a directory of bam files must be given.Exiting the program."
		usage
		exit 1
	fi
	
	# If the bamDir is empty, then the user must have supplied -I bam files
	if [ -z "$bamDir" ]; then
		for bam in "${bams[@]}"; do
			if [ -z "$bam" ] || [ ! -f "$bam" ]; then
				echo "The bam file $bam does not exist. Exiting the program."
				usage
				exit 1
			fi
		done
	else
		# The user must have supplied a bam directory
		bamsString=""
		for file in $bamDir/*.bam; do
			if [ -f "$file" ]; then
				bamsString+=" -I $file "
			fi
		done
		if [ -z "$bamsString" ]; then
			echo "No bam files in directory: $bamDir. Exiting the program."
			usage 
			exit 1
		fi
	fi
	
	if [ -z "$ref" ] || [ ! -f "$ref" ]; then
		echo "The reference file does not exist. Exiting the program."
		usage
		exit 1
	fi
	if [ -z "$gatk" ] || [ ! -f "$gatk" ]; then
		echo "Not a valid path to an executable gatk file. Exiting the program."
		usage
		exit 1
	fi
}

list_bams() {
	if [ -z "$bamDir" ]; then
		bamsString=""
		for bam in "${bams[@]}"; do
			bamsString+=" -I $bam "
		done
	fi
}

generate_pbs() {
	echo -n "#PBS -N $name
#PBS -l $memory
#PBS -l $cores
#PBS -l $time
#PBS -q $cluster
#PBS -j $writingOpts
#PBS -o $outputFile
#PBS -m $emailOpts
#PBS -M $email

cd \$PBS_O_WORKDIR
module purge
module load java

java -jar $gatk -R $ref -T HaplotypeCaller $bamsString -stand_call_conf $minPhred -o $gatkOut" > callVariants.pbs

if [ ! -z "$intList" ]; then
	echo -n " -L $intList" >> callVariants.pbs
fi
if [ ! -z "$exList" ]; then
	echo -n " -XL $exList" >> callVariants.pbs
fi
if [ ! -z "$bothList" ]; then
	echo -n " -L $bothList
java -jar $gatk -R $ref -T HaplotypeCaller $bamsString -stand_call_conf $minPhred -o $gatkOut -XL $bothList" >> callVariants.pbs
fi
}

main() {
	get_input "$@"
	check_files
	list_bams
	generate_pbs
	
	qsub callVariants.pbs
}


main "$@"
