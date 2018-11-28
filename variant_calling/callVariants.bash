#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Call the variants in the bam files

#Usage statement
usage () {
	echo "Usage: callVariants.bash <reference> <gatk_file_path> -I <bam_file>... [GATK_options] [pbs_options]
	reference: the reference genome
	gatk_file_path: Path to gatk executable
	
	GATK Options:
	-I bam file(s) to call (mandatory)
	-O name of output (default: out.vcf)
	-S minimum phred-scaled confidence threshold at which variants are called (default: 30)
	-L interval list in format \"-L interval_list\"
	
	PBS Options:
	-N Name of Job
	-l memory (default: mem=128gb)
	-t hh:mm:ss time needed, job will be killed if exceeded (default: walltime=72:00:00)
	-q specifies process queue (default: biocluster-6)
	-j controls what gets written to the ouputfile (default: oe)
	-o name of the job's outputfile (default: vcf.out)
	-m controls when email is sent to the submitter (default: abe)
	-M email of submitter (default: ggruenhagen3@gatech.edu)"
}

#Command-line options
get_input() {
	ref=$1
	gatk=$2
	shift
	shift
	
	gatkOut="out.vcf"
	minPhred="30"
	intList=""
	name="variant calling"
	memory="mem=128gb"
	time="walltime=72:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputFile="vcf.out"
	emailOpts="abe"
	email="ggruenhagen3@gatech.edu"
	while getopts "I:O:S:L:N:l:t:q:j:o:m:M:h" opt; do
		case $opt in
		I ) bams+=("$OPTARG");;
		O ) gatkOut=$OPTARG;;
		S ) minPhred=$OPTARG;;
		N ) name=$OPTARG ;;
		L ) intList=$OPTARG;;
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
	echo "finished getopts"
	echo
}

check_files() {
	if [ ${#bams[@]} -eq 0 ]; then
		echo "At least one bam file must be supplied.Exiting the program."
		usage
		exit 1
	fi
	
	for bam in "${bams[@]}"; do
		if [ -z "$bam" ] || [ ! -f "$bam" ]; then
			echo "The bam file $bam does not exist. Exiting the program."
			usage
			exit 1
		fi
	done
	
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
	bamsString=""
	for bam in "${bams[@]}"; do
		bamsString+="-I $bam"
	done
	echo "$bamsString"
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

cd \$PBS_O_WORKDIR
module purge
module load java

java -jar $gatk -R $ref -T HaplotypeCaller $bamsString -stand_call_conf $minPhred $intList -o $gatkOut" > callVariants.pbs
}

main() {
	get_input "$@"
	#check_files
	list_bams
	generate_pbs
	
	#qsub callVariants.pbs
}


main "$@"
