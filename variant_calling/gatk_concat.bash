#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Prepare the reference file

#Usage statement
usage () {
        echo "Usage: gatk_concat.bash <vcf_directory> <gatk_file_path> [options]
		
		vcf_directory: the directory containing all the vcf files to be concatenated
		gatk_file_path: the path to a gatk executable binary
		
		-O name of output file (default: out.vcf)
		-N Name of Job
		-l memory
		-t hh:mm:ss time needed, job will be killed if exceeded (default: walltime=10:00:00)
		-q specifies process queue
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
	
	dir=$1
	gatk=$2
	shift
	shift
	
	out="out.vcf"
	name="gatk_concat"
	memory="mem=1gb"
	time="walltime=10:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputFile="gatk_concat.out"
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
	if [ -z "$gatk" ] || [ ! -f "$gatk" ]; then
                echo "The reference file does not exist. Exiting the program."
                usage
                exit 1
	fi
}

generate_vcfString() {
	vcfString=""
	for file in `ls -v $dir/*.vcf`; do
		if [ -f "$file" ]; then
			vcfString+=" -I $file "
		fi
	done
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
$gatk GatherVcfs $vcfString -O $out
" > gatk_concat.pbs
}

main() {
	get_input "$@"
	check_files
	generate_vcfString
	generate_pbs
	
	qsub gatk_concat.pbs
}


main "$@"
