#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Use demuxlet to demultiplex individuals

#Usage statement
usage () {
        echo "Usage: bash demux.bash <cell_bam> <gt_vcf> [options]
		
		cell_bam: bam file of reads from cells/nuclei (filtered to keep only good reads)
		gt_vcf:   vcf file originating from genotyped individuals

  		-n no sbatch (ie run the script as a bash script instead of launching sbatch jobs)
		-O basename for output demuxlet files (default: demux_out)
		-N Name of Job
		-t hh:mm:ss time needed, job will be killed if exceeded (default walltime: 1:00:00)
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

 	no_sbatch=""
	out="demux_out"
	name="demux"
	memory="mem=64gb"
	time="walltime=1:00:00"
	writingOpts="oe"
	outputFile="demux.out"
	emailOpts="abe"
	email="ggruenhagen3@gatech.edu"
	while getopts "n:O:N:l:t:j:o:m:M:h" opt; do
		case $opt in
  		n ) no_sbatch=$OPTARG ;;
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
        if [ -z "$gt_vcf" ]; then
                echo "The vcf file does not exist: $gt_vcf"
                usage
                exit 1
        fi
}


generate_pbs() {
	pbs_str="#PBS -A GT-js585
#PBS -N $name
#PBS -l nodes=2:ppn=4
#PBS -l $time
#PBS -j $writingOpts
#PBS -o $outputFile
#PBS -m $emailOpts
#PBS -M $email

cd \$PBS_O_WORKDIR
"

	bash_str="demuxlet --sam $cell_bam --vcf $gt_vcf --out $out --field GT
"

	if [ -z "$no_sbatch" ]; then
		echo $pbs_str $bash_str > demux.pbs
  		qsub demux.pbs
  	else
   		echo $bash_str > tmp.bash
     		bash tmp.bash
       		rm tmp.bash
   	fi
}

main() {
	get_input "$@"
	check_files
	generate_pbs
}


main "$@"
