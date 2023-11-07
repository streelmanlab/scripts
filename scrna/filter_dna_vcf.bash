#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Use demuxlet to demultiplex individuals

#Usage statement
usage () {
        echo "Usage: bash demux.bash <dna_vcf> <num_ind> [options]
		
		dna_vcf: vcf file of genotyped individuals from DNA sequencing
		num_ind: Number of individuals (numeric)"
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
	
	dna_vcf=$1	
	num_ind=$2
	shift
	shift

	while getopts "h" opt; do
		case $opt in
		h ) usage; exit 0;
		esac
	done
}

check_files() {
	if [ -z "$dna_vcf" ]; then
		echo "The vcf file does not exist: $dna_vcf"
		usage
		exit 1
	fi
}


make_edit_vcf_awk() {
  echo "#! /usr/bin/awk -f

{OFS="\t"}
{ 
" > edit_vcf.awk
  for n in $(seq 1 $num_ind); do
      "new_str$n='./.'" >> edit_vcf.awk
  done
  echo "}" >> edit_vcf.awk
}

main() {
	get_input "$@"
	check_files
  	make_edit_vcf_awk
}


main "$@"
