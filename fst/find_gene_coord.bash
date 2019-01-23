#!/usr/bin/env bash

#Author: George Gruenhagen
#Purpose: Find the coordinates of the input genes

#Usage statement
usage () {
	echo "Usage: bash find_gene_coord.bash <treefam_annot> [options]
	
	treefam_annot: treefam annotation
	"
}

#Command-line options
get_input() {
	check_for_help "$1"
	
	tree=$1
	shift
	
	genesStr=""
	while getopts "G:" opt; do
		case $opt in
		G ) genesStr+=$OPTARG"," ;;
		h ) usage; exit 0;
		esac
	done
}

check_for_help() {
	if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
		usage;
		exit 0;
	fi
}

find_coord() {
echo "$genesStr"
IFS=',' read -ra genes <<< "$genesStr"
echo -n "c( "
for gene in "${genes[@]}"; do
	awk -v gene=$gene 'BEGIN {ORS=", "} { if ($2 == gene) { print "c(" $4 ", " $5 ", " "\"" $6 "\")"; exit; } }' $tree
done
echo " )"
}

main() {
	check_for_help
	get_input "$@"
	find_coord
}


main "$@"