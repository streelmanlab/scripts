#!/bin/bash

#Author: George Gruenhagen
#Purpose: Call the variants in the bam files

#Usage statement
usage () {
	echo "Usage: haploCalling.bash <reference> <gatk_file_path> [-I <bam_file>... or -D <bam_dir>] [GATK_options] [pbs_options]
	reference: the reference genome
	gatk_file_path: Path to gatk executable
	-I bam file(s) to call (must use -I or -D)
	-D directory of bam files (must use -I or -D)
	
	GATK Options:
	-O name of output directory (default: vcfDir)
	-S minimum phred-scaled confidence threshold at which variants are called (default: 30)
	-L interval list
	-X intervals to exclude
	-P run with intervals and again without excluding intervals
	-B batch size: how many jobs run together (default: 32)
	
	PBS Options:
	-N Name of Job
	-l memory (default: mem=128gb)
	-t hh:mm:ss time needed, job will be killed if exceeded (default: walltime=60:00:00)
	-q specifies process queue (default: biocluster-6)
	-j controls what gets written to the ouputfile (default: oe)
	-o name of the job's output directory where log files are stored (default: logs)
	-m controls when email is sent to the submitter (default: abe)
	-M email of submitter (default: ggruenhagen3@gatech.edu)"
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
	
	gatkOut="vcfDir"
	minPhred="30"
	batch_size=32
	intList=""
	exList=""
	bothList=""
	name="variant calling"
	memory="mem=128gb"
	cores="nodes=1:ppn=1"
	time="walltime=60:00:00"
	cluster="biocluster-6"
	writingOpts="oe"
	outputDir="logs"
	emailOpts="abe"
	email="ggruenhagen3@gatech.edu"
	while getopts "I:D:O:S:N:L:X:P:B:l:t:q:j:o:m:M:c:h" opt; do
		case $opt in
		I ) bams+=("$OPTARG");;
		D ) bamDir=$OPTARG;;
		O ) gatkOut=$OPTARG;;
		S ) minPhred=$OPTARG;;
		N ) name=$OPTARG;;
		L ) intList=$OPTARG;;
		X ) exList=$OPTARG;;
		P ) bothList=$OPTARG;;
		B ) batch_size=$OPTARG;;
		l ) memory=$OPTARG;;
		t ) time=$OPTARG;;
		q ) cluster=$OPTARG ;;
		j ) writingOpts=$OPTARG;;
		o ) outputDir=$OPTARG;;
		m ) emailOpts=$OPTARG;;
		M ) email=$OPTARG;;
		c ) cores=$OPTARG;;
		h ) usage; exit 0;
		esac
	done
}

check_files() {
	mkdir -p "$gatkOut"
	
	mkdir -p "$outputDir"
	
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

generate_jobs() {
	echo -n "" > jobs.txt
	k=0
	while read line; do
		echo "$gatk --java-options \"-Xmx4g\" HaplotypeCaller -R $ref  $bamsString -stand-call-conf $minPhred -O $gatkOut/out$k"".vcf -L \"$line\"" >> jobs.txt
		k=$((k+1))
	done < $list
}

generate_xl_jobs() {
	# One job for everything not on a chromosome
	if [ ! -z "$bothList" ]; then
		# Convert from chr1:1-100 format to bed format (separated by tabs)
		awk 'BEGIN {OFS="\t"} { p=index($0,":"); pos=substr($0,p+1);chrom=substr($0, 0,p-1); dash=index(pos,"-"); start=substr(pos,0,dash-1); stop=substr(pos,dash+1); print chrom, start, stop }' $bothList > pos.bed
		
		# Append the unplaced contig job at the end of all the chromosome jobs
		echo "$gatk --java-options \"-Xmx4g\" HaplotypeCaller -R $ref  $bamsString -stand-call-conf $minPhred -O $gatkOut/xl"".vcf -XL pos.bed" >> jobs.txt
	else
		# Convert from chr1:1-100 format to bed format (separated by tabs)
		awk 'BEGIN {OFS="\t"} { p=index($0,":"); pos=substr($0,p+1);chrom=substr($0, 0,p-1); dash=index(pos,"-"); start=substr(pos,0,dash-1); stop=substr(pos,dash+1); print chrom, start, stop }' $exList > pos.bed
		
		# Unplaced contigs are the only jobs, make a new file
		echo "$gatk --java-options \"-Xmx4g\" HaplotypeCaller -R $ref  $bamsString -stand-call-conf $minPhred -O $gatkOut/xl"".vcf -XL pos.bed" > jobs.txt
	fi
}

generate_multi_pbs() {
	echo $time
	echo '#PBS -N multi-paralleljob
#PBS -q biocluster-6 
#PBS -l '"$time"'
#PBS -l nodes=4:ppn=8
#PBS -j oe
#PBS -o '"$outputDir"'/out.$PBS_JOBID


cd $PBS_O_WORKDIR
NP=$(wc -l < $PBS_NODEFILE)
ls ../bin/gatk-4.0.11.0

# add all modules needed here
module load java/1.8.0_25
module load gnuparallel/20150422

#JOBFILE, BATCHSIZE, and BATCHNUM should be set in the environment
#If they are not, use some defaults.
# By setting BATCHSIZE to a default of the length of the jobfile we only require one of these jobs.
# The user can submit multiple jobs and split up the batchcount to work on multiple nodes.
JOBFILE=${JOBFILE:-jobs.txt}

if [ ! -f $JOBFILE ]; then echo "File $JOBFILE does not exist. Exiting"; exit 0; fi

BATCHSIZE=${BATCHSIZE:-$(wc -l < $JOBFILE)}
BATCHNUM=${BATCHNUM:-0}

JOBCOUNT=$(wc -l < $JOBFILE)

ENDLINE=$(($BATCHSIZE*$BATCHNUM + $BATCHSIZE))

if [ $ENDLINE -gt $JOBCOUNT ]
then

  if [ $(($ENDLINE-$BATCHSIZE)) -gt $JOBCOUNT ]
  then
    echo "Given \"BATCHNUM\" is greater than the number of possible batches. Exiting..."
    exit 0
  fi

  DIFFERENCE=$(($ENDLINE-$JOBCOUNT))
  REMAININGJOBCOUNT=$(($BATCHSIZE-$DIFFERENCE))

fi

BATCHSIZE=${REMAININGJOBCOUNT:-$BATCHSIZE}

head -n $ENDLINE $JOBFILE | tail -n $BATCHSIZE | /usr/local/pacerepov1/gnuparallel/20110822/bin/parallel -j $NP -k 
' > paralleljob.txt
}

make_batch() {
	for ((i=0;i<$jobs;i+=$batch_size)); do
		qsub -vBATCHSIZE=$batch_size,BATCHNUM=$((i / batch_size)) paralleljob.txt
	done
	
	rem=$(( jobs % batch_size ))
	
	if [ $rem != 0 ]; then
		# This is just so that way each batch has as many jobs as nodes*ppn
		for ((j=0;j<$(( batch_size - rem));j+=1)); do
			echo "ls $gatk" >> jobs.txt
		done
	fi
}

main() {
	get_input "$@"
	check_files
	list_bams
	
	if [ ! -z "$bothList" ] || [ ! -z "$intList" ] || [ ! -z "$exList" ]; then
		if [ ! -z "$bothList" ]; then
			echo "Generating jobs that include the provided regions and jobs that exclude the provided regions"
			list=$bothList
			generate_jobs
			generate_xl_jobs
		elif [ ! -z "$exList" ]; then
			echo "Generating jobs exclude the provided region list"
			generate_xl_jobs
		else
			echo "Generating jobs include the provided region list"
			list=$intList
			generate_jobs
		fi
		generate_multi_pbs
		
		jobs=$(wc -l < jobs.txt)
		make_batch
	else
		echo "Please provide an interval list to call variants on."
		usage
		exit 1
	fi
}


main "$@"
