#!/bin/bash
#SBATCH -A gts-js585
#SBATCH -J /storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/dna_seq/cc/cc.bam
#SBATCH --mem=128gb
#SBATCH -N 2 --ntasks-per-node=4
#SBATCH -t 40:00:00
#SBATCH -q inferno
#SBATCH -o /storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/dna_seq/cc/cc.bam.out
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ggruenhagen3@gatech.edu

cd $PBS_O_WORKDIR
module purge
module load anaconda3
module load bwa/0.7.17
conda activate r4
bwa mem -M /storage/home/hcoda1/6/ggruenhagen3/scratch/m_zebra_ref/GCF_000238955.4_M_zebra_UMD2a_genomic.fna /storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/dna_seq/cc/filtered/CC_filtered_R1_001.fastq.gz /storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/dna_seq/cc/filtered/CC_filtered_R2_001.fastq.gz > /storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/dna_seq/cc/cc.bam.sam     #creates the SAM
samtools view -bS /storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/dna_seq/cc/cc.bam.sam > /storage/home/hcoda1/6/ggruenhagen3/scratch/d_tooth/data/dna_seq/cc/cc.bam       #Converts to BAM
