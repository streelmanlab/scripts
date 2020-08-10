# Scripts used in variant calling
## Introduction:
These scripts can be used to auto-generate and call pbs scripts on Georgia Tech's PACE server. The scripts contained in this folder are used to call variants. The input to the start of this pipeline are fastq files and the output is a vcf file.

For scripts that require file paths to a binary, please provide the full path (including the binary in the path).

## Workflow:
```
bash filter_fastq.bash <fastq_1> <fastq_2> [options]
```
Download [GATK](https://software.broadinstitute.org/gatk/download/)
```
bash prepare_ref.bash <reference> <gatk_file_path> [options]

bash filtered_fastq_to_bam.bash <filtered_fastq_1> <filtered_fastq_2> <reference> <outputFile> [options]
```
Download the picard.jar file from the Latest Release of [Picard](https://broadinstitute.github.io/picard/)
```
bash bam_prep.bash <bamFile> <picard_file_path> [options]

bash haploCalling.bash <reference> <gatk_file_path> [-I <bam_file>... or -D <bam_dir>] [GATK_options] [pbs_options]
```
Note that the most convient way to use haploCalling.bash with lots of files is to use an interval list for the chromosomes and another interval list for the unmapped contigs. The interval list must be in this format "contig:start-stop".

Download [tabix](https://sourceforge.net/projects/samtools/files/tabix/), bzip will come with the tabix download. Uncompress and ```make``` the folder. Add the bgzip and tabix binaries to your path environment variable.
 ```
 bash gatk_concat.bash <vcf_directory> <gatk_file_path> [options]
 
 bash filter.bash <vcf_file> [options]
 ```
### ASE:
Download [ASEr](https://github.com/TheFraserLab/ASEr). Follow their instructions to install, the steps are simple.

This follows the same steps as above with an additional first step and last step. The additional first step is to mask the reference genome at sites of heterozygous SNPs. The point of doing this and the realigning everything is to remove allelic imbalance. The alt allele will always have a lower mapping score than the reference allele if the reference genome is not masked. The additional last step is to calculate read counts per allele, the output from this will be used in R.
```
bash prepare_ref_ase.bash <reference> <ASEr_file_path> [options]

bash filtered_fastq_to_bam.bash <filtered_fastq_1> <filtered_fastq_2> <reference> <outputFile> [options]

bash bam_prep.bash <bamFile> <picard_file_path> [options]

bash haploCalling.bash <reference> <gatk_file_path> [-I <bam_file>... or -D <bam_dir>] [GATK_options] [pbs_options]

bash ase_count.bash <reference> <vcf_file_path> <gatk_file_path> [-I <bam_file>... or -D <bam_dir>] [options]
```
