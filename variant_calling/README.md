# Scripts used in variant calling
## Workflow:
```
bash filter_fastq.bash fastq_1 fastq_2 [options]
```
Download [GATK](https://software.broadinstitute.org/gatk/download/)
```
bash prepare_ref.bash reference gatk_file_path [options]

bash filtered_fastq_to_bam.bash filtered_fastq_1 filtered_fastq_2 reference outputFile [options]
```
Download the picard.jar file from the Latest Release of [Picard](https://broadinstitute.github.io/picard/)
```
bash bam_prep.bash bamFile picard_file_path [options]

bash haploCalling.bash <reference> <gatk_file_path> [-I <bam_file>... or -D <bam_dir>] [GATK_options] [pbs_options]
```
Note that the most convient way to use haploCalling.bash with lots of files is to use an interval list for the chromosomes and another interval list for the unmapped contigs. The interval list must be in this format "contig:start-stop".

Download [tabix](https://sourceforge.net/projects/samtools/files/tabix/), bzip will come with the tabix download. Uncompress and make the folder.
