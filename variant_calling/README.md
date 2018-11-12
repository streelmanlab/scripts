# Scripts used in variant calling
## Workflow:
```
filter_fastq.bash fastq_1 fastq_2 [options]

prepare_ref.bash reference gatk_file_path [options]
```
Download [GATK](https://software.broadinstitute.org/gatk/download/)
```
filtered_fastq_to_bam.bash filtered_fastq_1 filtered_fastq_2 reference outputFile [options]
```
