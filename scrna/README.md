# Scripts used for single cell/nuclei RNA-seq analysis
## Introduction:

## Splitting pools of individuals into single indviduals
### With Genotype Information (demuxlet)
Variants from genotyped inviduals will be searched for in reads from cells/nuclei. First, the reads from ```cellranger counts``` need to be filtered out to keep the confidently mapped reads that cellranger actually uses.
```
bash filter_cr_bam.bash <cell_raw_bam> <reference> [options]
```
Then, call variants in the genotyped inviduals following the variant calling pipeline from the variant_calling folder in this repository. TODO (Clean this stuff up): the files we got from the core for our genotyped inviduals ended up creating vcf files with 2 columns per 1 invidual. Here are a few commands to get these 2 columns into 1 column. This is setup only to handle 4 inviduals at a time and isn't in a nice format to handle lots different inputs. The output file from the variant calling pipeline is called new_b1.vcf. The important lines that need to be run no matter what are the last few where I take sites where all the inviduals don't have the same genotype.
```
bcftools view --max-alleles 2 new_b1.vcf > new_b1_bi.vcf
head -n 1714 new_b1_bi.vcf > new_b1_for_demux.vcf
echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t1B11\t1B18\t1B4\t1B5" >> new_b1_for_demux.vcf
grep -v "^#" new_b1_bi.vcf | ./edit_vcf_4.awk >> new_b1_for_demux.vcf
```
Now, it's time to run demuxlet.
```
bash demuxlet.bash <cell_bam> <gt_vcf> [options]
```
The output file called demux_out.best (or if you're using a different basename, the output file called *.best) contains a columned called 'SNG.1ST', which is the most probable individual.

### Without Genotype Information (souporcell)
Souporcell will categorize cells into groups based on the profile of SNPs in their reads. These groups or clusters represent inviduals. However, without genotyping information, we cannot say which group represents which individual. 
```
bash filter_cr_bam.bash <cell_raw_bam> <reference> [options]
soup.bash <cell_bam> <reference> <gatk> <barcodes> <soup_py> <num_individuals> [options]
```
The output file called soup_out_pred.tsv (or if you're using a different basename, the output file called *_pred.tsv) contains the ___. In addition, the output file called soup_out_dbl.tsv contains information on doublets in the data.
