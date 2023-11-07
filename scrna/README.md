# Scripts used for single cell/nuclei RNA-seq analysis
## Introduction:

## Choosing the optimal clustering paramaters in Seurat (ChooseR)
[ChooseR](https://github.com/rbpatt2019/chooseR) performs clustering at multiple resolutions and finds the optimal resolution that results in the highest number of clusters while the median silhoutte score of bootstraps is still high. I have modified the code to follow this same logic for multiple clustering parameters, not just resolution. The additional parameters that my code will evaulate are min.dist in Seurat's RunUMAP and k.param (n.neighbors) in Seurat's FindNeighbors. Additionally, I modified the chooseR code to follow the clustering workflow we use in the Streelmanlab.

The first step is to use renv as described in ChooseR's [github](https://github.com/rbpatt2019/chooseR). I believe this creates a local R environment in the working directory. Next, modify the line in hb_chooser.R that loads a Seurat object, to load the correct object that you want. Also modify the results path to be what you want. Next call ```bash hb_params.bash```, which will find all the combos of min.dist and n.neighbors, then it will create job scripts for each combo that calls hb_chooseR.R. Then, hb_chooseR.R will do clustering with bootstraps on multiple resolutions using the values of min.dist and n.neighbors it was given. It will save clustering results to intermediate rds files in some directories. 

When all jobs are complete, run ```hb_chooser_collect.R``` to read those results from the rds files and see which clustering parameters are optimal (Note that you will probably have to modify the paths in hb_chooser_collect.R). ChooseR chooses the clustering paramaters by first finding the bottom 95th percentile of the bootstraps for all parameters, then it finds the highest of all of those, then it finds the clustering parameters with the most number of clusters where the median is above that highest bottom 95th percentile.

## Splitting pools of individuals into single indviduals
### With Genotype Information (demuxlet)
**Goal**: in pools that contain multiple indviduals, idenitify which individual a cell came from using demuxlet

**Prerequisites**: complete the cellranger pipeline (ie cellranger mkfastq and count) on a pool

**-- Note to users not on HPC --**: use the ```-n``` flags for ```filter_cr_bam_sbatch.bash``` command (this will execute the commands using a bash script instead launching jobs with a job manager on an HPC)

Demuxlet compares the variants present in the RNA reads of single cells to the DNA reads of genotyped individuals in order to determine which cells came from which individuals. First, the reads from ```cellranger count``` need to be filtered to keep only good quality reads. We will do this with the command below, but first for the sake of this tutorial, let's briefly use some arbitrary names like a real-world example. In this tutorial, let's say the output folder created from ```cellranger count``` is called ```sample1``` and let's say that the reference genome we created from ```cellranger mkref``` is called ```my_ref```. Then, to filter our reads we would execute the command below.
```
bash filter_cr_bam_sbatch.bash sample1/outs/possorted_genome_bam.bam my_ref [options] [pbs_options]
```
The command above just created a filtered bam file called ```filtered_final.bam``` that we can use as input to demuxlet. Before we do that though, variant calling needs to performed on the genotyped individual (see the [variant calling scripts in this repository](../variant_calling/README.md)). Once the variant calling is complete, we will have a vcf file containing a column for every individual. For the sake of this tutorial, let's call that file ```dna.vcf```.

This next step is optional, but generally recommended -> Lets filter the vcf files to include only biallelic SNPs. 
```
vcftools dna.vcf --remove-indels --min-alleles 2 --max-alleles 2 --recode --stdout > dna_bi.vcf
```
Now, let's filter the vcf file from the genotyped individual to only contain SNPs where not all the individuals have the same genotypes.
```
bash filter_dna_vcf.bash dna_bi.vcf <num_ind>
```
Finally, let's run demuxlet!
```
bash demux.bash sample1/outs/filtered_final.bam dna_bi.vcf [options]
```
The output file called demux_out.best (or if you're using a different basename, the output file called *.best) contains a column called 'SNG.1ST', which is the most probable individual.

### Without Genotype Information (souporcell)
Souporcell will categorize cells into groups based on the profile of SNPs in their reads. These groups or clusters represent inviduals. However, without genotyping information, we cannot say which group represents which individual. The ```<cell_raw_bam>``` file is the ```possorted_genome_bam.bam``` file in the cellranger output folder. Before this script can ran, the ```filtered_feature_bc_matrix/barcodes.tsv.gz``` needs to be gunzipped to ```barcodes.txt``` in the cellranger output folder.
```
bash filter_cr_bam.bash <cell_raw_bam> <reference> [options]
soup.bash <cell_bam> <reference> <gatk> <barcodes> <soup_py> <num_individuals> [options]
```
The output file called soup_out_pred.tsv (or if you're using a different basename, the output file called *_pred.tsv) contains the ___. In addition, the output file called soup_out_dbl.tsv contains information on doublets in the data.
