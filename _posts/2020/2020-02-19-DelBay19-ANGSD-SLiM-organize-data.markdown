---
layout: post
title: Data Organization
date: '2020-02-19 12:00'
tags:
  - ANGSD
  - SLiM
  - data organization
categories:
  - goals
---

I ended up using the minimum reads output per lane (660M reads) as a target number to calculate how many reads we need for each rerun individual. 

The list is posted here [rerun list](https://docs.google.com/spreadsheets/d/10V7vTdNp7oagq4SlPPfOGA-kgmrmh4x6m4olKCdzB6E/edit#gid=1271490464)

It is time to pause the WGS data analysis (as we need to wait for rerun) and 

1) make myself familar with PCA, Fst, and population structure/admixature analysis using ANGSD

[useful link to compute percentiles and plot distributions of quality scores and depths](https://github.com/mfumagalli/ngsTools/blob/master/Scripts/plotQC.R)

2) renew my knowledge about SLiM and build the model for oyster 

3) organize the raw/trimmed data in BioHPC, the space was eaten up during the analysis. I need to compress the fastq file using gzip, which allow me to save 3/4 space. Note, DO NOT compress the bam files as they are already compressed.

```shell

for file in *.fastq; do
gzip $file
done 

```

---

#### RESULTS

The space was reduced a lot by doing gzip for the trimmed fastq files. I also archived the bam files (before and after MQ20 trimming) into the lab Mac. Now my storage space in the BioHPC has been dropped to 528 GB.

For the future running, I may do the initial trimming and mapping in HP. To make sure that the results are consistent, I sent an request to update the Bowtie2 and angsd.   
