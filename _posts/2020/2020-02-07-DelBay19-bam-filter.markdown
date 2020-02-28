---
layout: post
title: DelBay19 Bamfile Trimming
date: '2020-02-07 12:04'
tags:
  - mapping
  - MAPQ
  - bowtie2
categories:
  - WGS data analysis
---

This step will mark duplicate reads in the bam file (nor removed). 

```shell

nohup ./sort_dedup_index.sh >& sort.dedup.log &

pwd

# /home/hz269/workdir/hz269/dedup_bam

#!/bin/bash

#!/bin/bash
## this script filters down .bam files to remove MapQ < 20, effectively removing non-unique
## but not identical filter to removing reads with XS flag
## run from within bamfiles folder ./filter_bam_MQ.sh >& filterbam.log &

for B in *.bam; do 
    OUTBASE=`echo $B | cut -d '_' -f 1-2`  #remove full base name from ".bam"
    echo $OUTBASE       #show which file is being processed
    samtools view -b -q 20 $B > $OUTBASE'.MQ20.bam'
    #samtools sort $OUTBASE'.MQ20.bam' -o $OUTBASE'.sorted.MQ20.bam'
    #samtools index $OUTBASE'.sorted.MQ20.bam'
done
```
---

This wasn't very quick, for a total of 339 samples (I did not exclude "rerun" samples at this step) it took ~ 14 hours to finish 
