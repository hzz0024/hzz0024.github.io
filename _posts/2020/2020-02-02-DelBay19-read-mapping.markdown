---
layout: post
title: DelBay19 Reads Mapping
date: '2020-02-02 12:05'
tags:
  - WGS
  - mapping
  - bowtie2

categories:
  - Data analysis
---
After trimming the raw reads using FastQC and Trimmomatic, I started mapping reads of pool1 samples using bowtie2.

Before mapping, creating a directory named bamfiles
```shell
makedir /bamfiles 	#needs to precede bamfile write command
```
Building the bowtie reference index
```shell
bowtie2-build 6565_refCV_10fixchr.fa cv30
```
Preparing the samplelist, this step is a little bit complicated, I decided to use Excel to make a samplelist. Note that I also edited the trimmed fastq filenames start with sample_ID
```
# A subset of the pipeline used in the pipeline described in Therkildsen and Palumbi: "Practical low-coverage genome-wide sequencing of hundreds of individually barcoded samples for population and evolutionary genomics in non-model species" published in Molecular Ecology Resources. 
# SAMPLELIST looks like (Seq_ID, library ID, PU, sample_ID)		

head -n 5 pool1.list
#113514_HVMLMBGXC_DelBay19_1_	DelBay19_1	113514_HVMLMBGXC	HC0419_017
#113515_HVMLMBGXC_DelBay19_1_	DelBay19_1	113515_HVMLMBGXC	HC0419_023
#113516_HVMLMBGXC_DelBay19_1_	DelBay19_1	113516_HVMLMBGXC	NB0419_016
#113517_HVMLMBGXC_DelBay19_1_	DelBay19_1	113517_HVMLMBGXC	NB0419_017
#113518_HVMLMBGXC_DelBay19_1_	DelBay19_1	113518_HVMLMBGXC	NB0419_019
```

```shell
#!/bin/bash

# running at medium gen1

./Bowtie2_mph.sh 6565_refCV_10fixchr.fa cv30 pool1.list /home/hz269/workdir/hz269/DelBay19_raw1_trim/bamfiles >& pool1_bt2.log &
```
In the future I will explore the difference between VIEW -Q 20 and XS:I filter parameters (using filter_bam_MQ.sh), and perhaps start to explore the three alternative way for mapping evaluation 
1. try bowtie2 with single-copy orthologs only
2. try Stampy or other mappers that seem fine-tuned for high polymorphism cases
3. dont worry about it and do the best I can with maf and depth filters, also HWE filter at pop level.

Good to know this [How does bowtie2 assign MAPQ scores?](http://biofinysics.blogspot.com/2014/05/how-does-bowtie2-assign-mapq-scores.html)

---

#### RESULTS
still wait for the results

