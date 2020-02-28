---
layout: post
title: DelBay19 read coverage estimate
date: '2020-02-11 12:00'
tags:
  - coverage
  - test
  - samtools
  - depth
  - WGS
categories:
  - WGS data analysis
---

I am looking for the answer to "what should be our minimum read count needed from each sample?", this could be done by estimating the sample coverage using filtered bam files (-MQ20), and determining the relative proportions of reads needed for rerun.

I followed the post from here [how to How do you calculate the genome coverage for a sequenced genome?](http://seqanswers.com/forums/showthread.php?t=17725)

```shell
#!/bin/bash

#loop over the trimmed bam files
for sample in *.sorted.cv30.bam; do
cov_depth=$(samtools depth $sample | awk '{sum+=$3;cnt++}END{print sum/cnt" "sum}')
echo $sample: $cov_depth
done

```
---

#### RESULTS

see google excel link here [coverage report](https://docs.google.com/spreadsheets/d/10V7vTdNp7oagq4SlPPfOGA-kgmrmh4x6m4olKCdzB6E/edit#gid=1728449447)

