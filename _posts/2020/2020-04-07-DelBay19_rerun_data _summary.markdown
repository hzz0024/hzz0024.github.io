---
comments: true
title: DelBay19 rerun data trimming & mapping summary
date: '2020-04-07 23:00'
tags:
  - rerun
  - mapping
  - trimming
  - WGS
categories:
  - WGS data analysis
---

The DelBay19 rerun samples (pool5) are trimmed and merged (with bam files after bt2 mapping) with previous processed datasets. The overall summary includes the number of reads, reads retained after trimming, genome coverage, which are listed [here](https://docs.google.com/spreadsheets/d/14A5CuNT15jhAgE89HAAXj6ddB3THszFZTpwYOKmu05M/edit?usp=sharing). Some featured values are listed below,

- the whole dataset includes 343 inviduals. Among them, one sample with extreme low mapping rate (SR0719ch_327, 1.7%) and nine individuals with relative low amount of reads (< 2 million) were removed from the final bamlist during angsd run (but I still keep the bamfiles of these samples). 

| Population | pool1 | pool2 | pool3 | pool4 | rerun & merged | sum |
| -----------|-------|-------|-------|-------|----------------|-----|
|     ch1    |   0   |   1   |   1   |   1   |       9        |  12 |
|     ch2    |   0   |   3   |   0   |   0   |       9        |  12 |
|     ch3    |   1   |   1   |   1   |   0   |       9        |  12 |
|     ch4    |   0   |   2   |   3   |   2   |       7        |  14 |
|     ref    |   9   |   6   |   0   |   1   |      32        |  48 |
|     HC     |  14   |   6   |   2   |   2   |      24        |  48 |
|     ARN    |   0   |   3   |   8   |  17   |      19        |  47 |
|     COH    |   3   |   6   |   7   |  12   |      16        |  44 |
|     SR     |   5   |   7   |  11   |  11   |      14        |  48 |
|     NB     |   4   |   8   |  11   |   1   |      24        |  48 |
| -----------|-------|-------|-------|-------|----------------|-----|
|            |       |       |       |       |                | 333 |

- the reads per sample ranges from 1.98 to 40.31 million, with a mean of 12.47 million reads (after trimming)

- the alignment rate ranges from 60.89% to 80.76%, with a mean value of 79.27%

- the genome coverage ranges from 0.15 to 2.87, with a mean value of 0.91   
note: here I used the same coverage definition in Gemma's Atlantic cod paper, which did not take the unmapped regions into account.





