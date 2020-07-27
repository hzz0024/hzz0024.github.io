---
comments: true
title: DelBay19 single generation selection (SGS) test
date: '2020-07-27 12:00'
tags:
  - DelBay19
  - deltap
  - ouliter
  - null model
  - theta
  - drift
  - WGS
categories:
  - WGS data analysis
---

### Test for single-generation selection

In [Within-Generation Polygenic Selection Shapes Fitness-Related Traits across Environments in Juvenile Sea Bream](https://www.mdpi.com/2073-4425/11/4/398/htm#app1-genes-11-00398), the author estimated the expected distribution of allele frequency changes (deltap) due to finite sample sizes in the absence of selection to generate null expectations.

### Conceptional steps:

1) a panmictic common gene pool of finite size (i.e. the real population)

2) two samples of size N1 (challenge + control, n = 98) and N2 (challenge, n = 50) are drawn within the same generation

3) compare the observed allele frequency difference - ΔP=P2-P1 with the distribution of ΔP expected from random sampling effects (i.e. due to finite sample size effects) and obtain the quantile outliers or p-values.    

### Perform the SGS test for DelBay19 challenge experiment

I used the *sample* function in the R to randomly draw two p (i.e. allele frequency) values based on allele frequency probability distribution in sample N1 and N2. Then I simply calculate the differences between two allele frequency values - i.e. deltap.

This time I incorporate the window-sized theta values into the SGS test. It is designed that each SNP in a window (in the test below I used the 200 bp as window size) will use the same theta value, and estimate the probability distribution of allele frequency in sample N1 and N2. A potential SNP outlier will be identified based on 99.9% quantile of the delta p distribution (two-side, positive or negative).

```R
# SGS test
Rscript --verbose SGS_DelBay19.R 

# check the result

cat p_values_local.txt | wc -l
# the number of potential outliers
4125 
# extract the allele frequency values from N1 and N2
python3 extract.py
# calculate the actual delta p values for 4125 SNPs
Rscript deltaP_act.R -d /Users/ryan/Documents/Ryan_workplace/DelBay19_HG/10_SGS/SGS/local_theat_SGS_results -p CHR_maf0.05_pctind0.7_cv30.mafs.extracted -q CH_maf0.05_pctind0.7_cv30.mafs.extracted -t 4125 -o obs_deltap.output
```
<img src="https://hzz0024.github.io/images/SGS/allele_1.jpeg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/SGS/allele_2.jpeg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/SGS/deltap.jpeg" alt="img" width="800"/>

Remained question: What is appropriate window size? It may be difficult to know without trying and comparing several (e.g. on one chromosome). How does StDev(theta) vary with window size. Maybe we want the window size where variance of theta is greatest?
  

