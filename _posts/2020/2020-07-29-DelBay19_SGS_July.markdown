---
comments: true
title: DelBay19 single generation selection (SGS) test - corrected
date: '2020-07-29 12:00'
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

1) The posterior probability distribution of 𝑍 in the N1 are generated based on allele counts (k), sample size (n), and theta (sliding with 200 window-size) produced from Angsd. 

2) I use the *sample* function in the R to randomly draw a *p* (i.e. allele frequency) value based on posterior probability distribution of Z in sample N1. This *p* is used to predict the allele count in the sample size n2 (i.e. n2_k = *p* X n2 X 2, n2 = challenge group with 50 samples). 

3) The posterior probability distribution of Z in sample N2 are generated based on n2_k, n2, and the same theta.

4) The expected allele freqeuncy for N1 (p1) and N2 (p2) are drew 1000 times using the *sample* function again. The distribution of allele frequency differences (delta_p) between the two samples is calculated with delta_p = p2 - p1 and recorded for potential outliers.  

5) The expected null distribution is compared to the observed value of Δp with quantile function or P-value estimation (not completed yet). A potential SNP outlier will be identified based on 99.9% quantile of the delta p distribution (two-side, positive or negative).

### Important updates compared to SGS post early this week [https://hzz0024.github.io/2020/07/27/DelBay19_SGS_theta.html](https://hzz0024.github.io/2020/07/27/DelBay19_SGS_theta.html)

1) I discard the *hash* function in model as the same expected delta_p distribution will be repeatedly used for outlier detection. That is, once (n, k, n2_k, n2, theta) are same for two different SNPs, their null-distribution of delta_p will be identical (althought this is less likely to happen). This *hash* is initially designed to save the computational time. Now I just divide the total SNP list into 10 parts and run them parallelly for effeciency.

2) I double check the whole input and output files and found an important error - the observed delta_p file for outlier detection is wrongly calculated. It should be challenge - control but what I used is control - challenge. It would be easy to make mistake like this if I using a another R script (get_deltaP.R) for observed delta_p calculation. To prevent such mistakes in the future, I calculate the observed delta_p directly in my SGS code. Now the inputs for SGS model are just a theta and two maf (challenge and control) files. After these corrections I ran the model again. The results are shown below. 

```R
# SGS test
Rscript --verbose SGS_DelBay19.R 

# check the result

cat p_values_local.txt | wc -l
# the number of potential outliers
1182
# extract the allele frequency values from N1 and N2
python3 extract.py
# calculate the actual delta p values for 1182 SNPs
Rscript deltaP_act.R -d ~/SGS/local_theat_SGS_results -p CHR_maf0.05_pctind0.7_cv30.mafs.extracted -q CH_maf0.05_pctind0.7_cv30.mafs.extracted -t 1182 -o obs_deltap.output
# number of outliers with postive change 
length(dp[dp>0])
> 549
# number of outliers with postive change 
length(dp[dp<0])
> 633
```
<img src="https://hzz0024.github.io/images/SGS/allele_0729.jpeg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/SGS/delta_p_0729.jpeg" alt="img" width="800"/>

Remained question: 

1) What is the repeatability of this model? Is 10000 times of drawing enough for outlier detection? 

So far I have performed this model twice using the same data, only 86 potential outliers overlapped in the two results. It is not unexpected as the posterior probability distribution of Z in sample N2 is largely depended on the k2, whereas k2 is a one-time draw from the N1 distribution. I will perform more test to identify the effect of drawing times.

2) What is appropriate window size? It may be difficult to know without trying and comparing several (e.g. on one chromosome). How does StDev(theta) vary with window size. Maybe we want the window size where variance of theta is greatest?

3) What is the genetic differntiation or diversity indice patterns along these potential outliers? 

4) Are there any common shared outliers between SGS results and wild transect comparsions? How to determine the potential SNPs of selection in the wild transect? - could be association analysis or the same SGS model for outlier identification in wild transect.
  

