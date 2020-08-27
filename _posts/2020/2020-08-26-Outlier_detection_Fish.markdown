---
comments: true
title: DelBay19 Fishers' exact test
date: '2020-08-26 12:00'
tags:
  - DelBay19
  - ouliter
  - WGS
  - Fisher's exact
categories:
  - WGS data analysis
---

### Fisher's exact tests (two draws from a common pool approach)

I performed the tests in both challenge vs. refernce and wild contrasts. 

First is the challenge (CH) vs. reference (REF) group, 

P-value histogram for challenge (CH) vs. reference (REF)

<img src="https://hzz0024.github.io/images/outlier/fisher_exact_chr_pvalue_hist.jpeg" alt="img" width="800"/>

An example of wild transect comparison is the HC (Hope Creek) vs New Beds (NB) group.

P-value histogram for HC (Hope Creek) vs New Beds (NB) 

<img src="https://hzz0024.github.io/images/outlier/fisher_exact_wild_pvalue_hist.jpeg" alt="img" width="800"/>

The detailed scripts are located in the DelBay_project/R_scripts/Fisher_exact/Fish_exact_HG.R

--- 

### Results 

|     Method       | No. SNPs | p-value < 0.05 | fdr < 0.2 | fdr < 0.1  | fdr < 0.05 |
|------------------|----------|----------------|------------|-----------|------------|
| CH vs. REF       | 1732036  |   152978       |   3829     |    478    |    96      |
| HC vs. ARN       | 2083615  |   168062       |   1779     |    193    |    54      |
| HC vs. COH       | 2083615  |   165714       |   1779     |    171    |    62      |
| HC vs. SR        | 2083615  |   169239       |   2320     |    271    |    2       |
| HC vs. NB        | 2083615  |   175197       |   2796     |    442    |    61      |
| ARN vs. COH      | 2083615  |   164884       |   1334     |    102    |    14      |
| ARN vs. SR       | 2083615  |   167113       |   1796     |    141    |    20      |
| ARN vs. NB       | 2083615  |   172733       |   2624     |    318    |    35      |
| COH vs. SR       | 2083615  |   169168       |   1710     |    285    |    63      | 
| COH vs. NB       | 2083615  |   172417       |   2656     |    403    |    79      | 
| SR vs. NB        | 2083615  |   174558       |   3238     |    417    |    33      | 

### Shared outliers between challenge-control and wild contrasts

|Group compared    | fdr < 0.2  | fdr < 0.1  | fdr < 0.05 |
|------------------|------------|------------|------------|
|CH_REF - HC_COH   |      5     |      0     |      0     |    
|CH_REF - HC_ARN   |      4     |      0     |      0     | 
|CH_REF - HC_SR    |      9     |      0     |      0     | 
|CH_REF - HC_NB    |      4     |      0     |      0     | 
|CH_REF - ARN_COH  |      3     |      0     |      0     | 
|CH_REF - ARN_SR   |      4     |      0     |      0     | 
|CH_REF - ARN_NB   |      1     |      0     |      0     | 
|CH_REF - COH_SR   |      2     |      0     |      0     | 
|CH_REF - COH_NB   |      4     |      0     |      0     | 
|CH_REF - SR_NB    |      8     |      0     |      0     | 

### Plot the delta_p against start p

- CH vs. REF outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 3829)

<img src="https://hzz0024.github.io/images/Fish/CH_REF_fdr02.jpg" alt="img" width="800"/>

- CH vs. REF randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/CH_REF_fdr02_random_sample.jpg" alt="img" width="800"/>

In the figures above I put some lines to illustrate why delta_p against p0 have shown such patterns.     
1) it makes sense that all the dots should be located above the y= - x line. That is, for a given p0 with negative delta_p, p0 + delta_p >= 0.    
2) it is weird to see that dots are constrained along the y= -2x + 1 and y -2x + 0.1 in the randomly sampling result. Theoretically, there is no constrain for the delta_p distribution. Let us take a look at the wild contrasts then. 

- HC vs. NB outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 2796)

<img src="https://hzz0024.github.io/images/Fish/HC_NB_fdr02.jpg" alt="img" width="800"/>

- HC vs. NB randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/HC_NB_fdr02_random_sample.jpg" alt="img" width="800"/>   

- ARN vs. COH outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 1334)

<img src="https://hzz0024.github.io/images/Fish/ARN_COH_fdr02.jpg" alt="img" width="800"/>

- ARN vs. COH randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/ARN_COH_fdr02_random_sample.jpg" alt="img" width="800"/> 

In the wild transect comparsions, however, there is no such constrains for the data distribution. This issue may come from two reasons. First is the MAF setting in the Angsd run. I set -doMaf as 0.05 for initial SNP calling across the whole populations (which creates a SNP list for downstream individual population running). The retalivaly high MAF setting may bias the allele frequency distribution in the individual population. Second is that I intially created the SNP lists seperately for the CH-REF and wild contrasts. Given that CH-REF has smaller sample size (N=98), I expect that more SNPs in the challenge or reference population are constrained by this MAF setting.

### Relationship between p0 and p1 in potential outliers

- CH vs. REF outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 3829)

<img src="https://hzz0024.github.io/images/Fish/CH_REF_fdr02_p0_p1.jpg" alt="img" width="800"/>

- CH vs. REF randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/CH_REF_fdr02_p0_p1_random_sample.jpg" alt="img" width="800"/>

- HC vs. NB outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 3829)

<img src="https://hzz0024.github.io/images/Fish/HC_NB_fdr02_p0_p1.jpg" alt="img" width="800"/>

- HC vs. NB randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/HC_NB_fdr02_p0_p1_random_sample.jpg" alt="img" width="800"/>

- HC vs. ARN outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 3829)

<img src="https://hzz0024.github.io/images/Fish/HC_ARN_fdr02_p0_p1.jpg" alt="img" width="800"/>

- HC vs. ARN randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/HC_ARN_fdr02_p0_p1_random_sample.jpg" alt="img" width="800"/>

- ARN vs. COH outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 3829)

<img src="https://hzz0024.github.io/images/Fish/ARN_COH_fdr02_p0_p1.jpg" alt="img" width="800"/>

- ARN vs. COH randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/ARN_COH_fdr02_p0_p1_random_sample.jpg" alt="img" width="800"/>

- COH vs. SR outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 3829)

<img src="https://hzz0024.github.io/images/Fish/COH_SR_fdr02_p0_p1.jpg" alt="img" width="800"/>

- COH vs. SR  randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/COH_SR_fdr02_p0_p1_random_sample.jpg" alt="img" width="800"/>

- SR vs. NB outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 3829)

<img src="https://hzz0024.github.io/images/Fish/SR_NB_fdr02_p0_p1.jpg" alt="img" width="800"/>

- SR vs. NB randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/SR_NB_fdr02_p0_p1_random_sample.jpg" alt="img" width="800"/>