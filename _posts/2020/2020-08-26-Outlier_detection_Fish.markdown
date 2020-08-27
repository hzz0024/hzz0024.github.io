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




