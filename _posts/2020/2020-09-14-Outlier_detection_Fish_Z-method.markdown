---
comments: true
title: DelBay19 outlier Fisher’s approach vs weighted Z-method
date: '2020-09-14 12:00'
tags:
  - DelBay19
  - ouliter
  - WGS
  - Fisher's exact
  - 
categories:
  - WGS data analysis
---

In this post I made some modifications in Fisher's exact tests. I also compared the results from Fisher's approach vs weighted Z-method.  

1) Modification on Fisher's approach. 

Previously I performed Fisher's approach on DelBay datasets. The replicates are set as challenge vs reference plus wild contrast: either Hope Creek (HC) vs Shell Rock (SR) or Hope Creek (HC) vs New Bed (NB). I used REF vs SR + COH vs ARN as the control group.  

However, I found an important error when I double check the codes in [Dixon et al. (2015) Genomic determinants of coral heat tolerance across latitudes](https://science.sciencemag.org/content/348/6242/1460). This happens when they tried to combine the p-values from replicate tests:

```sh
#FUNCTION fisher.method():
#performs fisher's method to combine a set of p values
#ARGUMENTS: ps = a vector of p values you want combined
fisher.method=function(ps) {
  ftest=-2*sum(log(ps))
  df=length(2*ps)
  pv=1-pchisq(q=ftest,df=df)
} 
```

Here the *df* sould be *df=2 x length(ps)* but not *df=length(2 x ps)*. The detailed explanation can be found in this post: [https://brainder.org/2012/05/11/the-logic-of-the-fisher-method-to-combine-p-values/]. I am not sure how this might impact the results in coral heat selection results, but this does impact the number of outliers in my datasets.

2) According to [Whitlock (2005)](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1420-9101.2005.00917.x), the weighted Z-test (also called ‘Stouffer’s method) is more power and more precision than does Fisher’s test. It favours symmetric rejection and is less sensitive to a single low p-value, requiring more consistently low p-values to yield a low combined p-value.

Here I used the *combinePValues* function in the [scran: Methods for Single-Cell RNA-Seq Data Analysis](https://rdrr.io/bioc/scran/man/combinePValues.html) package to obtain the combined p-values.

---

### Fisher's exact tests

The R scripts of Fisher's exact test for each repeated study are located in GitHub/DelBay_project/R_scripts/Fisher_exact/Fish_repeat

The R scripts of Fisher’s approach and Z-method are located in GitHub/DelBay_project/R_scripts/Fisher_exact/Fish_combine_p

---

### SNPs and Allele Frequency Files (mafs.gz)

Again I used the regenerated dataset (doMajorMinor 3 + doMaf 1) that includes 1,934,038 SNPs from global Angsd calling. This SNP list is used for indivudal population SNP dataset creation.

---

### Fisher’s approach

|   Test       | fdr < 0.1 | fdr < 0.05| fdr < 0.01|
|--------------|-----------|-----------|-----------|
|REF-CH-SR-HC  |     41    | 10        |      1    |
|REF-CH-NB-HC  |     32    | 16        |      6    |
|SR-REF-COH-ARN|     20    | 0         |      0    | 

|Group compared           | fdr < 0.1  | 
|-------------------------|-------------|
|CH vs. REF & HC vs. SR   |      41     |
|CH vs. REF & HC vs. NB   |      32     |
|Shared                   |      10     |

Manhattan plot for CH vs. REF + HC vs. SR. Red dashed line indicates 10% FDR threshold (No. of outliers = 41)

<img src="https://hzz0024.github.io/images/Fish/manhattan_REF_CH_SR_HC_fish.jpeg" alt="img" width="800"/>

Manhattan plot for CH vs. REF + HC vs. NB. Red dashed line indicates 10% FDR threshold (No. of outliers = 32)

<img src="https://hzz0024.github.io/images/Fish/manhattan_REF_CH_NB_HC_fish.jpeg" alt="img" width="800"/>

Manhattan plot for CH vs. REF + HC vs. NB. Red dashed line indicates 10% FDR threshold (No. of outliers = 20) 

<img src="https://hzz0024.github.io/images/Fish/manhattan_SR-REF-COH-ARN_fish.jpeg" alt="img" width="800"/>

---

### Fisher’s approach

|   Test       | fdr < 0.1 | fdr < 0.05| fdr < 0.01|
|--------------|-----------|-----------|-----------|
|REF-CH-SR-HC  |     41    | 10        |      1    |
|REF-CH-NB-HC  |     32    | 16        |      6    |
|SR-REF-COH-ARN|     20    | 0         |      0    | 

|Group compared           | fdr < 0.1  | 
|-------------------------|-------------|
|CH vs. REF & HC vs. SR   |      41     |
|CH vs. REF & HC vs. NB   |      32     |
|Shared                   |      10     |

Manhattan plot for CH vs. REF + HC vs. SR. Red dashed line indicates 10% FDR threshold (No. of outliers = 41)

<img src="https://hzz0024.github.io/images/Fish/manhattan_REF_CH_SR_HC_fish.jpeg" alt="img" width="800"/>

Manhattan plot for CH vs. REF + HC vs. NB. Red dashed line indicates 10% FDR threshold (No. of outliers = 32)

<img src="https://hzz0024.github.io/images/Fish/manhattan_REF_CH_NB_HC_fish.jpeg" alt="img" width="800"/>

Manhattan plot for SR vs. REF + COH vs. ARN. Red dashed line indicates 10% FDR threshold (No. of outliers = 20) 

<img src="https://hzz0024.github.io/images/Fish/manhattan_REF_CH_NB_HC_fish.jpeg" alt="img" width="800"/>

---

### Z or Stouffer’s approach

|   Test       | fdr < 0.1 | fdr < 0.05| fdr < 0.01|
|--------------|-----------|-----------|-----------|
|REF-CH-SR-HC  |     11    |  4        |      0    |
|REF-CH-NB-HC  |     31    |  8        |      1    |
|SR-REF-COH-ARN|      0    |  0        |      0    | 

|Group compared           | fdr < 0.1  | 
|-------------------------|-------------|
|CH vs. REF & HC vs. SR   |      11     |
|CH vs. REF & HC vs. NB   |      31     |
|Shared                   |       2     |

Manhattan plot for CH vs. REF + HC vs. SR. Red dashed line indicates 10% FDR threshold (No. of outliers = 11)

<img src="https://hzz0024.github.io/images/Fish/manhattan_REF_CH_SR_HC_z.jpeg" alt="img" width="800"/>

Manhattan plot for CH vs. REF + HC vs. NB. Red dashed line indicates 10% FDR threshold (No. of outliers = 31)

<img src="https://hzz0024.github.io/images/Fish/manhattan_REF_CH_NB_HC_z.jpeg" alt="img" width="800"/>

Manhattan plot for SR vs. REF + COH vs. ARN. Red dashed line indicates 10% FDR threshold (No. of outliers = 0) 

<img src="https://hzz0024.github.io/images/Fish/manhattan_SR-REF-COH-ARN_z.jpeg" alt="img" width="800"/>









