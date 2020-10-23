---
comments: true
title: Bootstrap parameter search
date: '2020-10-23 12:00'
tags:
  - DelBay19
  - ouliter
  - bootstrap
  - WGS
  - Fisher's exact
  - 
categories:
  - WGS data analysis
---

Previously, I tried to use the p-values AFTER FDR adjustment as the target values for bootstrap analysis and produce some preliminary results. However, this bootstrap approach is supposed to caputre not only the outlier window arond Fisher outliers, but also the windows harbouring loci with consistent large delta_p changes. For this purpose, I test different parameters that may impact the results.

- Threshold for adjusted p-value

First is the different cut threshold for adjusted p-value. Here the p-value is obtained by combining the p-values from each independent test. I performed a FDR adjustment on the p-value list (for 1.9 M SNPs). During bootstrap analysis, a window of SNPs (i.e. 150) were randomly picked up from the whole SNP list. To obtain the null distribution, I count how many SNPs with adjusted p-values less than 0.05 (figure on the left) or 0.1 (figure on the right). After that, I count how many SNPs with target p-values less than 0.05 or 0.1 in the observated window. To obtain the p-value for bootstrap analysis, I measured how many counts in the null distribution that are larger than the observed count and divide that by the totol number of null data (currently it is 10,000).

<img src="https://hzz0024.github.io/images/Fish_boot/cmp1.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/Fish_boot/cmp2.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/Fish_boot/cmp3.jpg" alt="img" width="800"/>

From the figures above we can see that although increasing the cutoff threshold could lead to discovery of some window outside the single-snp outliers, it could also lead to the occurance of window-based outliers in the control group. 

- Threshold for p-value in the bootstrap

