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

First is the different cut threshold for adjusted p-value. 

<img src="https://hzz0024.github.io/images/Fish_boot/cmp1.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/Fish_boot/cmp2.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/Fish_boot/cmp3.jpg" alt="img" width="800"/>




