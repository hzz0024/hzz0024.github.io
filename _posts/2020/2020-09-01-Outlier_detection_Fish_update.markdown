---
comments: true
title: DelBay19 Fishers' exact test update
date: '2020-09-01 12:00'
tags:
  - DelBay19
  - ouliter
  - WGS
  - Fisher's exact
categories:
  - WGS data analysis
---

### Fisher's exact tests

In this post I made some modifications for Fisher's exact tests. 

1) The SNP lists in the last post [DelBay19 Fishers' exact test](https://github.com/hzz0024/hzz0024.github.io/blob/master/_posts/2020/2020-08-26-Outlier_detection_Fish.markdown) were created seperatly for the challenge and wild groups. I suspect this will reduce the number of common shared outliers in the test results. Here I create a single SNP list by combining all the populations in the Angsd global run, and generated the allele frequency dataset for each population.

2) I am still exploring the best way to obtain the allele frequency values, should be either doMajorMinor 3 plus doMaf 1 or doMajorMinor 5 plus doMaf 2. I made a comparsion in this post.

3) The most important change is the Fisher's test. Following the coral paper publisted by [matzlab](https://matzlab.weebly.com/), I performed Fisher’s exact tests on two "replicate" datasets. One is the challenge vs reference contrast, another one is the wild contrast from Cohansey (COH) vs Shell Rock (SR). P-values from the replicates were then combined using *Fisher’s Combined Probability method*. This approach requires the use of one-tailed tests so that information on the direction of effect is retained, with alternative hypothesis derived from the average of delta_p. More importantly, the p-value is by 2 to compensate for this post hoc use of the data to inform one-tailed tests. 

Whitlock made a good explaination for this method, which can be found in this paper [Combining probability from independent tests: The weighted Z-method is superior to Fisher’s approach](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1420-9101.2005.00917.x). This method, however, will treats large and small p-values asymmetrically. That is, the combined p-value is sensitive to small p-values compared to large p-values, and is likely to reject the null hypothesis in favour of contradictory one-tailed alternative hypothesis. Anyway, given that the code is already developed for this method, I modified the code to fit into my data and tried to identify the poential outliers from the "replicates". 



































