---
comments: true
title: CVreseq outlier detection
date: '2020-04-21 12:00'
tags:
  - CVreseq
  - Fst
  - ouliter
  - WGS
categories:
  - WGS data analysis
---

Inspired by Sutherland et al. 2020, I'd like to identify the poential outliers under domestication selection using the thinned vcf data. Four methods will be used for data analysis. They are summarized below

| 			 |Key parameters |   Status     |  
| -----------|----------|--------------|
|   pcadapt  | K        |              |
|   BayeScan | pr_odds      |     done     |
|    RDA     | standard deviations|              |
| Percentile | percentile|     done     |
 

---
### Pcadapt

---
### Bayescan

- Parameter settings:  
-n 5000 Number of outputted iterations   
-nbp 20 Number of pilot runs   
-pilot 5000 Length of pilot runs   
-burn 50000 Burn-in length   
-pr_odds 10 or 100 Prior odds for the neutral model   

The key parameter of Bayescan running is the pr_odds, which indicates how much more likely we think the neutral model is compared to the model with selection. For example a pr_odds = 10 indicates that we think the neutral model is 10 times more likely than the model with selection. With low pr_odds settings, the number of false positives is likely to be large unless we use a very stringent FDR
threshold (see detailed parameter explaination and tutorial [here](http://evomicsorg.wpengine.netdna-cdn.com/wp-content/uploads/2016/01/BayeScan_BayeScEnv_exercises.pdf))

- Results interpreting

In order to decide if a locus is a good candidate for being under the influence of
selection one need to check the q_value in the *_fst.txt* file



---
### RDA

---
### Percentile
1. first estimate the Fst values using VCFtools 
```sh
./vcftools --vcf vcf_file1.vcf --weir-fst-pop individual_list_1.txt --weir-fst-pop individual_list_2.txt
```
see vcftools [manual](https://vcftools.github.io/man_latest.html) and biostars [post](https://www.biostars.org/p/46858/) for detailed parameter explanation. Here the Fst values are estimated using Weir and Cockerham’s (1984) method

---
### Results: number of outlier SNPs and Fst values



---
### Conclusion and questions


