---
comments: true
title: sampling design comparsion
date: '2020-06-23 12:00'
tags:
  - design
  - sample
  - challenge
categories:
  - WGS data analysis
---

In this post I will compare the different sampling strategies and experimental design in Gompert et al. 2014 and Rey et al. 2020


- Rey et al. 2020 Sea bream paper

<img src="https://hzz0024.github.io/images/samples/sea_bream_sample.jpg" alt="img" width="800"/>

Following that, all samples (256 including the larvae and juveniles) were send for RAD-seq. During snp calling, they used Stacks parameters shown below,

1) a minimum read depth (−m) of 5 × per individual per allele   
2) maximum three mismatches (−M)     
3) a minor allele frequency (maf) threshold of 1%     
4) no missing data in at least 90% of the individuals within each of the three samples     
5) Hardy–Weinberg equilibrium within at least one sample using a P-value threshold of 10^-3 in PLINK. Even strong spatially varying selection can be compatible with HWE proportions    
6) a custom script to detect systematic bias in read counts    

These information are helpful for my data filtering.

In the single-generation selection (SGS) test, the authors considered a panmictic common gene pool from which two samples of size N1 and N2 are drawn within the same generation, either at two different times or in two different environments.

- Gompert et al. 2014 stick insect paper

<img src="https://hzz0024.github.io/images/samples/stick_insect_sample.jpg" alt="img" width="800"/>

Description of the experiment. Blue circles represent Ceanothus and orange circles Adenostoma. A) A total of 500 individuals were transplanted into five replicate blocks. Numbers within circles denote number of recaptured individuals (out of 50 released). B) Map of the experiment. 

Field experiment details

1) individuals initially colleced from Far Hill Adenostoma (FHA, N=500), with both striped and green phenotypes       
2) placed in 500-ml containers with 50/container -- 10 populations      
3) randomly assign individuals to 10 bushes (5 Adenostoma and 5 Ceanothus), each individual has a portion of leg removed      
4) recapture after 8-day selection (mostly due to bird predation)     

GBS 

1) DNA from 500 released and 140 recaptured individuals were used for GBS.      
2) Burrow's composite measure of Hardy-Weinberg was calculated        
3) Linkage disequilibrium (LD) were estimated for each pair of SNPs          
4) only SNPs with maf > 5% were used for deltap calculation        

- DelBay19 study

<img src="https://hzz0024.github.io/images/samples/DelBay19.jpg" alt="img" width="800"/>

In general, the DelBay19 experimental design is similar to the stick insect paper, which follows capture-challenge(viability selection)-recapture. However, what we sequenced are not the entire gene pool but only the individuals after selection (recatpure part) and a ref population sourced from wild. This is more similar to the sea bream paper. That is, the selection and reference populations are independent from each other but from a big gene pool. 

Below is another stickleback paper that used similar sampling stratege.

<img src="https://hzz0024.github.io/images/samples/stickback_sample.jpg" alt="img" width="800"/>






   

