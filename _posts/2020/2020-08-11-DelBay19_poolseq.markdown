---
comments: true
title: DelBay19 poolseq analyses
date: '2020-08-11 12:00'
tags:
  - DelBay19
  - poolseq
  - ouliter
  - null model
  - WGS
categories:
  - WGS data analysis
---

### Reference about poolseq procedure

In [Evolutionary origins of genomic adaptations in an invasive copepod](https://www.nature.com/articles/s41559-020-1201-y), Stern and Lee (2020) used poolseq datasets to describe how biologial invaders, copepod, rapidly adapt to novel habitats. 

Their procedure for pool-seq data analysis:

1) raw read filtering and trimming -- BBDuk (BBTools package)    
2) genome mapping -- BEW-MEM v.0.7.17    
3) supplementary single-end reads mapping -- NextGenMap v.0.5.5    
4) duplicate reads removing -- Picard v. 2.18.27    
5) realign the regions around insertions and deletions -- GATK v3.8    
6) BAM file into mpileup format with Q>20 -- Samtools v.1.3.1    
7) sites within 3bp of an insertion or deletion were removed and the filtered mpileup conversion to sync format -- PoPoolation2    
8) bi-allelic SNP calling -- poolfstat v.1.0 (R), with parameters shown below,     
   MAF 0.05    
   at least four reads for a base call    
   a minimum of 20 and a maximum of 200 total read counts for all populations    


