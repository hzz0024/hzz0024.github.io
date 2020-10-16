---
comments: true
title: Trans-specific polymorphism in oyster
date: '2020-10-12 12:00'
tags:
  - TSP
  - trans-specific
  - polymorphism
  - C. virginica
  - C. gigas
  - WGS
categories:
  - Resequencing data analysis
---

This is an initial attempt for the discovery of tran-specific polymorphism (TSP) in oyster species. A simple definition of TSP is the positions that are polymorphic in both *C. virginica* and *C. gigas*.

I will start by literature review, trying to find some TSP introduction and technical difficults of TSP discovery.

---
What is the TSP patterns in the context of balancing selection?

*"Balancing selection is expected to produce molecular and phylogenetic footprints not consistent with neutrality (Fijarczyk and Babik 2015). Molecular footprints include: enrichment of old alleles (e.g., trans-species polymorphisms; TSPs), elevated genetic variation (high π), deficit of rare alleles (D > 0),
excess SNPs at medium allele frequencies, reduced divergence around the balanced locus (low FST), as well as the accumulation of non-synonymous variation in the vicinity of balanced polymorphisms, a phenomenon known as sheltered load (Uyenoyama 2005)."

- Ecological load and balancing selection in curcumboreal barnacles by Nunez et al. 2020. MBE [link](https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msaa227/5902836)

---
What aspects should we focus during TSP identification?

*"Because the calling of trans-specific SNPs (tsSNPs) is particularly sensitive to mismapping errors in repetitive sequences, we applied a set of stringent filters, resulting in 74% of the C. rubella reference genome remaining accessible to base calling in both species, with almost half (47%) of the masked sites in the repeat rich pericentromeric regions. "*

- Long-term balancing selection drives evolution of immunity genes in Capsella by Koenig et al. 2019. eLIFE [link](https://elifesciences.org/articles/43606)

---
What can we obtain from TSP analysis?

1) how is TSP occured across genome (e.g. coding regions, introns, promoters, 5'UTRs, 3'UTRs, and intergenic regions).     
2) an enrichment test to show which category is enriched by TSP (e.g. coding regions in the case of barnacle)  
3) perhaps check the exons bearing TSPs using Tajim's D, pi, Fst and so on.      
4) annotate the TSP and check nonsynonymous vs synonymous changes

---
Is there other options for balancing selection footprints detection?

Yes, a recent paper - [Detection of Shared Balancing Selection in the Absence of Trans-Species Polymorphism](https://academic.oup.com/mbe/article/36/1/177/5150441) by Cheng and DeGiorgio proposed a method for balancing selection detection withou TSP. Two related softwares, MULLET (MULti-species LikElihood Tests) and MuteBaSS (MUlTi-spEcies BAlancing Selection Summaries) were developed from this study [link here](http://degiorgiogroup.fau.edu/software.html).

- Detection of Shared Balancing Selection in the Absence of Trans-Species Polymorphism by Xiaoheng Cheng, Michael DeGiorgio. 2019. MBE [link](https://academic.oup.com/mbe/article/36/1/177/5150441)

Another software, [BalLeRMix](https://github.com/bioXiaoheng/BalLeRMix) is developed to include multiple Mixture Model approach for the robust detection of balancing selection. 

- Flexible Mixture Model Approaches That Accommodate Footprint Size Variability for Robust Detection of Balancing Selectio by Xiaoheng Cheng, Michael DeGiorgio. MBE. 2020 [link](https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msaa134/5848012) 





