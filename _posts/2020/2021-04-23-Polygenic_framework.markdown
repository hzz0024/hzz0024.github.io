---
comments: true
title: Polygenic framework
date: '2021-04-23 12:00'
tags:
  - DelBay
  - Wild 
  - WGS
  - C-score
  - polygenic
categories:
  - WGS data analysis
--- 

### What is polygenic selection?

In general, I think the concept of polygenic selection of adaptation is synonyms for genetic redundancy, which infers "situation where more than one combination of genetic variants produces the same phenotype (or overall fitness) over an evolutionary timescale". See Láruson et al. [The Importance of Genetic Redundancy in Evolution](https://www.cell.com/trends/ecology-evolution/fulltext/S0169-5347(20)30116-6). 

### Polygenic framework

The reference paper that helps me develop this polygenic framework is Ehrlich et al. 2020. Using the teleost Fundulus heteroclitus as a model species, they have shown that "populations inhabiting distinct environmental niches exhibit subtle, yet significant changes at many positions in the genome within a single generation, and that these changes lead to genetic differentiation among niches. Subtle changes at multiple genes of small effect may allow organisms to adapt to specific niches over just a single generation."

In this paper they did not apply stringent significance during the early stage of outlier detection. Instead, they tried to combine the raw p-values from each of the outlier identification approaches 1) a permutation analysis (using Fst), 2) a simulation approach, and 3) Barnard’s exact test (similar to Fisher’s exact), and do FDR/Bonferroni correction afterward. In fact, only two SNPs within pond 1 remained significant after multiple test collection. After p-value combination, they only identify three outliers (< FDR 10%) in concordant allele frequency tests (i.e. parallel selection), but found a lot more significant SNPs due to subpopulation polygenic selection.

Below is my polygenic framework,

Populations used for this polygenic framework:

Del19 challenge: REF19-CHR19       
Del20 challenge: REF20-CHR20       
Wild contrasts: SR-HC and NB-HC       

1. perform single-generation selection on each of the population contrast, no FDR correction        
2. perform Barnard's exact test, no FDR correction         
3. check the Spearman's rank correlation coefficient         
4. combine p-value from the tests above for each of the SNP using Z-method, retain SNPs with p-value < 0.05, and SNPs with FDR < 0.05              
5. check the outlier overlaps among contrasts
6. calculate the C-score metric to quantify redundancy
7. perform gene annotation for potential candidates
8. perform random forest on significant outliers in each contrast (perhaps SNPs with combined p-value < 0.05)




