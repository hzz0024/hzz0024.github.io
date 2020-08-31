---
comments: true
title: MAF evaluation
date: '2020-08-30 12:00'
tags:
  - DelBay19
  - MAF
  - WGS
categories:
  - WGS data analysis
---

In this post I collected some Angsd parameters from the related literatures. The main purpose of doing it to figure out why allele frequency values in the mafs file do not go beyond 0.8. Is it due to some technical problems? Or just a refection to the data itself?

--- 

### Angsd parameters in the literatures

- Stahlke et al. 2020. Historical and contemporary signatures of selection in response to transmissible cancer in the Tasmanian Devil (Sarcophilus harrisii)

<img src="https://hzz0024.github.io/images/MAF/TableS2.jpg" alt="img" width="800"/>

- Castruita et al. 2020. Analyses of key genes involved in Arctic adaptation in polar bears suggest selection on both standing variation and de novo mutations played an important role [link](https://link.springer.com/article/10.1186/s12864-020-06940-0)

"We used a genotype likelihood approach to construct the PCAs: input genotype likelihood files were constructed using ANGSD v0.929, with the SAMtools genotype likelihood algorithm (−GL 1), and specifying the following parameters: remove reads that have multiple mapping best hits (−unique_only), remove reads with a flag above 255/secondary hits (−remove_bads), include only read pairs with both mates mapping correctly (−only_proper_pairs), adjust mapQ for reads with excessive mismatches (−C 50), adjust quality scores around indels (−baq 1), a minimum mapping quality of 20 (−minMapQ 20), a minimum base quality of 20 (−minQ 20), discard sites where there is no data on at least 95% of the individuals (−minInd), skip tri-allelic sites (−skipTriallelic), and remove SNP sites with a p-value larger than 1e− 6 (−SNP_pval 1e-6). The ANGSD output beagle file was run through PCAngsd v0.95"

"Genotypes were called using ANGSD, specifying the same parameters as the PCA analyses with the additional parameters: write major and minor alleles and the genotype directly (−doGeno 5), estimate the posterior genotype probability based on the allele frequency as a prior (−doPost 1), use the reference allele as the major allele (−doMajorMinor 4), output as beagle likelihood file (−doGlf 2), and calculate allele frequencies assuming a fixed major allele and an unknown minor allele (−doMaf 2)"

- Fuller et al. 2020. Population genetics of the coral Acropora millepora: Toward genomic prediction of bleaching [link](https://science.sciencemag.org/content/sci/suppl/2020/07/15/369.6501.eaba4674.DC1/aba4674-Fuller-SM.pdf)

This science paper shows some interesting analyses with low-coverage sequencing data (~ 1.43 X). Maybe useful for my downstream analyses.

Not very informative in terms of Angsd running though.
 
- Calfee et al 2020 (Graham Coop lab). Selection and hybridization shaped the Africanized honey bee invasion of the Americas [link](https://www.biorxiv.org/content/10.1101/2020.03.17.994632v2.full.pdf)

"Using the software ANGSD, we identified a set of SNPs with minor allele frequency ≥ 5% in the combined sample based on read counts (-doMajorMinor 2 -doCounts 1). We excluded unplaced scaffolds (<5Mb total) and applied standard quality filters for SNP calling (base quality ≥ 20, mapping quality ≥ 30, total read depth ≤ 5500 (∼2x mean), and coverage across individuals ≥ 50%). "




