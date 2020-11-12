---
comments: true
title: Outlier detection updates
date: '2020-11-12 12:00'
tags:
  - Reseq
  - vcf
  - pcadapt
  - outlier 
categories:
  - WGS data analysis
---

#### Vcf preparation

- Populations:

SL - Louisiana wild line     
OBOYS2 - Louisiana selected line     
NEH - Delaware Bay selected NEH line      
DEBY - Chesapeake Bay selected line (initially from DB)       
CS - Cape Shore (Delaware Bay) wild line  

Among them, DEBY originally from Delaware Bay (source material), brought to Chesapeake Bay (CB) for selection in 1980s, some bottlenecks in hatchery, DEBY kept in different environment in CB since 80s. CS is sourced from Cape Shore at Delaware Bay (wild high salinity).     

- Features of the vcf file:

1) no missing genotyps;     
2) markers with maf < 0.05 are removed;      
3) only include 30 domestic vs wild contrast samples
4) after filtering, kept 5,533,417 out of 6,413,937 sites

|                                 Vcf file                                         |   No. SNPs    |
|----------------------------------------------------------------------------------|---------------|
|SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing|   5,533,417   |

2) Convert vcf to plink

```sh
plink --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.recode.vcf --double-id --make-bed --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing
```

- Process the bed file using [bigsnpr](https://github.com/privefl/bigsnpr)

```R
library(bigsnpr)

bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.LDclumping.dom_wild.bed"
# this will create a .rds file
snp_readBed(bedfile)
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.LDclumping.dom_wild.rds")
# See how it looks like
str(obj.bigSNP, max.level = 2, strict.width = "cut")

# Get aliases for useful slots
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
```

#### Outlier detection with thinned SNPs (best practice)

See Supplementary Material of Katies's paper for method details [The effect of neutral recombination variation on genome scans for selection](https://gsajournals.figshare.com/articles/Supplemental_Material_for_Lotterhos_2019/7973438)

In short, for each SNP, the PCAdapt computes a Mahalanobis distance on a vector of z-scores that corresponds to the z-scores obtained when regressing a SNP by the K PCs. The neutralparameterization in this case are the PC axes that describe population genetic structure. P-values for each locus are obtained from a chi-squared distribution with K degrees of freedom......The best practice was implemented by using the quasi-independent thinned set of SNPs to estimate the PC axes, and then using these estimate −log10P-values from the regression for all SNPs with MAF > 0.05*

See Katie's Github for code details [https://github.com/TestTheTests/TTT_RecombinationGenomeScans/blob/master/src/b_Proc_Sims.R](https://github.com/TestTheTests/TTT_RecombinationGenomeScans/blob/master/src/b_Proc_Sims.R)

Figures below are PCA plots produced by first few PCs in the thinned dataset

<img src="https://hzz0024.github.io/images/outflank/PC1_PC2.jpeg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/outflank/PC3_PC4.jpeg" alt="img" width="800"/>


|  PC used         | No. outliers (bonferroni)| 
|------------------|--------------------------|
|  1               |           32             |
|  1-3             |          155             | 
|  1-5             |           83             | 

- Best practice with K=3

<img src="https://hzz0024.github.io/images/pcadapt/Mahattan_BP_PC1-3.jpg" alt="img" width="800"/>

- Best practice with K=5

<img src="https://hzz0024.github.io/images/pcadapt/Mahattan_BP_PC1-5.jpg" alt="img" width="800"/>

- Best practice with PC1 only

<img src="https://hzz0024.github.io/images/pcadapt/Mahattan_BP_PC1.jpg" alt="img" width="800"/>

#### Outlier detection with all SNPs (naive method)

The nave approach was evaluated using the −log10P-values that resulted from running the algorithm on all SNPs with MAF > 0.05. The p-value peaks in the naive results should occur at the inversion regions in the genome (chr 5 and 6).

|  PC used         | No. outliers (bonferroni)| 
|------------------|--------------------------|
|  1               |           2512           |
|  1-3             |           27408          |
|  1-5             |           46243          |

- Naive with K=3

<img src="https://hzz0024.github.io/images/pcadapt/Mahattan_naive_PC1-3.jpg" alt="img" width="800"/>

- Naive with K=5

<img src="https://hzz0024.github.io/images/pcadapt/Mahattan_naive_PC1-5.jpg" alt="img" width="800"/>

- Naive with PC1 only

<img src="https://hzz0024.github.io/images/pcadapt/Mahattan_naive_PC1.jpg" alt="img" width="800"/>
     
---

The results above is generated from the whole dom-wild popultions (5 pop with 30 individuals). 

Now focus on the pairwise comparsions.

```sh
vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.vcf --maf 0.05 --max-missing 1 --keep CS_NEH --recode --recode-INFO-all --out CS_NEH
vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.vcf --maf 0.05 --max-missing 1 --keep SL_OBOYS2 --recode --recode-INFO-all --out SL_OBOYS2
vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.vcf --maf 0.05 --max-missing 1 --keep CS_DEBY --recode --recode-INFO-all --out CS_DEBY
```
