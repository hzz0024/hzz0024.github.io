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

- Features of the vcf file:

1) no missing genotyps;     
2) markers with maf < 0.05 are removed;      
3) only include 30 domestic vs wild contrast samples
4) after filtering, kept 5,533,417 out of 6,413,937 sites

|                                 Vcf file                                         |   No. SNPs  |
|----------------------------------------------------------------------------------|-------------|
|SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing|   5533417   |

2) Convert vcf to plink

```sh
plink --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.recode.vcf --double-id --make-bed --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing
```

- Process the bed file using (bigsnpr)[https://github.com/privefl/bigsnpr]

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

#### Outlier detection with outflank      

Outflank needs to use estimates of Fst that have not been corrected for sample size adjustments. In order to properly average Fst over loci, Outflank also needs the numerator and denominator of the Fst calculation.

When run Outflank, some parameters need to be set. 

1) LeftTrimFraction and RightTrimFraction (5% by default). This setting ensure that the lowest and highest Fst values to remove before estimating the shape of the neutral Fst distribution through likelihood.     
       
2) Hmin is the threshold for expected heterozygosity (10% by default). Loci with low He have distribution of Fst very different from loci with higher He and need to be removed.

3) NumberOfSaples is the number of populations sampled.

4) qthreshold sets the threshold for whether a locus is deemed an "outlier" (5% by default). These q-values are calculated based on the right-tail p-values for each locus.     


