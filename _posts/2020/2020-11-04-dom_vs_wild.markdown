---
comments: true
title: Resequencing data formating and examination
date: '2020-11-04 12:00'
tags:
  - Reseq
  - vcf
  - Outflank
  - outlier 
categories:
  - WGS data analysis
---

#### LD pruning

- Prepare the masked dataset 

```sh
grep -v "^#" SNP.ORIGINAL.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf |wc -l
# Masked: 6413937 (contracted vcf file size: 17G)

vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.vcf --remove LM --recode --recode-INFO-all --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.noLM.vcf

#VCFtools - 0.1.15
#(C) Adam Auton and Anthony Marcketta 2009
#
#Parameters as interpreted:
#    --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.vcf
#    --exclude LM
#    --recode-INFO-all
#    --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.noLM.vcf
#    --recode

#Excluding individuals in 'exclude' list
#After filtering, kept 85 out of 90 Individuals
#Outputting VCF file...
#After filtering, kept 6413937 out of a possible 6413937 Sites
#Run Time = 2397.00 seconds
```
- Convert vcf to plink

```sh
plink --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.noLM.recode.vcf --recode --double-id --make-bed --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.noLM
```

- Process the bed file using (bigsnpr)[https://github.com/privefl/bigsnpr]

```R
library(bigsnpr)

bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.noLM.bed"
# this will create a .rds file
obj.bed <- bed(bedfile)
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.noLM.rds")
# See how it looks like
str(obj.bigSNP, max.level = 2, strict.width = "cut")
# List of 3
#  $ genotypes:Reference class 'FBM.code256' [package "bigstatsr"] with 15 fields
#   ..and 26 methods, of which 12 are  possibly relevant:
#   ..  add_columns, as.FBM, bm, bm.desc, check_dimensions, check_write_permissions, copy#envRefClass, initialize, initialize#FBM, save,
#   ..  show#envRefClass, show#FBM
#  $ fam      :'data.frame':  85 obs. of  6 variables:
#   ..$ family.ID  : chr [1:85] "CL_1" "CL_2" "CL_3" "CL_4" ...
#   ..$ sample.ID  : chr [1:85] "CL_1" "CL_2" "CL_3" "CL_4" ...
#   ..$ paternal.ID: int [1:85] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ maternal.ID: int [1:85] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ sex        : int [1:85] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ affection  : int [1:85] -9 -9 -9 -9 -9 -9 -9 -9 -9 -9 ...
#  $ map      :'data.frame':  6413937 obs. of  6 variables:
#   ..$ chromosome  : int [1:6413937] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ marker.ID   : chr [1:6413937] "1_18601" "1_18606" "1_18626" "1_18646" ...
#   ..$ genetic.dist: int [1:6413937] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ physical.pos: int [1:6413937] 18601 18606 18626 18646 18712 18716 18718 18733 18746 18778 ...
#   ..$ allele1     : chr [1:6413937] "T" "T" "A" "T" ...
#   ..$ allele2     : chr [1:6413937] "C" "C" "G" "G" ...
#  - attr(*, "class")= chr "bigSNP"

# Get aliases for useful slots
G   <- obj.bigSNP$genotypes

# Note that most of the algorithms of this package don’t handle missing values. 
# I used snp_fastImputeSimple() to impute missing values of genotyped variants, methods are
# "random" (sampling according to allele frequencies),
# "mean0" (rounded mean), 
# "mean2" (rounded mean to 2 decimal places),
# "mode" (most frequent call)
G2 <- snp_fastImputeSimple(G, method = c("random"))
G[, 10224] # some missing values
# [1]  1  0  1  0  0  0  1  1  1  0  0  1  1  1  1  1  1  0  1  1  1  0  1  0  0  0  1  1  1  1  0  0  0  0  0  2  1  0 NA  2  0  0  1  1  0  1
#[47]  0  1  1  1  0  1  1  0  1  0  0  2  2  2  2  0  2  1  1  0  0  1  0  1  0  1  0  1  0  1  1  0  1  0  1  1  1  2  2
G2[, 10224] # no missing values anymore
# [1] 1 0 1 0 0 0 1 1 1 0 0 1 1 1 1 1 1 0 1 1 1 0 1 0 0 0 1 1 1 1 0 0 0 0 0 2 1 0 0 2 0 0 1 1 0 1 0 1 1 1 0 1 1 0 1 0 0 2 2 2 2 0 2 1 1 0 0 1 0 1
# [71] 0 1 0 1 0 1 1 0 1 0 1 1 1 2 2
```

- LD clumping


