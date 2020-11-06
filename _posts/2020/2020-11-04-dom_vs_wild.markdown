---
comments: true
title: Outlier detection in dom vs wild constrasts
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

1) Prepare the masked domestic-wild constrasts vcf file. 

Total number of populations in masked vcf file: 16

HI - Maine_HI     
SM - Maine_SM    
HC - DelBay_HC    
CS - DelBay_CS    
CLP - ChesBay_CLP    
HCVA - ChesBay_HC_VA    
SL - Louisiana_SL    
CL - Louisiana_CL    
LM - Mexico_LM    
UMFS - Maine_Sel_UMFS         
NEH - DelBay_Sel_NEH             
DEBY - Ches_Sel_DEBY      
LOLA - Ches_Sel_LOLA         
OBOYS2 - Louisiana_Sel      
HG - Inbred_RU_HG           
NG - Inbred_REU_NG

Five populations used for the initital test:

SL - Louisiana wild line     
OBOYS2 - Louisiana selected line     
NEH - Delaware Bay selected NEH line      
DEBY - Chesapeake Bay selected line (initially from DB)       
CS - Cape Shore (Delaware Bay) wild line        


```sh
cat dom_wild.txt
# SL_1
# SL_2
# SL_3
# SL_4
# SL_5
# SL_6
# OBOYS2_1
# OBOYS2_2
# OBOYS2_3
# OBOYS2_4
# OBOYS2_5
# OBOYS2_6
# CS_1
# CS_2
# CS_3
# CS_5
# CS_6
# CS_7
# NEH_1
# NEH_2
# NEH_3
# NEH_4
# NEH_5
# NEH_6
# DEBY_1
# DEBY_2
# DEBY_3
# DEBY_4
# DEBY_5
# DEBY_6

grep -v "^#" SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.vcf |wc -l
# Masked: 6413937 (contracted vcf file size: 17G)

vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.vcf --keep dom_wild.txt --recode --recode-INFO-all --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.ind30

# check the genotype missing status of this file
vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.LDclumping.dom_wild.ind30.recode.vcf --max-missing 1 --recode --recode-INFO-all --out test

# filter out snps with low maf (0.05)
vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.ind30.recode.vcf --maf 0.05 --recode --recode-INFO-all --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing
#After filtering, kept 5533417 out of a possible 6413937 Sites
```

Features of the output vcf file:

1) no missing genotyps;     
2) markers with maf < 0.05 are removed;      
3) include 30 domestic vs wild contrast samples
4) after filtering, kept 5533417 out of 6413937 sites

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

# Note that most of the algorithms of this package don’t handle missing values. 
# I used snp_fastImputeSimple() to impute missing values of genotyped variants, methods are
# "random" (sampling according to allele frequencies),
# "mean0" (rounded mean), 
# "mean2" (rounded mean to 2 decimal places),
# "mode" (most frequent call)
#G2 <- snp_fastImputeSimple(G, method = c("random"))
#G[, 10224] # a missing value in the whole dataset (NA here)
# [1]  1  0  1  0  0  0  1  1  1  0  0  1  1  1  1  1  1  0  1  1  1  0  1  0  0  0  1  1  1  1  0  0  0  0  0  2  1  0 NA  2  0  0  1  1  0  1
#[47]  0  1  1  1  0  1  1  0  1  0  0  2  2  2  2  0  2  1  1  0  0  1  0  1  0  1  0  1  0  1  1  0  1  0  1  1  1  2  2
#G2[, 10224] # no missing values anymore in the whole dataset
# [1] 1 0 1 0 0 0 1 1 1 0 0 1 1 1 1 1 1 0 1 1 1 0 1 0 0 0 1 1 1 1 0 0 0 0 0 2 1 0 0 2 0 0 1 1 0 1 0 1 1 1 0 1 1 0 1 0 0 2 2 2 2 0 2 1 1 0 0 1 0 1
# [71] 0 1 0 1 0 1 1 0 1 0 1 1 1 2 2
```

3) LD clumping

See(here)[Why clumping should be preferred over pruning] for why I used clumping over pruning. In short, pruning may remove SNPs in a way that leave regions of the genome with no representative SNP at all.

```sh
# loop over different window size (in kb) for snp clumping
SIZE <- c(20, 50, 100, 200, 500)
ind_keeps = list()
# perform snp_clumping on masked dataset
for(i in seq(length(SIZE))){
  size_ = SIZE[i]
  ind.keep_ <- snp_clumping(
    G2,
    infos.chr = CHR,
    infos.pos = POS,
    ind.row = rows_along(G2),
    S = NULL,
    thr.r2 = 0.2,
    size = size_,
    exclude = NULL,
    ncores = 1
  )
  ind_keeps[[i]] = ind.keep_
}
```

Key parameters here:

infos.chr: vector of integers specifying each SNP's chromosome.
infos.pos: vector of integers specifying the physical position on a chromosome (in base pairs) of each SNP.      
thr.r2: threshold over the squared correlation between two SNPs. Default is 0.2.      
size: for one SNP, window size around this SNP to compute correlations. Default is 100 / thr.r2 for clumping (0.2 -> 500; 0.1 -> 1000; 0.5 -> 200). If not providing infos.pos (NULL, the default), this is a window in number of SNPs, otherwise it is a window in kb (genetic distance).   

Note SNP number of each VCF file are:     

Original: 4734919 (contracted vcf file size: 13G)      
Masked: 6413937 (contracted vcf file size: 17G)      

SNPs after clumping (note I only used 5 domestic and wild populations here)

|  Size (in KB)    | 20     | 50     | 100    | 200    | 500   |
|------------------|--------|--------|--------|--------|-------|
|  SNP retained    | 452784 | 310805 | 227423 | 162055 | 99466 |
|  %   retained    | 7.06%  | 4.85%  | 3.55%  | 2.53%  | 1.55% |

For initial test purpose, I used 50k bp as the thin-window for dataset clumping and ended up with 310805 SNPs for outflank analysis.

#### Outlier detection with outflank      

Outflank needs to use estimates of Fst that have not been corrected for sample size adjustments. In order to properly average Fst over loci, Outflank also needs the numerator and denominator of the Fst calculation.

When run Outflank, some parameters need to be set. 

1) LeftTrimFraction and RightTrimFraction (5% by default). This setting ensure that the lowest and highest Fst values to remove before estimating the shape of the neutral Fst distribution through likelihood.     
       
2) Hmin is the threshold for expected heterozygosity (10% by default). Loci with low He have distribution of Fst very different from loci with higher He and need to be removed.

3) NumberOfSaples is the number of populations sampled.

4) qthreshold sets the threshold for whether a locus is deemed an "outlier" (5% by default). These q-values are calculated based on the right-tail p-values for each locus.     

1) first trial using 5 populations and 50k bp thinned dataset

```R
################# start outflank #################
str(obj.bigSNP)
# create a new matrix to store the genotypes in bigsnpr
newm = matrix(nrow=length(obj.bigSNP$genotypes[,1]), ncol=length(obj.bigSNP$genotypes[1,]))
for(i in seq(length(obj.bigSNP$genotypes[,1]))){
  newm[i,]=obj.bigSNP$genotypes[i,]
}
# add pop information
obj.bigSNP$pop <- rep(c(seq(1,5)), each = 6)
# calculate FST on all the loci in our dataset.
my_fst <- MakeDiploidFSTMat(newm, locusNames = obj.bigSNP$map$physical.pos, popNames = obj.bigSNP$pop)
# Using OutFLANK() function to estimate the parameters on the neutral FST distribution
out_trim <- OutFLANK(my_fst[ind_keeps[[2]],], NumberOfSamples=5, qthreshold = 0.05, Hmin = 0.1)
str(out_trim)
# Check the fit
OutFLANKResultsPlotter(out_trim, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)
```

<img src="https://hzz0024.github.io/images/outflank/Fst_distribution.jpeg" alt="img" width="800"/>

```R
# check the fit of right tail:
OutFLANKResultsPlotter(out_trim , withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         TRUE, RightZoomFraction = 0.15, titletext = NULL)

```

<img src="https://hzz0024.github.io/images/outflank/Fst_righttail.jpeg" alt="img" width="800"/>

```R
hist(out_trim$results$pvaluesRightTail)
```

<img src="https://hzz0024.github.io/images/outflank/pvaluerighttail.jpeg" alt="img" width="800"/>

```R
# Using estimated neutral mean FST and df to calculate P-values for all loci
P1 <- pOutlierFinderChiSqNoCorr(my_fst, Fstbar = out_trim$FSTNoCorrbar, 
                                dfInferred = out_trim$dfInferred, qthreshold = 0.05, Hmin=0.1)
length(P1$OutlierFlag[P1$OutlierFlag==TRUE])
```

Number of outliers: 712647 (unrealistic number)

- what if change the population labeling

```sh
obj.bigSNP$pop <- rep(c(1,2,2,2,1), each = 6)
```

2) Now I have wild population one group and the domestic ones another. Still the results are not promising.

- Check the Fst fit

<img src="https://hzz0024.github.io/images/outflank/Fst_dist_2pop.jpeg" alt="img" width="800"/>

- Check the Fst fit on the right tail

<img src="https://hzz0024.github.io/images/outflank/rightail_2pop.jpeg" alt="img" width="800"/>

- p-value distribution on the right tail

<img src="https://hzz0024.github.io/images/outflank/p_value_2pop.jpeg" alt="img" width="800"/>

number of SNPs flaged as outliers: 601800


#### Outlier detection with pcadapt

Prepare the dataset, again I used the thinned data with 50k bp as a window size.

```sh
# extract the thinnd list with 310805 SNPs
vcftools --gzvcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.recode.vcf.gz --snps SNP_thinned.txt --recode --recode-INFO-all --out pcadapt
# convert into bed format
plink --vcf pcadapt.recode.vcf --double-id --make-bed --out pcadapt
```

Perform pcadapt analysis

1) load the data

```sh
# reading genotype data (“pcadapt”, “lfmm”, “vcf”, “bed”, “ped”, “pool”)
path_to_file <- "./ALL.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
```

2) pcadapt

```sh
################# start pcadapt #################
install.packages("pcadapt")
library("pcadapt")
path_to_file <- "./pcadapt.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
# Scree plot
x <- pcadapt(input = filename, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
plot(x, option = "screeplot")
```

<img src="https://hzz0024.github.io/images/outflank/scree_plot.jpeg" alt="img" width="800"/>

```sh
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("Louisiana_SL", 6),rep("Louisiana_Sel", 6), rep("DelBay_Sel_NEH", 6),rep("Ches_Sel_DEBY", 6), rep("DelBay_CS", 6))
print(poplist.names)
plot(x, option = "scores", pop = poplist.names)
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)
```
    
<img src="https://hzz0024.github.io/images/outflank/PC1_PC2.jpeg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/outflank/PC3_PC4.jpeg" alt="img" width="800"/>

```sh
# Computing the test statistic based on PCA
x <- pcadapt(filename, K = 3)
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
plot(x , option = "manhattan")
# Q-Q Plot
plot(x, option = "qqplot")
```

<img src="https://hzz0024.github.io/images/outflank/PC3_PC4.jpeg" alt="img" width="800"/>

```sh
# Histograms of the test statistic and of the p-values
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
```

<img src="https://hzz0024.github.io/images/outflank/p-value.jpeg" alt="img" width="800"/>