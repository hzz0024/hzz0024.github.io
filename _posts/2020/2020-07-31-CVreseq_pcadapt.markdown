---
comments: true
title: DelBay19 single generation selection (SGS) test
date: '2020-07-27 12:00'
tags:
  - CVreseq
  - pcadapt
  - ouliter
  - WGS
categories:
  - WGS data analysis
---

### Pcadapt description

The package pcadapt is a user-friendly tool to detect signs of local adaptation in genetic data. The pcadapt method comes with four major steps:

1) computes the Principal Component Analysis (PCA) of a scaled genotype matrix   
2) regresses all variants onto the resulting PCs to get a matrix of Z-scores (i.e., one Z-score for each variant and each PC)   
3) computes robust Mahalanobis distances of these Z-scores to integrate all PCA dimensions in one multivariate distance for each variant      
4) these distances approximately follow a chi-squared distribution, which enables derivation of one P-value for each genetic variant   

Input: In pcadapt v4, the preferred format is now the PLINK “bed” format. Format “bed” is very compact, which stores each genotype using only 2 bits. The pcadapt also accecpt missing values (e.g. 9). I used the plink to convert the vcf to bed format. 

```sh
plink --vcf DB_1.sort.vcf --biallelic-only --maf 0.05 --geno 0.5 --mind 0.5 --make-bed --out DB_1
# the --geno 0.5 filter out snps with more than 50% missing genotypes
```

One intesting feature of pcadapt v4 is that it's able to perform the outlier detection with pooled sequencing data.

### Using pcadapt to detect local adaptation

1) load the data

```R
# reading genotype data (“pcadapt”, “lfmm”, “vcf”, “bed”, “ped”, “pool”)
path_to_file <- "./DB_1.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
```

2) Choosing the number K of Principal Components

```R
# 2.1 Scree plot
x <- pcadapt(input = filename, K = 5)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
plot(x, option = "screeplot")
plot(x, option = "screeplot", K = 2)

#2.2 Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.

poplist.names <- c(rep("POP1", 6),rep("POP2", 6))
print(poplist.names)
plot(x, option = "scores", pop = poplist.names)
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)
```

3) Computing the test statistic based on PCA

```R
x <- pcadapt(filename, K = 2)

summary(x)
>                Length Class  Mode   
>scores              24 -none- numeric
>singular.values      2 -none- numeric
>loadings        564134 -none- numeric
>zscores         564134 -none- numeric
>af              282067 -none- numeric
>maf             282067 -none- numeric
>chi2.stat       282067 -none- numeric
>stat            282067 -none- numeric
>gif                  1 -none- numeric
>pvalues         282067 -none- numeric
>pass            280464 -none- numeric

scores is a (n,K) matrix corresponding to the projections of the individuals onto each PC.   
singular.values is a vector containing the K ordered square root of the proportion of variance explained by each PC.   
loadings is a (L,K) matrix containing the correlations between each genetic marker and each PC.   
zscores is a (L,K) matrix containing the 𝑧-scores.   
af is a vector of size L containing allele frequencies of derived alleles where genotypes of 0 are supposed to code for homozygous for the reference allele.   
maf is a vector of size L containing minor allele frequencies.   
chi2.stat is a vector of size L containing the rescaled statistics stat/gif that follow a chi-squared distribution with K degrees of freedom.   
gif is a numerical value corresponding to the genomic inflation factor estimated from stat.   
pvalues is a vector containing L p-values.   
pass A list of SNPs indices that are kept after exclusion based on the minor allele frequency threshold.   
stat is a vector of size L containing squared Mahalanobis distances by default.   
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
> 1603
```

4) Graphical tools
```R
# 4.1 Manhattan Plot
plot(x , option = "manhattan")

# Q-Q Plot
plot(x, option = "qqplot")
```


This plot confirms that most of the p-values follow the expected uniform distribution. However, the smallest p-values are smaller than expected confirming the presence of outliers.

5) Histograms of the test statistic and of the p-values















```R
# SGS test
Rscript --verbose SGS_DelBay19.R 

# check the result

cat p_values_local.txt | wc -l
# the number of potential outliers
4125 
# extract the allele frequency values from N1 and N2
python3 extract.py
# calculate the actual delta p values for 4125 SNPs
Rscript deltaP_act.R -d /Users/ryan/Documents/Ryan_workplace/DelBay19_HG/10_SGS/SGS/local_theat_SGS_results -p CHR_maf0.05_pctind0.7_cv30.mafs.extracted -q CH_maf0.05_pctind0.7_cv30.mafs.extracted -t 4125 -o obs_deltap.output
```
<img src="https://hzz0024.github.io/images/SGS/allele_1.jpeg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/SGS/allele_2.jpeg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/SGS/deltap.jpeg" alt="img" width="800"/>

Remained question: What is appropriate window size? It may be difficult to know without trying and comparing several (e.g. on one chromosome). How does StDev(theta) vary with window size. Maybe we want the window size where variance of theta is greatest?

Are there any common shared outliers between SGS results and wild transect comparsions? How to determine the potential SNPs of selection in the wild transect?
  

