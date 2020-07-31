---
comments: true
title: CVreseq pcadapt
date: '2020-07-31 12:00'
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

Input: In pcadapt v4, the preferred format is now the PLINK “bed” format. Format “bed” is very compact, which stores each genotype using only 2 bits. One intesting feature of pcadapt v4 is that it's able to perform the outlier detection with pooled sequencing data. The pcadapt also accecpt missing values (e.g. 9). I used the plink to convert the vcf to bed format. 

### data description

We have three datasets, each with 12 individuals. They are:

LA-LSSL (SL) vs. LA-OBOY (OBOYS2)
Louisiana wild vs. selected line, source of oysters from selection from approx same environment (Oyster Bayou)

DB-HSCS (CS) vs. DB-NEHD (NEH)
Wild (Cape Shore, Delaware Bay wild high salinity) vs. selected (core line NEH)

DB-HSCS (CS) vs. CB-DEBY (DEBY)
Wild (Cape Shore, Delaware Bay wild high salinity) vs. Chesapeake Bay selected (initially from DB)

```sh
# extrac the 30 samples from vcf file
cat ALL
OBOYS2_1
OBOYS2_2
OBOYS2_3
OBOYS2_4
OBOYS2_5
OBOYS2_6
SL_1
SL_2
SL_3
SL_4
SL_5
SL_6
NEH_1
NEH_2
NEH_3
NEH_4
NEH_5
NEH_6
DEBY_1
DEBY_2
DEBY_3
DEBY_4
DEBY_5
DEBY_6
CS_1
CS_2
CS_3
CS_5
CS_6
CS_7

for pop in ALL; do
    vcftools --vcf CVreseq_chr.vcf --keep $pop --recode --recode-INFO-all --out $pop
    grep "^#" $pop'.recode.vcf' > $pop'.sort.vcf' && grep -v "^#" $pop'.recode.vcf' | \
    sort -V -k1,1 -k2,2n >> $pop'.sort.vcf'
done
```

### data conversion

```sh
# before conversion, a python script addID.py is used to add the SNP ID (CHR + "_" + POS) to the vcf file
fname = 'ALL.sort.vcf'
outname = fname + '.out'

idx = 0
with open(fname, 'r') as f, open(outname, 'w') as w:
    for l in f:
        if l.startswith('#'):
            pass
        else:
            idx += 1
            ss = l.split()
            chrom = ss[0]
            pos = ss[1]
            ID = chrom + '_' + pos
            assert ss[2] == '.'
            ss[2] = ID
            l = '\t'.join(ss)
            l += '\n'
        w.write(l)

python3 addID.py
# filter out SNPs with maf < 0.05 and call rate < 50%, for the ALL population pair
plink --vcf ALL.sort_ID.vcf --allow-extra-chr --biallelic-only --maf 0.05 --geno 0.5 --mind 0.5 --make-bed --out ALL_maf

log below,

PLINK v1.90b6.9 64-bit (4 Mar 2019)            www.cog-genomics.org/plink/1.9/
(C) 2005-2019 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to ALL_maf.log.
Options in effect:
  --allow-extra-chr
  --biallelic-only
  --geno 0.5
  --maf 0.05
  --make-bed
  --mind 0.5
  --out ALL_maf
  --vcf ALL.sort_ID.vcf

63991 MB RAM detected; reserving 31995 MB for main workspace.
--vcf: ALL_maf-temporary.bed + ALL_maf-temporary.bim + ALL_maf-temporary.fam
written.
334011 variants loaded from .bim file.
30 people (0 males, 0 females, 30 ambiguous) loaded from .fam.
Ambiguous sex IDs written to ALL_maf.nosex .
0 people removed due to missing genotype data (--mind).
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 30 founders and 0 nonfounders present.
Calculating allele frequencies... done.
0 variants removed due to missing genotype data (--geno).
13139 variants removed due to minor allele threshold(s)
(--maf/--max-maf/--mac/--max-mac).
320872 variants and 30 people pass filters and QC.
Note: No phenotypes present.
--make-bed to ALL_maf.bed + ALL_maf.bim + ALL_maf.fam ... done.

# LD pruning using plink
plink --bfile ALL_maf --indep-pairwise 50 5 0.5 --out ALL_tmp

log below,

PLINK v1.90b6.9 64-bit (4 Mar 2019)            www.cog-genomics.org/plink/1.9/
(C) 2005-2019 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to ALL_tmp.log.
Options in effect:
  --bfile ALL_maf
  --indep-pairwise 50 5 0.5
  --out ALL_tmp

63991 MB RAM detected; reserving 31995 MB for main workspace.
320872 variants loaded from .bim file.
30 people (0 males, 0 females, 30 ambiguous) loaded from .fam.
Ambiguous sex IDs written to ALL_tmp.nosex .
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 30 founders and 0 nonfounders present.
Calculating allele frequencies... done.
320872 variants and 30 people pass filters and QC.
Note: No phenotypes present.
Pruned 726 variants from chromosome 1, leaving 40803.
Pruned 715 variants from chromosome 2, leaving 34627.
Pruned 861 variants from chromosome 3, leaving 34244.
Pruned 816 variants from chromosome 4, leaving 31668.
Pruned 1432 variants from chromosome 5, leaving 61192.
Pruned 406 variants from chromosome 6, leaving 14338.
Pruned 487 variants from chromosome 7, leaving 19908.
Pruned 712 variants from chromosome 8, leaving 28229.
Pruned 954 variants from chromosome 9, leaving 39348.
Pruned 217 variants from chromosome 10, leaving 9189.
Pruning complete.  7326 of 320872 variants removed.
Marker lists written to ALL_tmp.prune.in and ALL_tmp.prune.out .

# extract the SNPs after pruning
plink --bfile ALL_maf --extract ALL_tmp.prune.in --make-bed --out ALL_prun

log below,

PLINK v1.90b6.9 64-bit (4 Mar 2019)            www.cog-genomics.org/plink/1.9/
(C) 2005-2019 Shaun Purcell, Christopher Chang   GNU General Public License v3
Logging to ALL_prun.log.
Options in effect:
  --bfile ALL_maf
  --extract ALL_tmp.prune.in
  --make-bed
  --out ALL_prun

63991 MB RAM detected; reserving 31995 MB for main workspace.
320872 variants loaded from .bim file.
30 people (0 males, 0 females, 30 ambiguous) loaded from .fam.
Ambiguous sex IDs written to ALL_prun.nosex .
--extract: 313546 variants remaining.
Using 1 thread (no multithreaded calculations invoked).
Before main variant filters, 30 founders and 0 nonfounders present.
Calculating allele frequencies... done.
313546 variants and 30 people pass filters and QC.
Note: No phenotypes present.
--make-bed to ALL_prun.bed + ALL_prun.bim + ALL_prun.fam ... done.
```

### Using pcadapt to detect local adaptation

1) load the data

```R
# reading genotype data (“pcadapt”, “lfmm”, “vcf”, “bed”, “ped”, “pool”)
path_to_file <- "./ALL.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
```

2) Choosing the number K of Principal Components

```R
# 2.1 Scree plot
x <- pcadapt(input = filename, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
plot(x, option = "screeplot")
plot(x, option = "screeplot", K = 2)
```

<img src="https://hzz0024.github.io/images/pcadapt/k_10.jpeg" alt="img" width="800"/>

k = 5 maybe the best k

```R
#2.2 Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.

poplist.names <- c(rep("POP1", 6),rep("POP2", 6), rep("POP3", 6),rep("POP4", 6), rep("POP5", 6))
print(poplist.names)
plot(x, option = "scores", pop = poplist.names)
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)
```

<img src="https://hzz0024.github.io/images/pcadapt/pca.jpeg" alt="img" width="800"/>

3) Computing the test statistic based on PCA

```R
x <- pcadapt(filename, K = 5)

summary(x)
                Length  Class  Mode   
scores              150 -none- numeric
singular.values       5 -none- numeric
loadings        1567730 -none- numeric
zscores         1567730 -none- numeric
af               313546 -none- numeric
maf              313546 -none- numeric
chi2.stat        313546 -none- numeric
stat             313546 -none- numeric
gif                   1 -none- numeric
pvalues          313546 -none- numeric
pass             313546 -none- numeric

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
[1] 0
```

4) Graphical tools
```R
# 4.1 Manhattan Plot
plot(x , option = "manhattan")
```

<img src="https://hzz0024.github.io/images/pcadapt/manhattan.jpeg" alt="img" width="800"/>

```R
# 4.2 Q-Q Plot
plot(x, option = "qqplot")
```

<img src="https://hzz0024.github.io/images/pcadapt/qqplot.jpeg" alt="img" width="800"/>

This plot confirms that most of the p-values follow the expected uniform distribution. However, the smallest p-values are smaller than expected confirming the presence of outliers.


```R
# 4.3 Histograms of the test statistic and of the p-values
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
```

<img src="https://hzz0024.github.io/images/pcadapt/qqplot.jpeg" alt="img" width="800"/>

An histogram of p-values confirms that most of the p-values follow an uniform distribution. The excess of small p-values indicates the presence of outliers.

```R
# The presence of outliers is also visible when plotting a histogram of the test statistic 𝐷𝑗.
plot(x, option = "stat.distribution")
```

<img src="https://hzz0024.github.io/images/pcadapt/qqplot.jpeg" alt="img" width="800"/>

The presence of outliers is also visible when plotting a histogram of the test statistic 𝐷𝑗 (not sure how to interpreter it yet)

5) Choosing a cutoff for outlier detection

```R
# q-values
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.01
outliers <- which(qval < alpha)
length(outliers)

[1] 2590
# Benjamini-Hochberg Procedure
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)

[1] 2590
# Bonferroni correction
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)

[1] 667
```


  

