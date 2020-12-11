---
comments: true
title: Outlier detection and zoom in
date: '2020-11-19 12:00'
tags:
  - Reseq
  - vcf
  - pcadapt
  - outlier
  - pi
  - dxy 
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
 
|                                 Vcf file                                             |   No. SNPs    |
|--------------------------------------------------------------------------------------|---------------|
|SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.vcf|   5,533,417   |

#### Outlier detection with thinned SNPs (Best Practice)

See Supplementary Material of Katies's paper for method details [The effect of neutral recombination variation on genome scans for selection](https://gsajournals.figshare.com/articles/Supplemental_Material_for_Lotterhos_2019/7973438)

In short, for each SNP, the PCAdapt computes a Mahalanobis distance on a vector of z-scores that corresponds to the z-scores obtained when regressing a SNP by the K PCs. The neutralparameterization in this case are the PC axes that describe population genetic structure. P-values for each locus are obtained from a chi-squared distribution with K degrees of freedom......The best practice was implemented by using the quasi-independent thinned set of SNPs to estimate the PC axes, and then using these estimate −log10P-values from the regression for all SNPs with MAF > 0.05*

See Katie's Github for code details [https://github.com/TestTheTests/TTT_RecombinationGenomeScans/blob/master/src/b_Proc_Sims.R](https://github.com/TestTheTests/TTT_RecombinationGenomeScans/blob/master/src/b_Proc_Sims.R)

```sh
plink --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.thin1K.vcf.recode.vcf --double-id --make-bed --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.thin1K
```

Figures below are PCA plots produced by first few PCs using thinned dataset 

- window size of 500K bp

<img src="https://hzz0024.github.io/images/pcadapt/500K_thin_screeplot.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/pcadapt/500K_PCs.jpg" alt="img" width="800"/>

- window size of 5K bp

<img src="https://hzz0024.github.io/images/pcadapt/5K_thin_screeplot.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/pcadapt/5K_PCs.jpg" alt="img" width="800"/>

- window size of 1K bp

<img src="https://hzz0024.github.io/images/pcadapt/1K_thin_screeplot.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/pcadapt/1K_PCs.jpg" alt="img" width="800"/>

Number of outlier identified with Best Practice (bonferroni p-value < 0.05)

|  PC used         | No. outliers (win = 500K)| No. outliers (win = 5K)| No. outliers (win = 1K)| 
|------------------|--------------------------|------------------------|------------------------|
|  1               |           32             |      2498              |      2243              |
|  1-3             |          155             |      1277              |      850               |
|  1-5             |           83             |      50                |      68                |


- Best practice with K=3 and window size of 5K bp

<img src="https://hzz0024.github.io/images/pcadapt/Mahattan_all_BP_outlier_PC1-3_5K.jpg" alt="img" width="800"/>

- Best practice with K=5 and window size of 5K bp

<img src="https://hzz0024.github.io/images/pcadapt/Mahattan_all_BP_outlier_PC1-5_5K.jpg" alt="img" width="800"/>

- Best practice with PC1 only and window size of 5K bp

<img src="https://hzz0024.github.io/images/pcadapt/Mahattan_all_BP_outlier_PC1_5K.jpg" alt="img" width="800"/>

- Best practice with K=3 and window size of 1K bp

<img src="https://hzz0024.github.io/images/pcadapt/Mahattan_all_BP_outlier_PC1-3_1K.jpg" alt="img" width="800"/>

- Best practice with K=5 and window size of 1K bp

<img src="https://hzz0024.github.io/images/pcadapt/Mahattan_all_BP_outlier_PC1-5_1K.jpg" alt="img" width="800"/>

- Best practice with PC1 only and window size of 1K bp

<img src="https://hzz0024.github.io/images/pcadapt/Mahattan_all_BP_outlier_PC1_1K.jpg" alt="img" width="800"/>

---

#### Explore the 1K window for outlier identification

Note: inversions.masked.bed is the the inversions called by [Delly](https://github.com/dellytools/delly) (the SV variant caller used for the CNVs)

```sh
# extract inversions only
vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.vcf --bed inversions.masked.bed --recode --recode-INFO-all --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.invers
plink --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.invers.recode.vcf --double-id --make-bed --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.invers
# extract non-inversion from the vcf file
vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.vcf --exclude-bed inversions.masked.bed --recode --recode-INFO-all --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers
plink --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers.recode.vcf --double-id --make-bed --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers
```

Code used to eliminate the chr5 and chr6 inversions

```sh
plink --bfile SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers --exclude range chr1.txt --double-id --make-bed --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers_nochr56invers

cat chr1.txt

5 60600000 80200000 chr5
6 29900000 44500000 chr6

# to keep the chr5 and 6 inversion
plink --bfile SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers --exclude range chr.txt --double-id --make-bed --out SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.noinvers_chr56invers

cat chr.txt

1 1 65644315 chr1
2 1 61751371 chr2
3 1 77053199 chr3
4 1 58754694 chr4
5 1 60600000 chr51
5 80200000 98632791 chr52
6 1 29900000 ch61
6 44500000 51169664 chr62
7 1 57829119 chr7
8 1 75928414 chr8
9 1 104144669 chr9
10 1 32637058 chr10
```

#### Zoom in outliers

1) Prepare the VCF files

```sh
for pop in SL_OBOYS2 CS_NEH CS_DEBY SL_OBOYS2_no_outlier; do
    vcftools --vcf SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.vcf --keep $pop --recode --recode-INFO-all --out $pop
done
```