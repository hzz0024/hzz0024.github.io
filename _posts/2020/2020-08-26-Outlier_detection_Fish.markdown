---
comments: true
title: DelBay19 Fishers' exact test
date: '2020-08-26 12:00'
tags:
  - DelBay19
  - ouliter
  - WGS
  - Fisher's exact
categories:
  - WGS data analysis
---

### Fisher's exact tests (two draws from a common pool approach)

I performed the tests in both challenge vs. refernce and wild contrasts. 

First is the challenge (CH) vs. reference (REF) group, 

P-value histogram for challenge (CH) vs. reference (REF)

<img src="https://hzz0024.github.io/images/outlier/fisher_exact_chr_pvalue_hist.jpeg" alt="img" width="800"/>

An example of wild transect comparison is the HC (Hope Creek) vs New Beds (NB) group.

P-value histogram for HC (Hope Creek) vs New Beds (NB) 

<img src="https://hzz0024.github.io/images/outlier/fisher_exact_wild_pvalue_hist.jpeg" alt="img" width="800"/>

The detailed scripts are located in the DelBay_project/R_scripts/Fisher_exact/Fish_exact_HG.R

--- 

### Results 

|     Method       | No. SNPs | p-value < 0.05 | fdr < 0.2 | fdr < 0.1  | fdr < 0.05 |
|------------------|----------|----------------|------------|-----------|------------|
| CH vs. REF       | 1732036  |   152978       |   3829     |    478    |    96      |
| HC vs. ARN       | 2083615  |   168062       |   1779     |    193    |    54      |
| HC vs. COH       | 2083615  |   165714       |   1779     |    171    |    62      |
| HC vs. SR        | 2083615  |   169239       |   2320     |    271    |    2       |
| HC vs. NB        | 2083615  |   175197       |   2796     |    442    |    61      |
| ARN vs. COH      | 2083615  |   164884       |   1334     |    102    |    14      |
| ARN vs. SR       | 2083615  |   167113       |   1796     |    141    |    20      |
| ARN vs. NB       | 2083615  |   172733       |   2624     |    318    |    35      |
| COH vs. SR       | 2083615  |   169168       |   1710     |    285    |    63      | 
| COH vs. NB       | 2083615  |   172417       |   2656     |    403    |    79      | 
| SR vs. NB        | 2083615  |   174558       |   3238     |    417    |    33      | 

### Shared outliers between challenge-control and wild contrasts

|Group compared    | fdr < 0.2  | fdr < 0.1  | fdr < 0.05 |
|------------------|------------|------------|------------|
|CH_REF - HC_COH   |      5     |      0     |      0     |    
|CH_REF - HC_ARN   |      4     |      0     |      0     | 
|CH_REF - HC_SR    |      9     |      0     |      0     | 
|CH_REF - HC_NB    |      4     |      0     |      0     | 
|CH_REF - ARN_COH  |      3     |      0     |      0     | 
|CH_REF - ARN_SR   |      4     |      0     |      0     | 
|CH_REF - ARN_NB   |      1     |      0     |      0     | 
|CH_REF - COH_SR   |      2     |      0     |      0     | 
|CH_REF - COH_NB   |      4     |      0     |      0     | 
|CH_REF - SR_NB    |      8     |      0     |      0     | 

### Plot the delta_p against start p

- CH vs. REF outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 3829)

<img src="https://hzz0024.github.io/images/Fish/CH_REF_fdr02.jpg" alt="img" width="800"/>

- CH vs. REF randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/CH_REF_fdr02_random_sample.jpg" alt="img" width="800"/>

In the figures above I put some lines to illustrate why delta_p against p0 have shown such patterns.     
1) it makes sense that all the dots should be located above the y= - x line. That is, for a given p0 with negative delta_p, p0 + delta_p >= 0.    
2) it is weird to see that dots are constrained along the y= -2x + 1 and y -2x + 0.1 in the randomly sampling result. Theoretically, there is no constrain for the delta_p distribution. Let us take a look at the wild contrasts then. 

- HC vs. NB outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 2796)

<img src="https://hzz0024.github.io/images/Fish/HC_NB_fdr02.jpg" alt="img" width="800"/>

- HC vs. NB randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/HC_NB_fdr02_random_sample.jpg" alt="img" width="800"/>   

- ARN vs. COH outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 1334)

<img src="https://hzz0024.github.io/images/Fish/ARN_COH_fdr02.jpg" alt="img" width="800"/>

- ARN vs. COH randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/ARN_COH_fdr02_random_sample.jpg" alt="img" width="800"/> 

In the wild transect comparsions, however, there is no such constrains for the data distribution. This issue may come from two reasons. First is the MAF setting in the Angsd run. I set -doMaf as 0.05 for initial SNP calling across the whole populations (which creates a SNP list for downstream individual population running). The retalivaly high MAF setting may bias the allele frequency distribution in the individual population. Second is that I intially created the SNP lists seperately for the CH-REF and wild contrasts. Given that CH-REF has smaller sample size (N=98), I expect that more SNPs in the challenge or reference population are constrained by this MAF setting.

### Relationship between p0 and p1 in potential outliers

- CH vs. REF outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 3829)

<img src="https://hzz0024.github.io/images/Fish/CH_REF_fdr02_p0_p1.jpg" alt="img" width="800"/>

- CH vs. REF randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/CH_REF_fdr02_p0_p1_random_sample.jpg" alt="img" width="800"/>

- HC vs. NB outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 2796)

<img src="https://hzz0024.github.io/images/Fish/HC_NB_fdr02_p0_p1.jpg" alt="img" width="800"/>

- HC vs. NB randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/HC_NB_fdr02_p0_p1_random_sample.jpg" alt="img" width="800"/>

- HC vs. ARN outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 1779)

<img src="https://hzz0024.github.io/images/Fish/HC_ARN_fdr02_p0_p1.jpg" alt="img" width="800"/>

- HC vs. ARN randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/HC_ARN_fdr02_p0_p1_random_sample.jpg" alt="img" width="800"/>

- ARN vs. COH outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 1334)

<img src="https://hzz0024.github.io/images/Fish/ARN_COH_fdr02_p0_p1.jpg" alt="img" width="800"/>

- ARN vs. COH randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/ARN_COH_fdr02_p0_p1_random_sample.jpg" alt="img" width="800"/>

- COH vs. SR outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 1710)

<img src="https://hzz0024.github.io/images/Fish/COH_SR_fdr02_p0_p1.jpg" alt="img" width="800"/>

- COH vs. SR  randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/COH_SR_fdr02_p0_p1_random_sample.jpg" alt="img" width="800"/>

- SR vs. NB outliers identified from the Fisher's exact test with fdr < 0.2 ( n = 3238)

<img src="https://hzz0024.github.io/images/Fish/SR_NB_fdr02_p0_p1.jpg" alt="img" width="800"/>

- SR vs. NB randomly sampling the same amount of SNPs 

<img src="https://hzz0024.github.io/images/Fish/SR_NB_fdr02_p0_p1_random_sample.jpg" alt="img" width="800"/>

---

Notes:

- Angsd command example

```sh
module load angsd/0.931
###this script will work on bamfiles by population and calculate saf  & maf
# maybe edit
target="CH"
NB_CPU=20 #change accordingly
REGIONS="-rf chr_list.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /scratch/hzz0024/DelBay19_July/01_scripts/01_config.sh
N_IND=$(wc -l $CH | cut -d " " -f 1)
MIN_IND=$(($N_IND*7/10))

echo "Ouput can be used for depth evaluation with all individuals listed in "$CH
echo "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"
echo "filter on allele frequency = "$MIN_MAF

angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doMajorMinor 3 \
      -anc $ANC -remove_bads 1 -minMapQ 30 -minQ 20 -b $CH \
      -sites CHR_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 \
      -out "/scratch/hzz0024/DelBay19_July/05_saf_maf_by_pop/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"PERCENT_IND"_cv30"
```

- Check outliers

Here I looked at one potential outlier in CH_REF and HC_NB comparison, and found four shared outliers with fdr < 0.2.

```sh
intersect(CH_REF_list, HC_NB_list)
[1] "NC_035780.1_32280271" "NC_035780.1_32280434" "NC_035784.1_18028707" "NC_035786.1_8772371"
```

Allele frequency information from Angsd maf files

```sh
>NC_035784.1_18028707
# CH vs REF
> CH
chromo  position  major minor anc knownEM nInd
NC_035784.1 18028707  T C T 0.366273  45
> REF
NC_035784.1 18028707  T C T 0.110646  37

# Wild contrasts
>HC
NC_035784.1 18028707  T C T 0.338600  41
>NB
NC_035784.1 18028707  T C T 0.092032  40

>NC_035786.1_8772371
# CH vs REF
> CH
chromo  position  major minor anc knownEM nInd
NC_035786.1 8772371 C A A 0.494092  40
> REF
NC_035786.1 8772371 C A A 0.204234  43

# Wild contrasts
>HC
NC_035786.1 8772371 C A A 0.130250  40
>NB
NC_035786.1 8772371 C A A 0.497266  36

>NC_035780.1_32280271
# CH vs REF
> CH
chromo  position  major minor anc knownEM nInd
NC_035780.1 32280271  C T C 0.496700  38
> REF
NC_035780.1 32280271  C T C 0.172002  35

# Wild contrasts
>HC
NC_035780.1 32280271  C T C 0.519852  40
>NB
NC_035780.1 32280271  C T C 0.194860  36
```
