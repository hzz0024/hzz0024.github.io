---
comments: true
title: Outlier identification with different maf filtering settings
date: '2020-10-14 12:00'
tags:
  - DelBay19
  - ouliter
  - WGS
  - MAF
  - Fisher
categories:
  - WGS data analysis
---

Different maf filtering parameters (0.01, 0.05 and 0.2) were set to reproduce the datasets. Below are the number of SNPs,

|               | maf < 0.01 | maf < 0.01| maf < 0.01| 
|---------------|------------|-----------|-----------|
| Total SNPs    | 566214907  | 566214907 |566214907  |
| Number of SNPs| 4969475    | 1934038   | 627518    | 

After then, the global snp list from each dataset was used to call SNP in each independent population. For example,

```sh
module load angsd/0.931
###this script will work on bamfiles by population and calculate saf  & maf
# maybe edit
target="CH"
NB_CPU=20 #change accordingly
REGIONS="-rf chr_list.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /scratch/hzz0024/DelBay19_Sep/01_scripts/01_config.sh

angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doMajorMinor 3 -anc $ANC -minQ 20 -b $CH $REGIONS -sites ALL_sites_all_maf0.01_pctind0.7_maxdepth3dv_snplist_4col_cv30 -out "/scratch/hzz0024/DelBay19_Sep/05_saf_maf_by_pop/"$target"_maf0.01_pctind"$PERCENT_IND"_cv30"
```


|   Test                 | fdr < 0.1 | fdr < 0.05| fdr < 0.01| fdr < 0.1 | fdr < 0.05| fdr < 0.01|
|------------------------|-----------|-----------|-----------|-----------|-----------|-----------|
|                        | Z method  | Z method  | Z method  | Fisher    | Fisher    |  Fisher   |
|REF-CH-SR-HC(maf<0.01)  |           |           |           |           |           |           |
|REF-CH-NB-HC(maf<0.01)  |           |           |           |           |           |           |
|SR-REF-COH-ARN(maf<0.01)|           |           |           |           |           |           |
|REF-CH-SR-HC(maf<0.05)  |     11    |  4        |      0    |     41    | 10        |      1    |
|REF-CH-NB-HC(maf<0.05)  |     31    |  8        |      1    |     32    | 16        |      6    |
|SR-REF-COH-ARN(maf<0.05)|      0    |  0        |      0    |     20    | 0         |      0    |
|REF-CH-SR-HC(maf<0.2)   |     44    |  19       |      2    |     106   | 30        |      0    |
|REF-CH-NB-HC(maf<0.2)   |     36    |  12       |      6    |     133   | 19        |      9    |
|SR-REF-COH-ARN(maf<0.2) |     13    |  0        |      0    |      48   | 13        |      0    |

Any shared loci between maf 0.05 and 0.01?

REF-CH-SR-HC(maf<0.05) and REF-CH-SR-HC(maf<0.2) 

|               | z     | fisher | 
|---------------|-------|--------|
| fdr < 0.1     | 10    |   10   |
| fdr < 0.05    | 4     |   4    |

REF-CH-NB-HC(maf<0.05) and REF-CH-NB-HC(maf<0.2)

|               | z     | fisher | 
|---------------|-------|--------|
| fdr < 0.1     | 18    |   18   |
| fdr < 0.05    | 8     |   8    |
