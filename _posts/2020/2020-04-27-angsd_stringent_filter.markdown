---
comments: true
title: DelBay19_stringent_filter
date: '2020-05-15 12:00'
tags:
  - DelBay19
  - angsd
  - filter
  - WGS
  - paralogs
categories:
  - WGS data analysis
---

In this post I'd like to evaluate how factors such as minInd, minMapQ, or paralogs impact the SNP counts and summary statistics produced from the low-coverage data. This post has two parts, in part I I only focus on the SNPs generated from chromosome 5 (chr 5), and examine how different parameters or paralogs impact the SNPs counts and fst values.  

In analysis part II, the GENOME-WIDE SNPs from chr 1-10 are used for two analyses: fst qqplot and PCA. I draw plots for fst distribution and see how does the data fit the exponential distribution. I also make several PCA plots based on different settings of minInd: 50%, 60%, 70% and 80%. By doing this I want to test whether the individual outliers are resulted from higher reads (therefore lead to higher genotype call rates) or not. 

---
PART I

The control run here used -minQ 20 -minMapQ 20 -minInd 25 (ch) or 24 (ref) - corresponding 50% minInd rate. 

### SNP counts summary

|     	     |Total sites (ch)|Sites retained (ch)|Total sites (ref)|Sites retained (ref)| Unweight Fst   |  Weight Fst       |  
| -----------|----------------|-------------------|-----------------|--------------------|----------------|-------------------|
|   minInd   |                |                   |                 |                    |                |                   |
|   50%      | 84791110       |    41680598       | 85218526        |    43137643        |    0.000733    |      0.001047     |
|   60%      | 84791110       |    35734721       | 85218526        |    37059402        |    0.000632    |      0.000910     |
|   70%      | 84791110       |    28608774       | 85218526        |    29617481        |    0.000492    |      0.000728     |
|   80%      | 84791110       |    18619559       | 85218526        |    18463132        |    0.000337    |      0.000547     |
|   90%      | 84791110       |     5569215       | 85218526        |     3689364        |    0.000317    |      0.000397     |
|  minMapQ   |                |                   |                 |                    |                |                   |
|   20       | 84791110       |    41680598       | 85218526        |    43137643        |    0.000733    |      0.001047     |
|   25       | 82304437       |    32023072       | 82728996        |    33802883        |    0.000403    |      0.001018     |
|   30       | 81435963       |    31389492       | 81971309        |    33181390        |    0.000393    |      0.001010     |
|remove_bads |  uniqueOnlys   | only_proper_pair  |                 |                    |                |                   |
|Not applied | 84791110       |    41680598       | 85218526        |    43137643        |    0.000733    |      0.001047     |
|  Applied   | 84791110       |    41680598       | 85218526        |    43137643        |    0.000766    |      0.001052     |
|paralogs    |                |                   |                 |                    |                |                   |
|Not filtered| 84791110       |    41680598       | 85218526        |    43137643        |    0.000733    |      0.001047     |
|filtered    | 84791110       |    41667727       | 85218526        |    43122665        |    0.000735    |      0.001069     |
|No. paralogs|                |       12871       |                 |       14978        |                |                   |

Note: the paralog filtering step is trying to filter out sites where heterozygotes likely comprise more than 50% of all genotypes (likely lumped paralogs)

### Different settings of minInd

- 50% 

<img src="https://hzz0024.github.io/images/DelBay19_fst/Mahattan_ch_ref_fold_minInd50.jpg" alt="img" width="800"/>

- 60%

<img src="https://hzz0024.github.io/images/DelBay19_fst/Mahattan_ch_ref_fold_minInd60.jpg" alt="img" width="800"/>

- 70% 

<img src="https://hzz0024.github.io/images/DelBay19_fst/Mahattan_ch_ref_fold_minInd70.jpg" alt="img" width="800"/>

- 80% 

<img src="https://hzz0024.github.io/images/DelBay19_fst/Mahattan_ch_ref_fold_minInd80.jpg" alt="img" width="800"/>

- 90% 

<img src="https://hzz0024.github.io/images/DelBay19_fst/Mahattan_ch_ref_fold_minInd90.jpg" alt="img" width="800"/>

---
### Different settings of minMapQ

- mapQ25

<img src="https://hzz0024.github.io/images/DelBay19_fst/Mahattan_ch_ref_fold_mapQ25.jpg" alt="img" width="800"/>

- mapQ30

<img src="https://hzz0024.github.io/images/DelBay19_fst/Mahattan_ch_ref_fold_mapQ30.jpg" alt="img" width="800"/>

---
### Effect of -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1

<img src="https://hzz0024.github.io/images/DelBay19_fst/Mahattan_ch_ref_fold_extra.jpg" alt="img" width="800"/>

---
### After removing potential paralogs

<img src="https://hzz0024.github.io/images/DelBay19_fst/Mahattan_ch_ref_fold_paralogs.jpg" alt="img" width="800"/> 

--- 

PART II

The results below are generated using genome-wide chr 1-10 SNP data.

### PCA plots

- challenge group with 50% minInd (after removing inversion regions)
<img src="https://hzz0024.github.io/images/DelBay19_fst/ch_50.pc1-2-1.jpg" alt="img" width="800"/>

- challenge group with 60% minInd (after removing inversion regions)
<img src="https://hzz0024.github.io/images/DelBay19_fst/ch_60.pc1-2-1.jpg" alt="img" width="800"/>

- challenge group with 70% minInd (after removing inversion regions)
<img src="https://hzz0024.github.io/images/DelBay19_fst/ch_70.pc1-2-1.jpg" alt="img" width="800"/>

- challenge group with 80% minInd (after removing inversion regions)
<img src="https://hzz0024.github.io/images/DelBay19_fst/ch_80.pc1-2-1.jpg" alt="img" width="800"/>

- wild group with 50% minInd (after removing inversion regions)
<img src="https://hzz0024.github.io/images/DelBay19_fst/wild_50.pc1-2-1.jpg" alt="img" width="800"/>

- wild group with 60% minInd (after removing inversion regions)
<img src="https://hzz0024.github.io/images/DelBay19_fst/wild_60.pc1-2-1.jpg" alt="img" width="800"/>

- wild group with 70% minInd (after removing inversion regions)
<img src="https://hzz0024.github.io/images/DelBay19_fst/wild_70.pc1-2-1.jpg" alt="img" width="800"/>

- wild group with 80% minInd (after removing inversion regions)
<img src="https://hzz0024.github.io/images/DelBay19_fst/wild_80.pc1-2-1.jpg" alt="img" width="800"/>

### fst distribution and qq-plot

- fst distribution for challenge group with 70% minInd and minMAPQ 25 
<img src="https://hzz0024.github.io/images/qqplot/fst_dis_ch_minInd70_mapq25.jpeg" alt="img" width="800"/>

- fst distribution for challenge group with 50% minInd and minMAPQ 20 
<img src="https://hzz0024.github.io/images/qqplot/fst_dis_ch_minInd50_mapq20.jpeg" alt="img" width="800"/>

- qq-plot for challenge group with 70% minInd and minMAPQ 25 
<img src="https://hzz0024.github.io/images/qqplot/qqplot_ch_minInd70_mapq25.jpeg" alt="img" width="800"/> 

- qq-plot for challenge group with 50% minInd and minMAPQ 20
<img src="https://hzz0024.github.io/images/qqplot/qqplot_ch_minInd50_mapq20.jpeg" alt="img" width="800"/> 

- qqplot for ch_ref_minInd70_mapq25 vs. ch_ref_minInd50_mapq20 
<img src="https://hzz0024.github.io/images/qqplot/50_vs_70.jpg" alt="img" width="800"/>

Conclusion of part I:

1. Weighted mean Fst and number of retained sites goes down linearly as the missingness threshold gets more stringent (from 50% to 90%). One potential explaination for this is that  `<the stringent missingness filter creates a data set enriched for loci suffering from paralogous mapping. This would presumably increase within population diversity and therefore decrease Fst>`. 

2.

