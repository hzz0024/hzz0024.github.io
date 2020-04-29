---
comments: true
title: DelBay19_stringent_filter
date: '2020-04-27 12:00'
tags:
  - DelBay19
  - angsd
  - filter
  - WGS
categories:
  - WGS data analysis
---

In this post I'd like to evaluate how parameters such as minInd and minMapQ affect the SNP number produced from the low-coverage data. 

The control run here used -minQ 20 -minMapQ 20 -minInd 25 (ch) or 24 (ref) - corresponding 50% missing rate.

### Set the MinDepth and MaxDepth based on global depth distribution

<img src="https://hzz0024.github.io/images/ngsLD/QC_DelBay19_12_maxD2000.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/ngsLD/QC_CVreseq_12_maxD2000.jpg" alt="img" width="800"/>

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
### Remove potential paralogs




