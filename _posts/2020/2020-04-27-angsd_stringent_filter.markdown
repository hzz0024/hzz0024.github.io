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

The control run here used -minQ 20 -minMapQ 20 -minInd 25 or 24 (50% missing rate).

### Different settings of minInd

|     	     |Total sites (ch)|Sites retained (ch)|Total sites (ref)|Sites retained (ref)|   
| -----------|----------------|-------------------|-----------------|--------------------|
|   50%      | 84791110       |    41680598       | 85218526        |    43137643        |
|   60%      | 84791110       |    35734721       | 85218526        |    37059402        |
|   70%      | 84791110       |    28608774       | 85218526        |    29617481        |
|   80%      | 84791110       |    18619559       | 85218526        |    18463132        |
|   90%      | 84791110       |     5569215       | 85218526        |     3689364        |

- 50% 
- 60% 
- 70% 
- 80% 
- 90% 
   
<img src="https://hzz0024.github.io/images/ngsLD/QC_DelBay19_12_maxD2000.jpg" alt="img" width="800"/>

---
### Different settings of minMapQ

---
### Effect of -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1

---
### Remove potential paralogs




