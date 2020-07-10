---
comments: true
title: DelBay19 ngsLD 
date: '2020-07-09 12:00'
tags:
  - DelBay19
  - deltap
  - ouliter
  - ngsLD
  - WGS
categories:
  - WGS data analysis
---

This post shows the results from ngsLD. The detailed steps for data creation is in post [Rebuild DelBay19 data](https://hzz0024.github.io/2020/07/08/DelBay_data_redo.html)

Briefly, I created ngsLD output for each population: CH, REF, HC, ARN, COH, SR, and NB. For each population, the LD decay curvy will be ploted for chromosome 1-10. 

Initially, I evaluated what is the optimal setting of --max_kb_dist. This means the maximum distance between SNPs (in Kb) to calculate LD. I created three input with 1, 5, and 10 kb as the mximum distance between SNPs for LD calculation.

- REF 1kb

<img src="https://hzz0024.github.io/images/ngsLD/REF_1k.jpg" alt="img" width="800"/>

- REF 5kb

<img src="https://hzz0024.github.io/images/ngsLD/REF_5k.jpg" alt="img" width="800"/>

- REF 10kb

<img src="https://hzz0024.github.io/images/ngsLD/REF_10k.jpg" alt="img" width="800"/>

It looks the curvy became flat after 5kb. I then create the LD decay curvy for each of the examined population with --max_kb_dist 5

- CH

<img src="https://hzz0024.github.io/images/ngsLD/CH_5k.jpg" alt="img" width="800"/>

- REF

<img src="https://hzz0024.github.io/images/ngsLD/REF_5k.jpg" alt="img" width="800"/>

- HC

<img src="https://hzz0024.github.io/images/ngsLD/HC_5k.jpg" alt="img" width="800"/>

- ARN

<img src="https://hzz0024.github.io/images/ngsLD/ARN_5k.jpg" alt="img" width="800"/>

- COH

<img src="https://hzz0024.github.io/images/ngsLD/COH_5k.jpg" alt="img" width="800"/>

- SR

<img src="https://hzz0024.github.io/images/ngsLD/SR_5k.jpg" alt="img" width="800"/>

- NB

<img src="https://hzz0024.github.io/images/ngsLD/NB_5k.jpg" alt="img" width="800"/>

