---
comments: true
title: DelBay adult challenge & wild transects data summary
date: '2021-02-19 12:00'
tags:
  - DelBay
  - Challenge
  - Wild 
  - WGS
  - QC
categories:
  - WGS data analysis
--- 

### Depth evaluation

Table below summarize the read depth distribution for each dataset. The last column is useful for Angsd -setMaxDepth setting.

| Data       | Mean | Deviation |  SD | Mean+3SD |
|------------|------|-----------|-----|----------|
| DelBay19   | 725  |   43615   | 209 |   1351   |
| DelBay20   | 437  |   34197   | 185 |   991    |
                 
DelBay19: dataset include DelBay19 challenge (n=103) and wild samples (n=236).             
DelBay20: dataset include DelBay20 challenge samples (n=102).                   
       
### Sample coverage



### Global SNP calling using different datasets 

Global SNPs:

|                                          | Del20                   | Del19                 |
|------------------------------------------|-------------------------|-----------------------|
| Total number of sites analyzed           | 503801179               | 514216357             |
| Number of sites retained after filtering | 5675253                 | 2505019               |
| Private sites in each batch              | 3400996 (59.93%)        | 230762 (9.21%)        |
| Shared sites                             | 2274257 (40.07%)        | 2274257 (90.79%)      |

### Diverstiy estimate

<img src="https://hzz0024.github.io/images/DelBay_adult/Diversity.jpg" alt="img" width="800"/>

### MDS

<img src="https://hzz0024.github.io/images/DelBay_adult/All_maf0.05_minq20_minmq25_pctind0.7_CV30_masked_noinvers_shared_sites.jpg" alt="img" width="800"/>

### PCA 

<img src="https://hzz0024.github.io/images/DelBay_adult/All_maf0.05_minq20_minmq25_pctind0.7_CV30_masked_noinvers_shared_sites_PC1-2.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/DelBay_adult/All_maf0.05_minq20_minmq25_pctind0.7_CV30_masked_noinvers_shared_sites_PC2-3.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/DelBay_adult/All_maf0.05_minq20_minmq25_pctind0.7_CV30_masked_noinvers_shared_sites_PC3-4.jpg" alt="img" width="800"/>

### Combined Fisher's exact test

