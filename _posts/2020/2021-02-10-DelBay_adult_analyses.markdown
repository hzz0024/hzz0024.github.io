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
| DelBay19   | 739  |   42977   | 207 |   1361   |
| DelBay20   | 440  |   32848   | 181 |   983    |
                 
DelBay19: dataset include DelBay19 challenge (n=103) and wild samples (n=236).             
DelBay20: dataset include DelBay20 challenge samples (n=102).                   
       
### Global SNP calling using different datasets 

Global SNPs:

|                                          | Del20                   | Del19                 |
|------------------------------------------|-------------------------|-----------------------|
| Total number of sites analyzed           | 442020518               | 451557395             |
| Number of sites retained after filtering | 4792076                 | 2021601               |
| Private sites in each batch              | 2957647                 | 187172                |
| Shared sites                             | 1834429                 |                       |

### Diverstiy estimate

<img src="https://hzz0024.github.io/images/DelBay_adult/Diversity.jpg" alt="img" width="800"/>

### PCA 

