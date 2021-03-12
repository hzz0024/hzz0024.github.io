---
comments: true
title: DelBay adult challenge & wild transects data summary
date: '2021-03-12 12:00'
tags:
  - DelBay
  - Challenge
  - Wild 
  - WGS
  - QC
categories:
  - WGS data analysis
--- 

`Downsampling 2020 data to account for batch effect`     
1) Calculate the individual downsampling scale (p) for different relative coverages ratios (0.5x, 0.8x, 1x, 1.2x, and 1.5x relative to 2019 average realized coverage, 2.49). The different ratio settings should be useful for batch effect evaluation - done.            
2) Determine -setMaxDepth for each datasets - done.            
3) Call SNP for each 2020 challenge dataset with different coverage settings - under running.            
4) Produce a global SNP list by checking the shared SNPs between 2020 (with different coverage settings) and 2019 datasets.           
5) Calling SNP again using a global SNP list and a combined 2019 (both challenge and wild, n=331) & 2020 (with different coverage settings - 0.5x, 0.8x, 1x, 1.2x, and 1.5x ) dataset.          
6) Perform PCA on combined datasets with different coverage settings.          

`Combined Fisher’s exact tests (after addressing batch effect issue)`    
1) Update the depth for both Del19 (n=97) and Del20 (n=101) - done      
2) Rerun SNP calling for Del19 and Del20     
3) Generate a global SNP list by extracting the common shared SNPs       
4) Calling allele frequency for each individual population using the global SNP list      
5) Perform Fisher’s exact tests      
6) Examine common shared SNPs      
7) SNP annotation and functional gene exploration       

`SGS test (not restricted by the shared global SNP list)`     
1) Determine the global theta value for each population (REF19, REF20, wild?) - need discuss with Matt      
2) Determine -setMaxDepth for each population (CHR19, REF19, CHR20, REF20, HC, NB, SR) - done.      
3) Produce global SNP list for each of the population contrasts (CHR19 vs REF19, CHR20 vs REF20, HC vs. NB, HC vs SR)       
4) Generate allele frequency data for each of the population (with global SNP list, set as -site)      
5) Perform SGS tests on each of the contrast      

### Depth evaluation

Table below summarize the read depth distribution for each dataset. The last column is useful for Angsd -setMaxDepth setting.

| Data       | Mean | Deviation |  SD | Mean+3SD |
|------------|------|-----------|-----|----------|
| DelBay19   | 721  |   41946   | 204 |   1333   |
| DelBay20   | 433  |   33698   | 184 |   985    |
                 
DelBay19: dataset include DelBay19 challenge (n=97) and wild samples (n=234).             
DelBay20: dataset include DelBay20 challenge samples (n=101).                   
       
### Sample coverage


### Downsampling

| Relative Coverage   | Mean | Deviation |  SD | Mean+3SD |
|---------------------|------|-----------|-----|----------|
| 0.5x                | 158  |   1603    | 40  |   278    |
| 0.8x                | 206  |   3933    | 63  |   395    |
| 1.2x                | 271  |   9378    | 97  |   562    |
| 1.5x                | 318  |   14972   | 122 |   685    |
| 1x                  | 239  |   6351    | 80  |   478    |

Note Relative Coverage is determined as the coverage of Del20 samples relative to the Del19 realized coverage: 2.49 (set as 1)


### Global SNP calling using different datasets 

Global SNPs (need edits):

|                                          | Del20                   | Del19                 |
|------------------------------------------|-------------------------|-----------------------|
| Total number of sites analyzed           | 501510962               | 514216357             |
| Number of sites retained after filtering | 5328491                 | 2505019               |
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


