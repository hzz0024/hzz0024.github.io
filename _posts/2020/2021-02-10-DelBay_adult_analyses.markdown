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

- Downsampling 2020 data to account for batch effect

1) Calculate the individual downsampling scale (p) for different relative coverages ratios (0.5x, 0.8x, 1x, 1.2x, and 1.5x relative to 2019 average realized coverage, 2.49). The different ratio settings should be useful for batch effect evaluation - done.            
2) Determine -setMaxDepth for each datasets - done.            
3) Call SNP for each 2020 challenge dataset with different coverage settings - done.            
4) Produce a global SNP list by checking the shared SNPs between 2020 (with different coverage settings) and 2019 datasets - done       
5) Calling SNP again using a global SNP list and a combined 2019 (both challenge and wild, n=331) & 2020 (with different coverage settings - 0.5x, 0.8x, 1x, 1.2x, and 1.5x ) dataset - done    
6) Perform PCA on combined datasets with different coverage settings - done
7) Conclusion: it seems the relative converage below 1.5x could largely eliminate the batch effect. 
8) Check why small number of SNPs could cause batch effect         

- Combined Fisher’s exact tests (after addressing batch effect issue)

1) Update the depth for both Del19 (n=97) and Del20 (n=101) - done      
2) Rerun SNP calling for Del19 and Del20     
3) Generate a global SNP list by extracting the common shared SNPs       
4) Calling allele frequency for each individual population using the global SNP list      
5) Perform Fisher’s exact tests      
6) Examine common shared SNPs      
7) SNP annotation and functional gene exploration       

- SGS test (not restricted by the shared global SNP list) 

1) Determine the global theta value for each population (REF19, REF20, wild?) - REF19 done, the other two under running   
2) Determine -setMaxDepth for each population (CHR19, REF19, CHR20, REF20, HC, NB, SR) - done.      
3) Produce global SNP list for each of the population contrasts (CHR19 vs REF19, CHR20 vs REF20, HC vs. NB, HC vs SR)       
4) Generate allele frequency data for each of the population (with global SNP list, set as -site)       
5) Perform SGS tests on each of the contrast      

- Map a GIS map for DelBay project


### Depth evaluation

Table below summarize the read depth distribution for each dataset. The last column is useful for Angsd -setMaxDepth setting.

| Data       | Mean | Deviation |  SD | Mean+3SD |
|------------|------|-----------|-----|----------|
| DelBay19   | 721  |   41946   | 204 |   1333   |
| DelBay20   | 433  |   33698   | 184 |   985    |
                 
DelBay19: dataset include DelBay19 challenge (n=97) and wild samples (n=234).             
DelBay20: dataset include DelBay20 challenge samples (n=101).                   
       
### Sample coverage

[DelBay19_summary](https://github.com/hzz0024/HG_Code_Bay/blob/master/DelBay_summary/DelBay19_summary_final.xlsx)           
[DelBay20_summary](https://github.com/hzz0024/HG_Code_Bay/blob/master/DelBay_summary/DelBay20_summary_final.xlsx)     

### Downsampling

Tool used for downsampling: [https://gatk.broadinstitute.org/hc/en-us/articles/360037056792-DownsampleSam-Picard-#--RANDOM_SEED](https://gatk.broadinstitute.org/hc/en-us/articles/360037056792-DownsampleSam-Picard-#--RANDOM_SEED)

| Relative Coverage   | Mean | Deviation |  SD | Mean+3SD |
|---------------------|------|-----------|-----|----------|
| 0.5x                | 158  |   1603    | 40  |   278    |
| 0.8x                | 206  |   3933    | 63  |   395    |
| 1.2x                | 271  |   9378    | 97  |   562    |
| 1.5x                | 318  |   14972   | 122 |   685    |
| 1x                  | 239  |   6351    | 80  |   478    |

Note Relative Coverage is determined as the coverage of Del20 samples relative to the Del19 realized coverage: 2.49 (set as 1)

|                                            |  Del19                  | Ratio     |  Del20                                      | Ratio     |
|--------------------------------------------|-------------------------|-----------|---------------------------------------------|-----------|
|                                            |  Del19                  |           |  Del20 (no downsampling)                    |           |
| Total number of sites analyzed             | 512960891               |           | 501510962                                   |           |
| Number of sites retained after filtering   | 2322712                 |           | 5289614                                     |           |
| Private sites                              | 273099                  | 11.76%    | 3240001                                     | 60.81%    |
| Shared sites                               | 2049613                 | 88.24%    | 2049613                                     | 38.47%    |
|                                            |  Del19                  |           |  Del20 (0.5x)                               |           |
| Total number of sites analyzed             | 512960891               |           | 472421980                                   |           |
| Number of sites retained after filtering   | 2322712                 |           | 1370058                                     |           |
| Private sites                              | 1231793                 | 53.03%    | 279139                                      | 20.37%    |
| Shared sites                               | 1090919                 | 46.97%    | 1090919                                     | 79.63%    |
|                                            |  Del19                  |           |  Del20 (0.8x)                               |           |
| Total number of sites analyzed             | 512960891               |           | 484548529                                   |           |
| Number of sites retained after filtering   | 2322712                 |           | 2963811                                     |           |
| Private sites                              | 372380                  | 16.03%    | 1013479                                     | 34.20%    |
| Shared sites                               | 1950332                 | 83.97%    | 1950332                                     | 65.80%    |
|                                            |  Del19                  |           |  Del20 (1x)                                 |           |
| Total number of sites analyzed             | 512960891               |           | 489141910                                   |           |
| Number of sites retained after filtering   | 2322712                 |           | 3600633                                     |           |
| Private sites                              | 290599                  | 12.51%    | 1568520                                     | 43.56%    |
| Shared sites                               | 2032113                 | 87.49%    | 2032113                                     | 56.44%    |
|                                            |  Del19                  |           |  Del20 (1.2x)                               |           |
| Total number of sites analyzed             | 512960891               |           | 492505988                                   |           |
| Number of sites retained after filtering   | 2322712                 |           | 4068394                                     |           |
| Private sites                              | 274375                  | 11.81%    | 2020057                                     | 49.65%    |
| Shared sites                               | 2048337                 | 88.19%    | 2048337                                     | 50.35%    |
|                                            |  Del19                  |           |  Del20 (1.5x)                               |           |
| Total number of sites analyzed             | 512960891               |           | 496183649                                   |           |
| Number of sites retained after filtering   | 2322712                 |           | 4576585                                     |           |
| Private sites                              | 269986                  | 11.62%    | 2523859                                     | 55.15%    |
| Shared sites                               | 2052726                 | 88.38%    | 2052726                                     | 44.85%    |

Next, I am wondering how many identical SNPs among the shared files: Del19_20_global_share_snps.list - no downsampling treatman and Del19_20_1x_share_snps.list - downsampling Del20 for 1x relative coverage. 

```sh
find_shared <- function(f1_name,f2_name) {
  file1 = f1_name
  df1 <- read.delim(file1, header = FALSE, sep='\t')
  df1$V5 <- paste(df1$V1, df1$V2, sep="_")
  message("File ", f1_name, " SNPs: ", length(df1$V5))
  file2 = f2_name
  df2 <- read.delim(file2, header = FALSE, sep='\t')
  df2$V5 <- paste(df2$V1, df2$V2, sep="_")
  message("File ", f1_name, " SNPs: ", length(df2$V5))
  share_snps <- df1[(df1$V5 %in% df2$V5),][,1:4]
  message("Shared SNPs: ", length(share_snps$V1))
  df1_pt <- df1[!(df1$V5 %in% df2$V5),][,1:4]
  message("Private SNPs in ", f1_name, ": ", length(df1_pt$V1)) 
  df2_pt <- df2[!(df2$V5 %in% df1$V5),][,1:4]
  message("Private SNPs in ", f2_name, ": ", length(df2_pt$V1))  
}

find_shared("Del19_20_global_share_snps.list", "Del19_20_1x_share_snps.list")

File Del19_20_global_share_snps.list SNPs: 2049613
File Del19_20_global_share_snps.list SNPs: 2032113
Shared SNPs: 1963240
Private SNPs in Del19_20_global_share_snps.list: 86373
Private SNPs in Del19_20_1x_share_snps.list: 68873
```

The example above shows that simple count of SNP difference between the two shared SNP lists is smaller (2049613-2032113=17500) than the private SNP sites in each file (86373 and 68873), suggesting that downsampling processes may create distinct datasets with non-overlapping SNPs. From the PCA plots, starting from 1.2x I no longer observe the distinct clusters from two sequencing batches. Perhaps checking the non-overlapping SNPs between 1.5x and 1.2x datasets may give me some hints about where the batch effect come from.

- PCA result without downsampling

<img src="https://hzz0024.github.io/images/DelBay_adult/All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_shared_sites_PC1-2.jpg" alt="img" width="800"/>

- PCA result with 1.5x relative coverage

<img src="https://hzz0024.github.io/images/DelBay_adult/All_1p5x_minq20_minmq30_CV30_masked_noinvers_shared_sites_PC1-2.jpg" alt="img" width="800"/>

- PCA result with 1.2x relative coverage

<img src="https://hzz0024.github.io/images/DelBay_adult/All_1p2x_minq20_minmq30_CV30_masked_noinvers_shared_sites_PC1-2.jpg" alt="img" width="800"/>

- PCA result with 1x relative coverage

<img src="https://hzz0024.github.io/images/DelBay_adult/All_1x_minq20_minmq30_CV30_masked_noinvers_shared_sites_PC1-2.jpg" alt="img" width="800"/>

- PCA result with 0.8x relative coverage

<img src="https://hzz0024.github.io/images/DelBay_adult/All_0p8x_minq20_minmq30_CV30_masked_noinvers_shared_sites_PC1-2.jpg" alt="img" width="800"/>

- PCA result with 0.5x relative coverage

<img src="https://hzz0024.github.io/images/DelBay_adult/All_0p5x_minq20_minmq30_CV30_masked_noinvers_shared_sites_PC1-2.jpg" alt="img" width="800"/>


### Global SNP calling using different datasets 

Global SNPs (need edits):

|                                          | Del19                 | Del20                   |
|------------------------------------------|-----------------------|-------------------------|
| Total number of sites analyzed           | 512960891             | 489141910               |
| Number of sites retained after filtering | 2322712               | 3600633                 |
| Private sites in each batch              | 290599 (12.51%)       | 1568520 (43.56%)        |
| Shared sites                             | 2032113 (87.49%)      | 2032113 (56.44%)        |

### Relatedness

<img src="https://hzz0024.github.io/images/DelBay_adult/relatedness_wtoutlier.jpg" alt="img" width="800"/> 

### Diverstiy estimate

<img src="https://hzz0024.github.io/images/DelBay_adult/Diversity.jpg" alt="img" width="800"/>

### MDS

<img src="https://hzz0024.github.io/images/DelBay_adult/All_maf0.05_minq20_minmq25_pctind0.7_CV30_masked_noinvers_shared_sites.jpg" alt="img" width="800"/>

### Combined Fisher's exact test


