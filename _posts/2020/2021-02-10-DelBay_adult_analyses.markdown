---hg111111111žq    
comments: true
title: DelBay adult challenge & wild transects data summary
date: '2021-03-25 12:00'
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

Table 1。 Summary of the read depth distribution for each dataset. The last column is useful for Angsd -setMaxDepth setting.

| Data       | Mean | Deviation |  SD | Mean+3SD |
|------------|------|-----------|-----|----------|
| DelBay19   | 721  |   41946   | 204 |   1333   |
| DelBay20   | 433  |   33698   | 184 |   985    |
                 
DelBay19: dataset include DelBay19 challenge (n=97) and wild samples (n=234).             
DelBay20: dataset include DelBay20 challenge samples (n=101).                   
       
### Sample coverage

[DelBay19_summary](https://github.com/hzz0024/HG_Code_Bay/blob/master/DelBay_summary/DelBay19_summary_final.xlsx)    

Mean realized depth:2.49 (SD = 0.68)      

[DelBay20_summary](https://github.com/hzz0024/HG_Code_Bay/blob/master/DelBay_summary/DelBay20_summary_final.xlsx)  

Mean realized depth:5.17 (SD = 1.40) 

### Downsampling

- Steps

1) Calculate the individual downsampling scale (p) for different relative coverages ratios (0.5x, 0.8x, 1x, 1.2x, and 1.5x relative to 2019 average realized coverage, 2.49). The different ratio settings should be useful for batch effect evaluation - done.            
2) Determine -setMaxDepth for each datasets - done.            
3) Call SNP for each 2020 challenge dataset with different coverage settings - done.            
4) Produce a global SNP list by checking the shared SNPs between 2020 (with different coverage settings) and 2019 datasets - done       
5) Calling SNP again using a global SNP list and a combined 2019 (both challenge and wild, n=331) & 2020 (with different coverage settings - 0.5x, 0.8x, 1x, 1.2x, and 1.5x ) dataset - done    
6) Perform PCA on combined datasets with different coverage settings - done
7) Conclusion: it seems the relative converage below 1.5x could largely eliminate the batch effect. 

Tool used for downsampling: [https://gatk.broadinstitute.org/hc/en-us/articles/360037056792-DownsampleSam-Picard-#--RANDOM_SEED](https://gatk.broadinstitute.org/hc/en-us/articles/360037056792-DownsampleSam-Picard-#--RANDOM_SEED)

Table 2. Summary of SNP depth for Angsd -setMaxDepth setting

| Relative Coverage   | Mean | Deviation |  SD | Mean+3SD |
|---------------------|------|-----------|-----|----------|
| 0.5x                | 158  |   1603    | 40  |   278    |
| 0.6x                | 173  |   2144    | 46  |   312    |
| 0.7x                | 189  |   2931    | 54  |   352    |
| 0.8x                | 206  |   3933    | 63  |   395    |
| 1.2x                | 271  |   9378    | 97  |   562    |
| 1.5x                | 318  |   14972   | 122 |   685    |
| 1x                  | 239  |   6351    | 80  |   478    |

```sh
cat All_1x_maf0.05_minmapq30_minq20_pctind0.7_CV30_masked_mean_sd.log

Mean  Deviation SD
962.438872709503  75772.7262815195  275.268462199213
```

Table 3. Summary of total SNPs analyzed, sites retained after filtering, private sites, and shared sites during downsampling processing. Note Relative Coverage is determined as the coverage of Del20 samples relative to the Del19 (all samples) realized coverage: 2.49 (set as 1)

|                                                 |  Total   SNPs | SNPs after filtering    |  Private | Ratio  | Shared  | Ratio  |
|-------------------------------------------------|---------------|-------------------------|----------|--------|---------|--------|
| Del19 (all sample)                              | 512960891     | 2322712                 | 273099   | 11.76% | 2049613 | 88.24% |
|     Del20 (original)                            | 501510962     | 5289614                 | 3240001  | 60.81% | 2049613 | 38.47% |
| Del19 (all sample)                              | 512960891     | 2322712                 | 1231793  | 53.03% | 1090919 | 46.97% |
|     Del20 (0.5x)                                | 472421980     | 1370058                 | 279139   | 20.37% | 1090919 | 79.63% |
| Del19 (all sample)                              | 512960891     | 2322712                 | 372380   | 16.03% | 1950332 | 83.97% |
|     Del20 (0.8x)                                | 484548529     | 2963811                 | 1013479  | 34.20% | 1950332 | 65.80% |
| Del19 (all sample)                              | 512960891     | 2322712                 | 290599   | 12.51% | 2032113 | 87.49% |
|     Del20 (1x)                                  | 489141910     | 3600633                 | 1568520  | 43.56% | 2032113 | 56.44% |
| Del19 (all sample)                              | 512960891     | 2322712                 | 274375   | 11.81% | 2048337 | 88.19% |
|     Del20 (1.2x)                                | 492505988     | 4068394                 | 2020057  | 49.65% | 2048337 | 50.35% |
| Del19 (all sample)                              | 512960891     | 2322712                 | 269986   | 11.62% | 2052726 | 88.38% |
|     Del20 (1.5x)                                | 496183649     | 4576585                 | 2523859  | 55.15% | 2052726 | 44.85% |

Table 4. Summary of total SNPs analyzed, sites retained after filtering, private sites, and shared sites during downsampling processing, with Del19 dataset constrained to challenge samples only. Note Relative Coverage is determined as the coverage of Del20 samples relative to the Del19 (all samples) realized coverage: 2.49 (set as 1)

|                                                 |  Total   SNPs | SNPs after filtering    |  Private | Ratio  | Shared  | Ratio  |
|-------------------------------------------------|---------------|-------------------------|----------|--------|---------|--------|
| Del19 (only challenge)                          | 482669010     | 2140437                 | 380304   | 17.77% | 1760133 | 82.23% |
|     Del20 (original)                            | 501510962     | 5289614                 | 3529481  | 66.72% | 1760133 | 33.28% |
| Del19 (only challenge)                          | 482669010     | 2140437                 | 1087430  | 50.80% | 1053007 | 49.20% |
|     Del20 (0.5x)                                | 472421980     | 1370058                 | 317051   | 23.14% | 1053007 | 76.86% |
| Del19 (only challenge)                          | 482669010     | 2140437                 | 430595   | 20.12% | 1709842 | 79.88% |
|     Del20 (0.8x)                                | 484548529     | 2963811                 | 1253969  | 42.31% | 1709842 | 57.69% |
| Del19 (only challenge)                          | 482669010     | 2140437                 | 379531   | 17.73% | 1760906 | 82.27% |
|     Del20 (1x)                                  | 489141910     | 3600633                 | 1839727  | 51.09% | 1760906 | 48.91% |
| Del19 (only challenge)                          | 482669010     | 2140437                 | 370648   | 17.32% | 1769789 | 82.68% |
|     Del20 (1.2x)                                | 492505988     | 4068394                 | 2298605  | 56.50% | 1769789 | 43.50% |
| Del19 (only challenge)                          | 482669010     | 2140437                 | 371884   | 17.37% | 1768553 | 82.63% |
|     Del20 (1.5x)                                | 496183649     | 4576585                 | 2808032  | 61.36% | 1768553 | 38.64% |

- Screenshot of different depth in IGV

<img src="https://hzz0024.github.io/images/batch_effect/igv_chr5_16552962.png" alt="img" width="800"/>

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

The example above shows that simple count of SNP difference between the two shared SNP lists is smaller (2049613-2032113=17500) than the private SNP sites in each file (86373 and 68873), suggesting that downsampling processes may create distinct datasets with non-overlapping SNPs. From the PCA plots, starting from 1.2x I no longer observe the divergent clusters from two sequencing batches. Perhaps checking the non-overlapping SNPs between 1.5x and 1.2x datasets may give me some hints about where the batch effect come from, but I will proceed with 1x SNP list for downstream analysis.

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

Table 5. Summary of SNPs used as a global SNP list (2032113 SNPs). Note that private and shared SNP number and the ratio are not the same between two datasets. I am working on other downsampling schemes (0.6 and 0.7x) to balance these numbers.

|                                          | Del19 (all samples)   | Del19 (only challenge)| Del20 (1x)              |
|------------------------------------------|-----------------------|-----------------------|-------------------------|
| Total number of sites analyzed           | 512960891             | 482669010             | 489141910               |
| Number of sites retained after filtering | 2322712               | 2140437               | 3600633                 |
| Private sites in each batch              | 290599 (12.51%)       | 379531 (17.73%)       | 1568520 (43.56%)        |
| Shared sites                             | 2032113 (87.49%)      | 1760906 (82.27%)      | 2032113 (56.44%)        |

### Relatedness

<img src="https://hzz0024.github.io/images/DelBay_adult/relatedness_wtoutlier.jpg" alt="img" width="800"/> 

### Diverstiy estimate

<img src="https://hzz0024.github.io/images/DelBay_adult/Diversity.jpg" alt="img" width="800"/>

### MDS

<img src="https://hzz0024.github.io/images/DelBay_adult/All_maf0.05_minq20_minmq25_pctind0.7_CV30_masked_noinvers_shared_sites.jpg" alt="img" width="800"/>

### Combined Fisher's exact test

1) Update the depth for both Del19 (n=97) and Del20 (n=101) - done      
2) Rerun SNP calling for Del19 and Del20 - done    
3) Generate a global SNP list by extracting the common shared SNPs - done    
4) Calling allele frequency for each individual population using the global SNP list - done    
5) Perform Fisher’s exact tests - done
6) Examine common shared SNPs - done     
7) SNP annotation and functional gene exploration  

In the p-value combination step, because the weighted Z-test (also called ‘Stouffer’s method) is more power and more precision than does Fisher’s test ([Whitlock 2005](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1420-9101.2005.00917.x)), I used Z method to combine the p-values from the parallel tests. It favours symmetric rejection and is less sensitive to a single low p-value, requiring more consistently low p-values to yield a low combined p-value. 

The combinePValues function in the [scran R package](https://rdrr.io/bioc/scran/man/combinePValues.html) was used to perform Z method.

Results:

| Contrasts                  | No. outliers |
|----------------------------|--------------|
| REF19_CHR19_NB_HC          | 13           |
| REF19_CHR19_SR_HC          | 4            |
| Shared                     | 1            |
| REF19_SR_ARN_COH (control) | 0            |
| REF20_CHR20_NB_HC          | 7            |
| REF20_CHR20_SR_HC          | 0            |
| Shared                     | 0            |
| REF20_SR_ARN_COH (control) | 3            |
| REF19_REF20_CHR19_CHR20    | 696          |

The detailed SNP lists are [here](https://docs.google.com/spreadsheets/d/1hDH_lp_BQC9grGAK2vbZZBW-tzpiZCvMvajyMe-rYUo/edit?usp=sharing)

Now take a look at the allele frequency changes at these potential outliers.

- REF19_CHR19_NB_HC (13 outliers)

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_REF19_CHR19_NB_HC.jpg" alt="img" width="800"/>

One of them, SNP NC_035784.1_16552716 is shared between REF19_CHR19_NB_HC and REF19_CHR19_SR_HC. This SNPs is located in the gene Actin-depolymerizing factor 1-like (LOC111134891). Below is the delta_p patterns for this SNP

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_NC_035784.1_16552716.jpg" alt="img" width="800"/>

- REF20_CHR20_NB_HC (7 outliers)

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_REF20_CHR20_NB_HC.jpg" alt="img" width="800"/>

- REF19_REF20_CHR19_CHR20 (696 outlies)

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_control.jpg" alt="img" width="800"/>

Zoom in on some SNPs at chromosome 5

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_control_chr5.jpg" alt="img" width="800"/>

Zoom in on some SNPs at chromosome 1

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_control_chr1.jpg" alt="img" width="800"/>

### Single-generation selection (SGS) 

1) Determine the global theta value for each population (REF19, REF20, wild?) - REF19 done, the other two under running   
2) Determine -setMaxDepth for each population (CHR19, REF19, CHR20, REF20, HC, NB, SR) - done.      
3) Produce global SNP list for each of the population contrasts (CHR19 vs REF19, CHR20 vs REF20, HC vs. NB, HC vs SR) - done      
4) Generate allele frequency data for each of the population (with global SNP list, set as -site) - done      
5) Perform SGS tests on each of the contrast 

Table 6. Summary of SNP depth for Angsd -setMaxDepth settings. Note CHR19-REF19, CHR20-REF20, HC-SR, and HC-NB are contrasts used for SGS tests.

|     Contrasts       | Mean | Deviation |  SD | Mean+3SD |
|---------------------|------|-----------|-----|----------|
| CHR19-REF19         | 189  |   2966    | 54  |   352    |
| CHR20-REF20         | 433  |   33698   | 184 |   983    |
| HC-SR               | 216  |   4487    | 67  |   417    |
| HC_NB               | 210  |   3994    | 63  |   400    |

Note: p-value distribution is initially odd, probably due the usage of quantile-based p-value. After consulting with CSCU, now I switched to 2-sides p-value estimation with naive exterme delta-p counts.

p-value distribution for CHR19-REF19 contrast using 2-side p-value

<img src="https://hzz0024.github.io/images/DelBay_adult/ps_2side.jpeg" alt="img" width="800"/>

Table 6. Number of outliers identifed from SGS test for each population contrast (FDR < 0.05).

|     Contrasts       | Outliers |
|---------------------|----------|
| CHR19-REF19         | 3265     |
| CHR20-REF20         | 180      |              
| HC-SR               | 2444     |                     
| HC-NB               | 2559     |              

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_SGS.jpg" alt="img" width="800"/>

Shared SNPs:

CHR20-REF20 & HC_NB: NC_035782.1_62395550

CHR19-REF19 & HC_NB or HC_SR: 12 SNPS

CHR19-REF19 & HC_NB & HC_SR: NC_035786.1_8772371

### Probabilistic Random Forest

Following the paper by Reis et al. 2018. [Probabilistic Random Forest: A machine learning algorithm for noisy datasets](https://arxiv.org/pdf/1811.05994.pdf). I am trying to incorporate the genotype likelihood into random forest test.

The initial trial is performed on 3006 outliers SNPs identified from SGS CHR19-REF19 contrasts (total 2032113 SNPS).

<img src="https://hzz0024.github.io/images/DelBay_adult/PRF_accuracy.jpg" alt="img" width="800"/>

- shared SNPs among repeat runs (1000, 500, 100, 50 are top SNPs based on importance)

1000
NC_035780.1_50001008 NC_035782.1_45223277 NC_035784.1_10248759 *NC_035784.1_16552962* NC_035784.1_18676441 NC_035784.1_83138159 NC_035786.1_33698717

500
NC_035780.1_50001008 NC_035782.1_45223277 *NC_035784.1_16552962* NC_035784.1_18676441 NC_035786.1_33698717

100
*NC_035784.1_16552962*

50
*NC_035784.1_16552962*

- Zoom-in for NC_035784.1_16552962 (red) and NC_035784.1_16552716 (lightgreen), which both located in Actin-depolymerizing factor 1-like (LOC111134891).

<img src="https://hzz0024.github.io/images/DelBay_adult/Mahattan_PRF.jpg" alt="img" width="800"/>

### Genotype-environment association


