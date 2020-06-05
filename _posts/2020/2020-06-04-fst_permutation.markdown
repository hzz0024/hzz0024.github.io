---
comments: true
title: DelBay19 fst outlier detection using permutation method
date: '2020-06-04 12:00'
tags:
  - DelBay19
  - Fst
  - ouliter
  - WGS
categories:
  - WGS data analysis
---

In this post I tried to perform the significance test for Fst outliers. This could be done by randomly drawing individauls from two populations and obtain a neutral dataset of Fst (or permutation test), and using that dataset for p-value calculation. The workframe is shown below,

1) For each SNP k: have the observed Fst_k

2) Repeated scrambled the subpopulation membership to get Fst_k_hat, say M = 100 times. For each Fst_k, we have 100 resampled Fst_k_hat values - this will have a mean and SE

3) Or p-value = proportion of resamples that are larger than Fst_k - WINNER! i.e. the proportion of resamples where the Ha is true. Assume that p-value =0.01 -> that means that only 1% of resamples are greater than Fst_k

4) we could run a one sample t-test (could be a one-side test or two-sided)    
   Ho: mean resampled Fst >= Fst_k    
   Ha: mean resampled Fst < Fst_k   

5) Apply FDR to p-values

### Data generation

-minQ 20
-minInd 70% 
-minMapQ 25

Step 1 Using all samples of two contrast groups (ch and ref here) for maf file generation. The maf file will be used for minor allele frequency filter. Here I used three filters, 0.05, 0.10, and 0.20 to generate the SNP list for angsd running.

```sh

files=$(ls *.mafs.gz)
for file in $files; do
    file=${file/_doMAF.mafs.gz/}
    zcat $file'_doMAF.mafs.gz' | tail -n +2 | awk '$6>0.05 {print ;}' | awk '{print $1,$2,$3,$4}' > $file'_maf05.snplist'
    zcat $file'_doMAF.mafs.gz' | tail -n +2 | awk '$6>0.10 {print ;}' | awk '{print $1,$2,$3,$4}' > $file'_maf10.snplist'
    zcat $file'_doMAF.mafs.gz' | tail -n +2 | awk '$6>0.20 {print ;}' | awk '{print $1,$2,$3,$4}' > $file'_maf20.snplist'
done

walltime used =  17648.00 sec
```

| 	 MAF	 | No. SNPs |  
| -----------|----------|
|  no filter | 34313641 | 
|    0.05    | 872160   | 
|    0.10    | 539901   |
|    0.20    | 291145   |


Step 2 Perform the angsd running and generate the saf file for downstream Fst calling. The SNP list generated from step 1 will be used here.

```

Number of SNPs for each saf

| 			 |Total number of sites analyzed | Number of sites retained after filtering |  
| -----------|----------|--------------|
|    REF     | 549584875| 189886237    |
|    CH      | 545117898| 181283004    |
   
Number of SNPs for each pairwise comparison

| 			 |No. sites | Fst unweight | Fst weight | 
| -----------|----------|--------------|------------|
|    CH_REF  | 170859421|  0.000664    |   0.000883 |
|    HC_ARN  | 185953051|  0.001151    |   0.000852 |
|    HC_COH  |          |              |            |
|    HC_SR   |          |              |            |
|    HC_NB   | 183396017|  0.000490    |   0.000827 |
|    ARN_COH | 189084056|  0.000815    |   0.000675 |
|    ARN_SR  |          |              |            |
|    ARN_NB  | 182211161|  0.001035    |   0.000827 |
|    COH_SR  |          |              |            |
|    COH_NB  |          |              |            |
|    SR_NB   |          |              |            |

---
### Fst plots

- CH_REF

Mahattan plot based on single snp (Fst ranges from 0 - 0.2)

<img src="https://hzz0024.github.io/images/DelBay19_fst/Mahattan_ch_ref_singlesnp_fold.jpg" alt="img" width="800"/> 

- HC_ARN 

Mahattan plot based on single snp (Fst ranges from 0 - 0.2)

<img src="https://hzz0024.github.io/images/DelBay19_fst/Mahattan_HC_ARN_singlesnp_fold.jpg" alt="img" width="800"/>

- ARN_COH

Mahattan plot based on single snp (Fst ranges from 0 - 0.2)

<img src="https://hzz0024.github.io/images/DelBay19_fst/Mahattan_ARN_COH_singlesnp_fold.jpg" alt="img" width="800"/>

- ARN_NB

Mahattan plot based on single snp (Fst ranges from 0 - 0.2)

<img src="https://hzz0024.github.io/images/DelBay19_fst/Mahattan_ARN_NB_singlesnp_fold.jpg" alt="img" width="800"/>

---

### Outlier detection (Percentile)

Python script for percentile calculation: percentile.py        
Python script for common shared loci detection: common.py

Salinity from low to high:   
Hope Creek (HC)        
Arnolds (ARN)   
Cohansey (COH)   
Shell Rock (SR)   
New Beds (NB)   

<img src="https://hzz0024.github.io/images/map1.jpg" alt="img" width="400"/>

- 99.9% Percentile Results 

| Group	     |Number of SNP bins (1kb)| Percentile   |  Fst threshold | Max Fst    |No. outlier|
| -----------|------------------------|--------------|----------------|------------|-----------|
|   CH_REF   |        277410          |     99.9     |    0.0217      |   0.0802   |   278     |
|   HC_NB    |        289111          |     99.9     |    0.0220      |   0.0839   |   290     |
|   HC_ARN   |        290807          |     99.9     |    0.0212      |   0.0674   |   291     |
|   ARN_COH  |        294114          |     99.9     |    0.0222      |   0.0546   |   295     |
|   ARN_NB   |        287305          |     99.9     |    0.0217      |   0.0912   |   288     |

- 99% Percentile Results

| Group	     |Number of SNP bins (1kb)| Percentile   |  Fst threshold | Max Fst    |No. outlier|
| -----------|------------------------|--------------|----------------|------------|-----------|
|   CH_REF   |        277410          |     99       |    0.0117      |   0.0802   |  2775     |
|   HC_NB    |        289111          |     99       |    0.0119      |   0.0839   |  2892     |
|   HC_ARN   |        290807          |     99       |    0.0117      |   0.0674   |  2908     |
|   ARN_COH  |        294114          |     99       |    0.0118      |   0.0546   |  2942     |
|   ARN_NB   |        287305          |     99       |    0.0114      |   0.0912   |  2874     |

- Comparing challenge vs. wild 

| Groups used for comparison     | No. shared outlier (99.9 percentile)|No. shared outlier (99 percentile)|
| -------------------------------|-------------------------------------|----------------------------------|
| CH_REF & HC_NB                 |    4                                |        65                        |
| CH_REF & HC_ARN                |    0                                |        54                        |
| CH_REF & ARN_COH               |    0                                |        63                        |
| CH_REF & ARN_NB                |    0                                |        60                        |

<img src="https://hzz0024.github.io/images/comp99.jpg" alt="img" width="800"/>

Most of the outliers are not in the inversion regions (NC_035784.1 60600000-80200000)

For example, in CH_REF & HC_NB comparison

5_8234000   
5_12612000   
5_12699000   
5_13141000   
5_16297000   
5_16551000   
5_18934000   
5_19565000   
5_31984000   
5_43506000   
__5_62134000__  
__5_64129000__   
__5_65798000__   
__5_65800000__   
__5_65822000__   
__5_70279000__   
__5_70908000__   
5_86165000   
5_86240000   
5_90503000   
