---
comments: true
title: DelBay19 percentile outlier detection
date: '2020-05-10 12:00'
tags:
  - DelBay19
  - Fst
  - ouliter
  - percentile
  - WGS
categories:
  - WGS data analysis
---

Some references that are helpful in understanding the outlier detections

[Uninformative polymorphisms bias genome scans for signatures of selection](https://bmcevolbiol.biomedcentral.com/articles/10.1186/1471-2148-12-94)

[tutorial on statistical methods for population association studies](nature.com/articles/nrg1916)

[Generating Manhattan Plots in R](https://genome.sph.umich.edu/wiki/Code_Sample:_Generating_Manhattan_Plots_in_R)

[Using Haplotype Information for Conservation Genomics](https://www.sciencedirect.com/science/article/pii/S0169534719303040?via%3Dihub)

### Number of SNPs and Fst values

Basic parameter used during filtering

Base dataset

-minQ 20
-minInd 50% 
-minMapQ 20

Stringent dataset

-minQ 20
-minInd 70% 
-minMapQ 25

Number of SNPs for each saf

| 			 |Total number of sites analyzed | Number of sites retained after filtering |Total number of sites analyzed | Number of sites retained after filtering | 
| -----------|----------|--------------|--------------|--------------|
|   		 |Base      |              | Stringent    |              |
|    HC      | 552206979| 201271276    | 528287677    | 94684245     |
|    ARN     | 547968598| 198211880    | 521758280    | 96654202     |
|    COH     | 548674038| 208674608    | 523144804    | 104599911    |
|    SR      | 550602834| 204514057    | 525334891    | 99762627     |
|    NB      | 546390064| 193778045    | 520958396    | 87044093     |
|    REF     | 549584875| 189886237    | 525678689    | 81957757     |
|    CH      | 545117898| 181283004    | 520698396    | 77282415     |
   
Number of SNPs for each pairwise comparison

| 			     |No. sites | Fst unweight | Fst weight | 
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
