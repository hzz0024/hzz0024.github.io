---
comments: true
title: DelBay adult PCA & MDS summary
date: '2021-04-01 12:00'
tags:
  - DelBay
  - Challenge
  - Wild 
  - WGS
  - QC
categories:
  - WGS data analysis
--- 


Table 1. Summary of total SNPs analyzed, sites retained after filtering, private sites, and shared sites during downsampling processing. Note Relative Coverage is determined as the coverage of Del20 samples relative to the Del19 (all samples) realized coverage: 2.49 (set as 1)

|                                                 |  Total   SNPs | SNPs after filtering    |  Private | Ratio  | Shared  | Ratio  |
|-------------------------------------------------|---------------|-------------------------|----------|--------|---------|--------|
| Del19 (all sample)                              | 512960891     | 2322712                 | 273099   | 11.76% | 2049613 | 88.24% |
|     Del20 (no downsampling)                     | 501510962     | 5289614                 | 3240001  | 60.81% | 2049613 | 38.47% |
| Del19 (all sample)                              | 512960891     | 2322712                 | 511985   | 22.04% | 1810727 | 77.96% |
|     Del20 (0.7x)                                | 481289979     | 2519717                 | 708990   | 28.14% | 1810727 | 71.86% |
| Del19 (all sample)                              | 512960891     | 2322712                 | 290599   | 12.51% | 2032113 | 83.97% |
|     Del20 (1x)                                  | 489141910     | 3600633                 | 1568520  | 43.56% | 2032113 | 56.44% |

Table 2. Summary of total SNPs analyzed, sites retained after filtering, private sites, and shared sites during downsampling processing, with Del19 dataset constrained to challenge samples only. Note Relative Coverage is determined as the coverage of Del20 samples relative to the Del19 (all samples) realized coverage: 2.49 (set as 1)

|                                                 |  Total   SNPs | SNPs after filtering    |  Private | Ratio  | Shared  | Ratio  |
|-------------------------------------------------|---------------|-------------------------|----------|--------|---------|--------|
| Del19 (only challenge)                          | 482669010     | 2140437                 | 380304   | 17.77% | 1760133 | 82.23% |
|     Del20 (no downsampling)                     | 501510962     | 5289614                 | 3529481  | 66.72% | 1760133 | 33.28% |
| Del19 (only challenge)                          | 482669010     | 2140437                 | 522567   | 24.41% | 1617870 | 75.59% |
|     Del20 (0.7x)                                | 481289979     | 2519717                 | 901847   | 35.79% | 1617870 | 64.21% |
| Del19 (only challenge)                          | 482669010     | 2140437                 | 379531   | 17.73% | 1760906 | 82.27% |
|     Del20 (1x)                                  | 489141910     | 3600633                 | 1839727  | 51.09% | 1760906 | 48.91% |

### PCA & MDS plots for each datasets

PCAngsd-PCA for shared SNPs (no downsampling, 2049613 SNPs)

<img src="https://hzz0024.github.io/images/DelBay_adult/All_maf0.05_minq20_minmq30_pctind0.7_CV30_masked_noinvers_shared_sites_PC1-2.jpg" alt="img" width="800"/>

Angsd-PCA for shared SNPs (no downsampling, 2049613 SNPs)

<img src="https://hzz0024.github.io/images/DelBay_adult/angsd_pca_no_downsampling.jpeg" alt="img" width="800"/>

Angsd-MDS for shared SNPs (no downsampling, 2049613 SNPs)

<img src="https://hzz0024.github.io/images/DelBay_adult/angsd_MDS_no_downsampling.jpeg" alt="img" width="800"/>

PCAngsd-PCA for shared SNPs (downsampling to 1x, 2032113 SNPs)

<img src="https://hzz0024.github.io/images/DelBay_adult/All_1x_minq20_minmq30_CV30_masked_noinvers_shared_sites_PC1-2.jpg" alt="img" width="800"/>

Angsd-PCA for shared SNPs (downsampling to 1x, 2032113 SNPs)

<img src="https://hzz0024.github.io/images/DelBay_adult/angsd_pca_1x.jpeg" alt="img" width="800"/>

Angsd-MDS for shared SNPs (downsampling to 1x, 2032113 SNPs)

<img src="https://hzz0024.github.io/images/DelBay_adult/angsd_mds_1x.jpeg" alt="img" width="800"/>

PCAngsd-PCA for shared SNPs (downsampling to 0.7x, 1810727 SNPs)

<img src="https://hzz0024.github.io/images/DelBay_adult/All_0p7x_minq20_minmq30_CV30_masked_noinvers_shared_sites_PC1-2.jpg" alt="img" width="800"/>

Angsd-PCA for shared SNPs (downsampling to 0.7x, 1810727 SNPs)

<img src="https://hzz0024.github.io/images/DelBay_adult/angsd_pca_0.7x.jpeg" alt="img" width="800"/>

Angsd-MDS for shared SNPs (downsampling to 0.7x, 1810727 SNPs)

<img src="https://hzz0024.github.io/images/DelBay_adult/angsd_MDS_0.7x.jpeg" alt="img" width="800"/>

### Conclusion for PCA

1. PCAngsd is more susceptible to batch effect caused by coverage difference than ANGSD.              
2. Angsd-PCA and Angsd-MDS found no batch effects in 1) no downsampling; 2) downsampling to 1x relative coverage, or 3) 0.7x relative coverage datasets.       
3. Homogenizing genetic differentiation was found among the populations.                    
4. It would be interesting to see how option of -doIBS 2 (Concensus base) change the Angsd-PCA and Angsd-MDS results (under running).       


### Global pairwise Fst

Question: Del19 and Del20 challenge group is souced from the same locale but at different years. Why did I identify so many outliers between two years?

To answer this question, I want to take a look at the pairwise global Fst among population contrasts.

Table 3 Pairwise global Fst among population contrasts using 2032113 SNPs (shared SNPs from 1x downsampling dataset). The weighted Fst for a region is the ratio between the sum of As and the sum of B. The unweighted is the mean of the persite ratios. A is the alpha from the reynolds 1983 (or Bhatia) and B is the alpha + beta.

| Pop_pair    |    Fst.Unweight     | Fst.Weight |
|-------------|---------------------|------------|
| Sur19_Ref19 | 2.42E-04            | 3.52E-04   |
| Sur20_Ref20 | 1.87E-04            | 2.72E-04   |
| HC_NB       | 3.67E-04            | 5.10E-04   |
| HC_SR       | 2.33E-04            | 3.19E-04   |
| ARN_COH     | 6.30E-05            | 1.85E-04   |
| Ref19_Ref20 | 1.00E-03            | 1.50E-03   |
| Sur19_Sur20 | 1.98E-03            | 2.62E-03   |

The global Fst among wild contrasts make sense. The smaller Fst in Sur20-Ref20 relative to Sur19-Ref19 is surprising to me, as we only used non-survivals as the reference, whereas both survivals and non-survivals were grouped as 2019 reference samples.

The Fst values between Sur19-Sur20 and Ref19-Ref20 show elevated Fst due to Temporal genetic structure.  



