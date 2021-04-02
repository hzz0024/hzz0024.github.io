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



