---
comments: true
title: uncover the relationship between delta_p and hyposalinity tolerance allele
date: '2020-10-01 12:00'
tags:
  - DelBay19
  - ouliter
  - WGS
  - delta_p
categories:
  - WGS data analysis
---

An trial R script is create to reveal the relationship between delta_p and start p (p0) for hyposalinity tolerance alleles. 

For example, from Fisher’s approach (post [here](https://hzz0024.github.io/2020/09/14/Outlier_detection_Fish_Z-method.html)) I identified 10 potential outliers from the REF-CH vs SR-HC constrasts (fdr < 0.05). By looking at the allele frequency changes I cannot tell which allele is favored by the low salinity, as mafs file only provide the "minor" allele frequencies (here "minor" is determined by the global SNP call at inital stage).

The assumption is that the upbay WILD population may harbors the hyposalinity tolerance allele while downbay population will favor the alternative allele.

Therefore, if I calculate the delta_p for each outlier in the WILD constrast, allele with positive delta_p (p1 = Hope creek, P0 = Shell rock, delta_p = p1-p0) may be the hyposalinity tolerance allele. In constrast, allele with negative delta_p may be the one favored by high salinity. In such case, I will use the alternative allele (i.e. major allele in the maf files) to show the relationship between delta_p and start p (p0) in Challenge vs Reference constrast.

- Wild constrast

| Outlier     | p1(HC)   | p0(SR)   | delta_p1  | Negative |
|-------------|----------|----------|-----------|----------|
| 1_58674134  | 0.109064 | 0.378928 | -0.269864 | Y        |
| 3_24092285  | 0.625825 | 0.232215 | 0.39361   | N        |
| 4_16409729  | 0.282566 | 0.495661 | -0.213095 | Y        |
| 4_18041038  | 0.492074 | 0.408282 | 0.083792  | N        |
| 5_13147035  | 0.19468  | 0.238715 | -0.044035 | Y        |
| 5_16552716  | 0.636006 | 0.46324  | 0.172766  | N        |
| 5_21404457  | 0.544794 | 0.240402 | 0.304392  | N        |
| 5_28997935  | 0.46095  | 0.098022 | 0.362928  | N        |
| 5_52873868  | 0.774964 | 0.427941 | 0.347023  | N        |
| 7_42579646  | 0.242026 | 0.580153 | -0.338127 | Y        |

- Challenge group before identifying the hyposalinity tolerance allele

| Outlier    | p1(CH)   | p0(REF)  | delta_p2  |
|------------|----------|----------|-----------|
| 1_58674134 | 0.090752 | 0.365972 | -0.27522  |
| 3_24092285 | 0.366359 | 0.206749 | 0.15961   |
| 4_16409729 | 0.301443 | 0.721825 | -0.420382 |
| 4_18041038 | 0.612028 | 0.175424 | 0.436604  |
| 5_13147035 | 0.391882 | 0.000003 | 0.391879  |
| 5_16552716 | 0.672417 | 0.240718 | 0.431699  |
| 5_21404457 | 0.643278 | 0.332726 | 0.310552  |
| 5_28997935 | 0.395473 | 0.239822 | 0.155651  |
| 5_52873868 | 0.749064 | 0.490743 | 0.258321  |
| 7_42579646 | 0.110729 | 0.418905 | -0.308176 |

- Challenge group after switching the target allele

| Outlier    | p1(CH)   | p0(REF)  | delta_p3  |
|------------|----------|----------|-----------|
| 1_58674134 | 0.909248 | 0.634028 | 0.27522   |
| 3_24092285 | 0.366359 | 0.206749 | 0.15961   |
| 4_16409729 | 0.698557 | 0.278175 | 0.420382  |
| 4_18041038 | 0.612028 | 0.175424 | 0.436604  |
| 5_13147035 | 0.608118 | 0.999997 | -0.391879 |
| 5_16552716 | 0.672417 | 0.240718 | 0.431699  |
| 5_21404457 | 0.643278 | 0.332726 | 0.310552  |
| 5_28997935 | 0.395473 | 0.239822 | 0.155651  |
| 5_52873868 | 0.749064 | 0.490743 | 0.258321  |
| 7_42579646 | 0.889271 | 0.581095 | 0.308176  |

- Challenge group before identifying the hyposalinity tolerance allele from the SR-HC wild contrasts

<img src="https://hzz0024.github.io/images/outlier/REF_CH_fdr1_p0_p11.jpg" alt="img" width="800"/>

- Challenge group after switching the target allele

<img src="https://hzz0024.github.io/images/outlier/REF_CH_fdr1_p0_p12.jpg" alt="img" width="800"/>

--- 

Above are delta_p plot for 10 outliers identified in REF-CH vs SR-HC contrasts. Now I'd like to look at the modified delta_p patterns for REF-CH vs NB-HC contrasts (16 outliers).

- Challenge group before identifying the hyposalinity tolerance allele from the NB-HC wild contrasts

<img src="https://hzz0024.github.io/images/outlier/REF_CH_fdr5_p0_p1_NB_HC_old.jpg" alt="img" width="800"/>

- Challenge group after switching the target allele

<img src="https://hzz0024.github.io/images/outlier/REF_CH_fdr1_p0_p1_NB_HC_new.jpg" alt="img" width="800"/>



