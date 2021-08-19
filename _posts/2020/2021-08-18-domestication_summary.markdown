---
comments: true
title: 600K Domestication Project summary
date: '2021-08-18 12:00'
tags:
  - 600K
  - SNP array
  - oyster
  - domestication 
  - WGS
  - summary
categories:
  - WGS data analysis
--- 

### 0. Overall population information

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/pop_info.jpeg" alt="img" width="800"/>

UMFS are selected groups from LIW+MEW+NEH, NYH is selected group from LIW, NEHs are selected groups from LIW, DBX1 is selected groups from DBW F1, DBX2 is highly inbred group many generations from DBW, DBX3 is from mixture of DBW and NEH. VIMS families are selected groups from DBW, Golf oysters and CBW. UNC groups are selected from NCW. The groups whose has “W” at the last in their names should be wild group.

### 1. VCF process

Filtering parameters: --maf 0.05, --max-missing 0.7, --chr 1-10

Number of SNPs in the original vcf: 300,446  
Number of total samples: 842    
Number of SNPs after filtering: 159,849      
Number of total samples: 514

Figure 1. SNP summary for vcf file. From top to bottom: Proportion of missing data per SNP; Proportion of missing data per individual; distribution of minor allele frequencies

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/vcf_summary.jpg" alt="img" width="800"/>

Our vcf data has a very promising profile for downstream analyses - 1) clearly most individuals have a call at almost every site; 2) the proportion of missing data per individual is low; 3) it is clear that a large number of variants have low frequency alleles.

### 2. Population divergence based on global *F*<sub>ST</sub>

Figure 2. Window-based *F*<sub>ST</sub> for each wild vs. selected line. *F*<sub>ST</sub> was calculated with 150 SNPs per non-overlapping windows. Solid Line represents the median value.

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/Fst_violin_plot.jpg" alt="img" width="800"/>

### 3. Diverstiy estimate

Figure 3.1 Heterozygosity comparison between wild (left four) and selected (right four) lines. Heterozygosity was estimated for each individual. Solid Line represents the median value.

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/Het_ind_plot.jpg" alt="img" width="800"/>

Figure 3.2 Nucleotide diversity (*pi*) comparison between wild (left four) and selected (right four) lines. Pi was estimated with 150 SNPs per non-overlapping windows using Simon martin's script (https://github.com/simonhmartin/genomics_general). Solid Line represents the median value.

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/Pi_150snp_per_window_plot.jpg" alt="img" width="800"/>

Figure 3.3 Nucleotide diversity (*pi*) comparison between wild (left four) and selected (right four) lines. Pi was estimated with 500K bp per non-overlapping windows using VCFtools. Solid Line represents the median value.

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/pi_vcftools_ind_plot.jpg" alt="img" width="800"/>

Note Figure 3.2 and 3.3 both measure the *pi* but their absolute values differences is over an order of magnitude. This is because VCFtools (results in Figure 3.2) underestimates the pi by dividing the total length of the window (that is 500K), whereas Simon martin's script (https://github.com/simonhmartin/genomics_general) only takes vcf input sites into account (see [https://github.com/simonhmartin/genomics_general/issues/22](https://github.com/simonhmartin/genomics_general/issues/22)). A third approach like Angsd is needed to calculate the *pi* .

### 4. Inbreeding coefficient

Figure 4 Inbreeding coefficient comparison between wild (left four) and selected (right four) lines. Inbreeding coefficient was estimated for each individual using the same amount of SNPs. Solid Line represents the median value.

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/Fis_ind_plot.jpg" alt="img" width="800"/>

### 5. PCA 

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/PCA2.jpeg" alt="img" width="800"/>

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/PCA3.jpeg" alt="img" width="800"/>

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/PCA6.jpeg" alt="img" width="800"/>

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/PCA5.jpeg" alt="img" width="800"/>

### 6. Pairwise *F*<sub>ST</sub> contrasts

