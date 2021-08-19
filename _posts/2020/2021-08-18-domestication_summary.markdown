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

### 6. Pairwise *ZF*<sub>ST</sub> contrasts

- contrast between Maine wild and selected line

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/MEW_MES_noinvers.sliding.zfst.hudson.jpg" alt="img" width="800"/>

- contrast between Long Island Sound wild and selected line

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/LIW_LIS_noinvers.sliding.zfst.hudson.jpg" alt="img" width="800"/>

- contrast between Delaware wild and selected line

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/DBW_DBS_noinvers.sliding.zfst.hudson.jpg" alt="img" width="800"/>

- contrast between North Carolina wild and selected line

<img src="https://github.com/hzz0024/hzz0024.github.io/blob/master/images/Dom/NCW_NCS_noinvers.sliding.zfst.hudson.jpg" alt="img" width="800"/>

Among them only Long Island Sound contrast and Delaware contrast have common shared ZFst window outliers, 

| Chromosome    | Start     | End       | No. genes   | Gene                                                                                              | Gene Id               |
|-------------  |---------- |---------- |-----------  |-------------------------------------------------------------------------------------------------- |---------------------  |
| NC_035784.1   | 90424000  | 90434000  | 2           | uncharacterized LOC111135324;uncharacterized LOC111135206                                         | 111135324;111135206   |
| NC_035784.1   | 91382000  | 91392000  | 1           | uncharacterized LOC111135518                                                                      | 111135518             |
| NC_035784.1   | 93210000  | 93226000  | 2           | WD repeat-containing protein 11-like;eukaryotic translation initiation   factor 3 subunit A-like  | 111134227;111134217   |
| NC_035784.1   | 93430000  | 93448000  | 1           | inositol 1,4,5-trisphosphate receptor type 3-like                                                 | 111138112             |
| NC_035784.1   | 94242000  | 94258000  | 1           | uncharacterized LOC111135873                                                                      | 111135873             |
| NC_035784.1   | 94682000  | 94692000  | 2           | F-box/WD repeat-containing protein 4-like; cilia- and flagella-associated   protein 43-like       | 111137357;111134972   |
| NC_035784.1   | 95324000  | 95340000  | 1           | uncharacterized LOC111134781                                                                      | 111134781             |
| NC_035784.1   | 95722000  | 95736000  | 2           | degenerin mec-10-like;uncharacterized LOC111132886                                                | 111132886;111132885   |
| NC_035784.1   | 95876000  | 95894000  | 2           | protein MAK16 homolog A-like;lymphocyte cytosolic protein 2-like                                  | 111136220;111099040   |
| NC_035784.1   | 98636000  | 98646000  | 1           | uncharacterized LOC111132793                                                                      | 111132793             |
| NC_035788.1   | 57378000  | 57386000  | 1           | multiple epidermal growth factor-like domains protein 10                                          | 111113010             |