---
layout: post
title: DelBay19 ANGSD trial
date: '2020-02-16 12:00'
tags:
  - SNP
  - test
  - angsd
  - filter
  - WGS
categories:
  - WGS data analysis
---

In the last few days I was looking for the answer to "what should be our minimum read count needed from each sample?", this could be done by

1) estimating the sample coverage using filtered bam files (MQ20), and determining the relative proportions of reads needed for rerun

2) practicing with ANGSD to give the insight about the relative proportions of reads needed for rerun.  

The second goal can be subdivided into multiple directions,

1) key questions is how many SNPs I get for a given filter set (like Matt used with NYC samples). I will bin the samples by read count and then run ANGSD on those bins of samples with a uniform filter set, to see how read count translates into relative post-filter SNP numbers. 

2) performaing initial test for PCA population structure and looking for inversions.

3) checking the SNP density, the ideal value should be >= 2% (2 SNPs every 100 bp). The density should be calculated from the mappable genome length, not the total genome length. 

below is the script for ANGSD running of sample bins,

```shell

# angsd run for 7 samples with reads counts at 5000k-1M, 1M-2M, 2M-3M, 3M-4M .....17-18M categories

/programs/angsd20191002/angsd/angsd -b 500K_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/500K_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_500K_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 1M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/1M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_1M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 2M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/2M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_2M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 3M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/3M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_3M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 4M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/4M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_4M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 5M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/5M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_5M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 6M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/6M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_6M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 7M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/7M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_7M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 8M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/8M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_8M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 9M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/9M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_9M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 10M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/10M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_10M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 11M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/11M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_11M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 12M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/12K_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_12K_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 13M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/13M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_13M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 14M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/14M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_14M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 15M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/15M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_15M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 16M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/16M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_16M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 17M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/17M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_17M_7_minMAF05.log &
/programs/angsd20191002/angsd/angsd -b 18M_7.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_str/18M_7_D4maxD60_minQ20_minMAF05_SNPe6_no2inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 4 -setMaxDepth 60 -minInd 2 -doMajorMinor 1 >& test_18M_7_minMAF05.log &

# running angsd for the wild and challenge vs. reference 

# angsd run for 227 wild samples
# the parameters allow each SNP has a 75% individual coverage, with some stringent settings (e.g. -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6). A maximum total depth threshold of 450 to remove loci from repetitive regions.

angsd -b wild_227.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_wild/wild_227_D170maxD450_minQ20_minMAF05_SNPe6_no227inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 16 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 170 -setMaxDepth 450 -minInd 171 -doMajorMinor 1 > minMaf05_wild_227_minMAF05.log

# angsd run for 97 challenge vs. reference samples

angsd -b challenge_97.bamlist -anc cv30.fa -ref cv30.fa -out minMaf05_challenge/challenge_97_D73maxD200_minQ20_minMAF05_SNPe6_no97inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 16 -minQ 20 -minMapQ 20 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 73 -setMaxDepth 200 -minInd 73 -doMajorMinor 1 > minMaf05_challenge_97_minMAF05.log


```
- useful links

[Best ANGSD tutorial](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md)

---

#### RESULTS

Goal 1: 

[sample coverage](https://docs.google.com/spreadsheets/d/10V7vTdNp7oagq4SlPPfOGA-kgmrmh4x6m4olKCdzB6E/edit#gid=1728449447)

Goal 2:
I generated a plot to reflect the SNP number in each sample bin (based on read counts). From the plot we can see that as reads increases, the SNP number steadily increases (in a linear fashing). Therefore, the SNP number data may not be useful to estimate how many reads we need for the rerun.  

<img src="https://hzz0024.github.io/images/2020-02-17.jpeg" alt="img" width="800"/>
