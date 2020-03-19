---
comments: true
title: DelBay1_angsd_SFS_3 
date: '2020-03-15 12:00'
tags:
  - angsd
  - SFS
  - fst
  - paralogs
  - WGS
categories:
  - WGS data analysis
---

This is an update about the DelBay1_angsd_SFS_FST post. I conducted a small test to explore the effects of different angsd options on SNP counts and the Fst.

Three parameters options are:

-setMaxDepth

-minMapQ

-ref

The test is only tested on the chromosome 5-6, the regions that may have genome inversion in our target species. I include 32 samples in this test, 16 from challenge and 16 from reference populations (run independently for beagle file and joined later for SFS estimates). 

code for testing:

```shell
# Base 
angsd -b ref_16.bamlist -anc cv30.fa -out Fst_test/ref_no16inv_minI8D8maxD16_MQ20 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -rf ch5_6.list

# double the maximum depth by -setMaxDepth 32
angsd -b ref_16.bamlist -anc cv30.fa -out Fst_test/ref_no16inv_minI8D8maxD32_MQ20 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 32 -minInd 8 -rf ch5_6.list

# not specifying the -minMapQ 20 
angsd -b ref_16.bamlist -anc cv30.fa -out Fst_test/ref_no16inv_minI8D8maxD16 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -rf ch5_6.list

# giving both reference and ancestral genome
angsd -b ref_16.bamlist -anc cv30.fa -ref cv30.fa -out Fst_test/ref_no16inv_minI8D8maxD16_MQ20_ref -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -rf ch5_6.list
```

---
Results

#### test the results between -setMaxDepth 16 and -setMaxDepth 32 (which reflects 1X and 2X maximum coverage)

| -setMaxDepth |   SNP No.  |Fst(unwighted/wighted)| SNPs shared by ch and ref in ch5| SNPs shared by ch and ref in ch6|
|--------------|------------|-----------------|---------------------------------|---------------------------------|
|    16 (base) |   5225855  |0.022380/0.014844| 			964146				  | 			245704				| 
|      32      |  25978065  |0.008321/0.009190| 	   13869518				  |			 3048976				|  

Time used for ref run =  1263.00 sec vs. 2888.00 sec
Time used for ch run = 1356.00 sec vs. 3467.00 sec

#### test what if remove the mapping quality control

|  -minMapQ    |   SNP No.  |Fst(unwighted/wighted)|SNPs shared by ch and ref in ch5| SNPs shared by ch and ref in ch6|
|--------------|------------|-----------------|--------------------------------|---------------------------------|
|   20 (base)  |   5225855  |0.022380/0.014844|			964146				  | 			245704				| 
|not specified |   5225855  |0.024870/0.014737|	        964146			|	   245704				| 
 
Time used for ref run =  1263.00 sec vs. 1180.00 sec
Time used for ch run = 1356.00 sec vs. 1311.00 sec

#### test what if both ancestral and reference genomes are both specified 

|  	-ref       |   SNP No.  |Fst(unwighted/wighted)|SNPs shared by ch and ref in ch5| SNPs shared by ch and ref in ch6|
|--------------|------------|-----------------|--------------------------------|---------------------------------|
|not specified (base)|   5225855  |0.022380/0.014844|			964146				  | 			245704				| 
|specified           |   5225855  |0.024743/0.014759|	        964146				  |			    245704				| 

Time used for ref run =  1263.00 sec vs. 1223.00 sec
Time used for ch run = 1356.00 sec vs. 1281.00 sec

#### Fst plots along chr 5-6

#### Base 
<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_single_SNP_base.jpg" alt="img" width="800"/>

#### double the maximum depth by -setMaxDepth 32
<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_single_SNP_depth.jpg" alt="img" width="800"/>

#### not specifying the -minMapQ 20 
<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_single_SNP_mq20.jpg" alt="img" width="800"/>

#### giving both reference and ancestral genome
<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_single_SNP_ref.jpg" alt="img" width="800"/>

---
Conclusion

1. parameter -setMaxDepth impacts a lot on the number of SNPs, and the average Fst values. By changing the -setMaxDepth from 16 (1x) to 32 (2x), we can see more high-Fst SNPs and clearer inversion patterns in the chr 6

2. -minMapQ has limited impacts on the Fst (differ by 0.0025 from Base result), but keep this option for data generation may increase the SNP quality

3. not sure how -ref will impact the results

