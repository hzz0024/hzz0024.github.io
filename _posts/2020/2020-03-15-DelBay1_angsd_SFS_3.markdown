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

THis is an update about the DelBay1_angsd_SFS_FST post. I conducted a small test to explore the effects of different angsd options on SNP counts and the Fst.

Three parameters options are:

-setMaxDepth

-minMapQ

-ref

The test is only focused on the chromosome 5-6, the regions that may have genome inversion in our target species.

code for testing:

```shell
angsd -b ref_16.bamlist -anc cv30.fa -out Fst_test/ref_no16inv_minI8D8maxD16_MQ20 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -rf ch5_6.list

angsd -b ref_16.bamlist -anc cv30.fa -out Fst_test/ref_no16inv_minI8D8maxD32_MQ20 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 32 -minInd 8 -rf ch5_6.list

angsd -b ref_16.bamlist -anc cv30.fa -out Fst_test/ref_no16inv_minI8D8maxD16 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -rf ch5_6.list

angsd -b ref_16.bamlist -anc cv30.fa -ref cv30.fa -out Fst_test/ref_no16inv_minI8D8maxD16_MQ20_ref -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -rf ch5_6.list
```

---
Results

#### SNP counts

1. test the results between -setMaxDepth 16 and -setMaxDepth 32 (which reflects 1X and 2X maximum coverage)

| -setMaxDepth |   SNP No.  | 
|--------------|------------|
|      16      |   5225855  | 
|      32      |  25978065  | 

cpu-time used =  1574.93 sec vs. 3290.04 sec

2. test what if remove the mapping quality control

|  -minMapQ  |   SNP No.  | 
|------------|------------|
|      20    |   5225855  | 
|      not specified     |   5225855  | 

cpu-time used =  1574.93 sec vs. 1570.83 sec

3. test what if both ancestral and reference genomes are both specified 

|      -ref         |   SNP No.  | 
|-------------------|------------|
|   not specified   |   5225855  | 
|   specified       |   5225855  | 

cpu-time used =  1574.93 sec vs. 1597.00 sec


---
Conclusion


