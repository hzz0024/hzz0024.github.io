---
comments: true
title: DelBay1_angsd_SFS_4 
date: '2020-03-19 12:00'
tags:
  - angsd
  - SFS
  - fst
  - paralogs
  - WGS
categories:
  - WGS data analysis
---

This is an update about the DelBay1_angsd_SFS_FST post. One major goals of this post is to check whether the inversions exist in both chr 5 and 6, or just in chr 6. This is done by removing the potential inversion sites (generated from NYC data) and checking the PCA results.

⋅⋅* extract the loci with no inversion at chr 6 and chr 5_6

```shell
# tail -n +2 will only remove the first header line
zcat ch_no16inv_minI8D8maxD16.mafs.gz | tail -n +2 > FILE.tmp && mv FILE.tmp ch_no16inv_minI8D8maxD16_base_snplist
awk '{print $1,$2,$3,$4}' ch_no16inv_minI8D8maxD16_base_snplist > ch_no16inv_minI8D8maxD16_base_snplist_4col
awk '(NR == 1 ); !(($1 == "NC_035785.1") && ($2 > 29900000) && ($2 < 44500000))' ch_no16inv_minI8D8maxD16_base_snplist_4col > ch_base_4col_nochr6invers.snplist
awk '(NR == 1 ); !(($1 == "NC_035784.1") && ($2 > 60600000) && ($2 < 80200000))' ch_base_4col_nochr6invers.snplist > ch_base_4col_nochr56invers.snplist

zcat ref_no16inv_minI8D8maxD16.mafs.gz | tail -n +2 > FILE.tmp && mv FILE.tmp ref_no16inv_minI8D8maxD16_base_snplist
awk '{print $1,$2,$3,$4}' ref_no16inv_minI8D8maxD16_base_snplist > ref_no16inv_minI8D8maxD16_base_snplist_4col
awk '(NR == 1 ); !(($1 == "NC_035785.1") && ($2 > 29900000) && ($2 < 44500000))' ref_no16inv_minI8D8maxD16_base_snplist_4col > ref_base_4col_nochr6invers.snplist
awk '(NR == 1 ); !(($1 == "NC_035784.1") && ($2 > 60600000) && ($2 < 80200000))' ref_base_4col_nochr6invers.snplist > ref_base_4col_nochr56invers.snplist
```

⋅⋅* run angsd using the extracted snp list

```shell
angsd -b ch_16.bamlist -anc cv30.fa -out Fst_invers/ch_no16inv_minI8D8maxD16_MQ20_nochr6_invers -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -sites ch_base_4col_nochr6invers.snplist -rf ch5_6.list

angsd -b ch_16.bamlist -anc cv30.fa -out Fst_invers/ch_no16inv_minI8D8maxD16_MQ20_nochr56_invers -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -sites ch_base_4col_nochr56invers.snplist -rf ch5_6.list

angsd -b ref_16.bamlist -anc cv30.fa -out Fst_invers/ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -sites ref_base_4col_nochr6invers.snplist -rf ch5_6.list

angsd -b ref_16.bamlist -anc cv30.fa -out Fst_invers/ref_no16inv_minI8D8maxD16_MQ20_nochr56_invers -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -sites ref_base_4col_nochr56invers.snplist -rf ch5_6.list
```

note: -rf need to specifed during angsd running, otherwise the output saf files will cause errors during the realSFS step

error message: Read block operation failed with error -1 after 5352952 of 18480688 bytes & Problem reading chunk in bgzf_read

---
### Results

⋅⋅* ch_no16inv_minI8D8maxD16_MQ20_nochr6_invers

Total number of sites analyzed: 109553033\
Number of sites retained after filtering: 6718359\
Time used =  1587.00 sec

⋅⋅* ch_no16inv_minI8D8maxD16_MQ20_nochr56_invers

Total number of sites analyzed: 109553033\
Number of sites retained after filtering: 5411893\
Time used =  1383.00 sec

⋅⋅* ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers

Total number of sites analyzed: 111768475\
Number of sites retained after filtering: 4728400\
Time used =  1503.00 sec

⋅⋅* ref_no16inv_minI8D8maxD16_MQ20_nochr56_invers

Total number of sites analyzed: 111768475\
Number of sites retained after filtering: 3876828\
Time used =  1353.00 sec

---
#### Fst plots used to check if the potential inversion sites are removed or not

#### Control 
<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_single_SNP_base.jpg" alt="img" width="800"/>

#### No chr 6 inversion
<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers.jpg" alt="img" width="800"/>

#### not chr 5-6 inversion 
<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_no16inv_minI8D8maxD16_MQ20_nochr56_invers.jpg" alt="img" width="800"/>

#### PCA plots for each angsd run (ch and ref seperate datasets)

#### ch_no16inv_minI8D8maxD16_MQ20_nochr6_invers.pc1-2, 2-3, and 3-4
<img src="https://hzz0024.github.io/images/ch_no16inv_minI8D8maxD16_MQ20_nochr6_invers.pc1-2.jpg" alt="img" width="800"/>
<img src="https://hzz0024.github.io/images/ch_no16inv_minI8D8maxD16_MQ20_nochr6_invers.pc2-3.jpg" alt="img" width="800"/>
<img src="https://hzz0024.github.io/images/ch_no16inv_minI8D8maxD16_MQ20_nochr6_invers.pc3-4.jpg" alt="img" width="800"/>

#### ch_no16inv_minI8D8maxD16_MQ20_nochr56_invers.pc1-2, 2-3, and 3-4
<img src="https://hzz0024.github.io/images/ch_no16inv_minI8D8maxD16_MQ20_nochr56_invers.pc1-2.jpg" alt="img" width="800"/>
<img src="https://hzz0024.github.io/images/ch_no16inv_minI8D8maxD16_MQ20_nochr56_invers.pc2-3.jpg" alt="img" width="800"/>
<img src="https://hzz0024.github.io/images/ch_no16inv_minI8D8maxD16_MQ20_nochr56_invers.pc3-4.jpg" alt="img" width="800"/>

#### ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers.pc1-2, 2-3, and 3-4
<img src="https://hzz0024.github.io/images/ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers.pc1-2.jpg" alt="img" width="800"/>
<img src="https://hzz0024.github.io/images/ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers.pc2-3.jpg" alt="img" width="800"/>
<img src="https://hzz0024.github.io/images/ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers.pc3-4.jpg" alt="img" width="800"/>

#### ref_no16inv_minI8D8maxD16_MQ20_nochr56_invers.pc1-2, 2-3, and 3-4
<img src="https://hzz0024.github.io/images/ref_no16inv_minI8D8maxD16_MQ20_nochr56_invers.pc1-2.jpg" alt="img" width="800"/>
<img src="https://hzz0024.github.io/images/ref_no16inv_minI8D8maxD16_MQ20_nochr56_invers.pc2-3.jpg" alt="img" width="800"/>
<img src="https://hzz0024.github.io/images/ref_no16inv_minI8D8maxD16_MQ20_nochr56_invers.pc3-4.jpg" alt="img" width="800"/>

---
Conclusion



