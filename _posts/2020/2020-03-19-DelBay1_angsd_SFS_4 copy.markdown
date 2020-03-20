---
comments: true
title: DelBay1_inversion 
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

The major goal of this post is to check whether the inversions exist in both chr 5 and 6, or just in chr 6. This is done by removing the potential inversion sites (generated from NYC data) and checking the PCA results. I used two seperate datasets, challenge (ch, 16 samples) and reference (ref, 16 samples), for SFS data creation and Fst plotting, and ran PCA on these two seperate dataset. After that, I ran angsd on the combined dataset (32 samples) and check the impact of inversion on PCA patterns.

- extract the loci with no inversion at chr 6 and chr 5_6

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

- run angsd using the extracted snp list

```shell
angsd -b ch_16.bamlist -anc cv30.fa -out Fst_invers/ch_no16inv_minI8D8maxD16_MQ20_nochr6_invers -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -sites ch_base_4col_nochr6invers.snplist -rf ch5_6.list

angsd -b ch_16.bamlist -anc cv30.fa -out Fst_invers/ch_no16inv_minI8D8maxD16_MQ20_nochr56_invers -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -sites ch_base_4col_nochr56invers.snplist -rf ch5_6.list

angsd -b ref_16.bamlist -anc cv30.fa -out Fst_invers/ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -sites ref_base_4col_nochr6invers.snplist -rf ch5_6.list

angsd -b ref_16.bamlist -anc cv30.fa -out Fst_invers/ref_no16inv_minI8D8maxD16_MQ20_nochr56_invers -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 8 -setMaxDepth 16 -minInd 8 -sites ref_base_4col_nochr56invers.snplist -rf ch5_6.list
```

note: -rf need to specifed during angsd running, otherwise the output saf files will cause errors during the realSFS step

error message: Read block operation failed with error -1 after 5352952 of 18480688 bytes & Problem reading chunk in bgzf_read

---
#### Results

- ch_no16inv_minI8D8maxD16_MQ20_nochr6_invers

Total number of sites analyzed: 109553033  
Number of sites retained after filtering: 6718359  
Time used =  1587.00 sec

- ch_no16inv_minI8D8maxD16_MQ20_nochr56_invers

Total number of sites analyzed: 109553033  
Number of sites retained after filtering: 5411893  
Time used =  1383.00 sec

- ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers

Total number of sites analyzed: 111768475  
Number of sites retained after filtering: 4728400  
Time used =  1503.00 sec

- ref_no16inv_minI8D8maxD16_MQ20_nochr56_invers

Total number of sites analyzed: 111768475  
Number of sites retained after filtering: 3876828  
Time used =  1353.00 sec

---
#### Fst plots used to check if the potential inversion sites are removed or not

- Control 
<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_single_SNP_base.jpg" alt="img" width="800"/>

- Without chr 6 inversion
<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers.jpg" alt="img" width="800"/>

- Withou chr 5-6 inversion 
<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_no16inv_minI8D8maxD16_MQ20_nochr56_invers.jpg" alt="img" width="800"/>

#### PCA plots for each angsd run (ch and ref seperate datasets)

- ch_no16inv_minI8D8maxD16_MQ20_nochr6_invers.pc1-2
<img src="https://hzz0024.github.io/images/ch_no16inv_minI8D8maxD16_MQ20_nochr6_invers.pc1-2.jpg" alt="img" width="800"/>

- ch_no16inv_minI8D8maxD16_MQ20_nochr56_invers.pc1-2
<img src="https://hzz0024.github.io/images/ch_no16inv_minI8D8maxD16_MQ20_nochr56_invers.pc1-2.jpg" alt="img" width="800"/>

- ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers.pc1-2
<img src="https://hzz0024.github.io/images/ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers.pc1-2.jpg" alt="img" width="800"/>

- ref_no16inv_minI8D8maxD16_MQ20_nochr56_invers.pc1-2
<img src="https://hzz0024.github.io/images/ref_no16inv_minI8D8maxD16_MQ20_nochr56_invers.pc1-2.jpg" alt="img" width="800"/>

---
Then I ran the angsd using 32 ch_ref individuals to check the pca patterns before/after inversion exclusion. This part I added some stringent parameters to filter out non-polymorphic and low-quality loci, these include -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -minMaf 0.05 -SNP_pval 1e-6.

- run angsd  
```shell
angsd -b ch_ref_32.bamlist -anc cv30.fa -out pca_invers/ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 16 -minQ 20 -minMapQ 20 -setMinDepth 16 -setMaxDepth 32 -minInd 16 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -minMaf 0.05 -SNP_pval 1e-6 -rf ch5_6.list
```

- remove potential inversions 

```shell

zcat ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6.mafs.gz | tail -n +2 > FILE.tmp && mv FILE.tmp ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_snplist

awk '{print $1,$2,$3,$4}' ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_snplist > ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_snplist_4col

awk '(NR == 1 ); !(($1 == "NC_035785.1") && ($2 > 29900000) && ($2 < 44500000))' ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_snplist_4col > ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_snplist_4col_nochr6invers.snplist

awk '(NR == 1 ); !(($1 == "NC_035784.1") && ($2 > 60600000) && ($2 < 80200000))' ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_snplist_4col_nochr6invers.snplist > ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_snplist_4col_nochr56invers.snplist

# index the snp sites

angsd sites index ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_snplist_4col_nochr6invers.snplist  
angsd sites index ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_snplist_4col_nochr56invers.snplist

cat ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_snplist_4col_nochr6invers.snplist | wc -l

# 117471 (117470 sites) 

cat ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_snplist_4col_nochr56invers.snplist | wc -l

# 92543 (92542 sites)

```

- run angsd again using extracted snp lists

```shell

# exclude the chr 6 inversion (ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_nochr6invers)

angsd -b ch_ref_32.bamlist -anc cv30.fa -out pca_invers/ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_nochr6invers -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 16 -minQ 20 -minMapQ 20 -setMinDepth 16 -setMaxDepth 32 -minInd 16 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -minMaf 0.05 -SNP_pval 1e-6  -sites ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_snplist_4col_nochr6invers.snplist -rf ch5_6.list

# exclude both chr 5 and 6 inversion (ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_nochr56invers)

angsd -b ch_ref_32.bamlist -anc cv30.fa -out pca_invers/ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_nochr56invers -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 -P 16 -minQ 20 -minMapQ 20 -setMinDepth 16 -setMaxDepth 32 -minInd 16 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -minMaf 0.05 -SNP_pval 1e-6  -sites ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_snplist_4col_nochr56invers.snplist -rf ch5_6.list

```


#### Results

- global run (ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6)

Total number of sites analyzed: 119963088  
Number of sites retained after filtering: 131325  
Time used =  1494.00 sec

- exclude the chr 6 inversion (ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_nochr6invers)

Total number of sites analyzed: 119963088  
Number of sites retained after filtering: 117470  
Time used =  1458.00 sec

- exclude both chr 5 and 6 inversion (ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_nochr56invers)

Total number of sites analyzed: 119963088  
Number of sites retained after filtering: 92541  
Time used =  1457.00 sec

#### PCA plots

- global run result PC 1-2 (ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6)
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6.pc1-2.jpg" alt="img" width="800"/>

- exclude the chr 6 inversion PC 1-2 (ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_nochr6invers)
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_nochr6invers.pc1-2.jpg" alt="img" width="800"/>

- exclude both chr 5 and 6 inversion PC 1-2 (ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_nochr56invers)
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD32_MQ20_minMAF05_SNPe6_nochr56invers.pc1-2.jpg" alt="img" width="800"/>

---

Conclusion

1. My previous PCA analyses using all challenge and reference samples have shown that the alternative homokaryote (BB) is rare (2 in 97 samples), most of the samples are either AA homokaryote or AB heterokaryote. This might explain the two major cluster I observed in the baseline PCA results.  

2. PCA plots without chr 5-6 inversions still show patterns of two genetic clusters, which may suggest the existence of additional inversion on chr 5 and 6. My next step is to blown up the chr 5 and identify the extra inversion regions.




