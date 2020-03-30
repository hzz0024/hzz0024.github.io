---
comments: true
title: DelBay19 check inversions by chromosomes 
date: '2020-03-29 12:00'
tags:
  - angsd
  - SFS
  - PCA
  - WGS
categories:
  - WGS data analysis
---

The major goal of this post is to check whether the inversions exist in ohter chromosomes except the chr 5 and 6. This is done by checking the PCA results and the Fst plots. The challen datasets, which includes 16 challenge (ch) and 16 reference (ref) is used for the test. 

### PCA plots for each chr

- beagle file generation

```sh
module load angsd/0.931

for i in {1..10}; do
angsd -b ch_ref_32.bamlist \
-anc cv30.fa \
-out 'invers_by_chr/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6'_$i \
-dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -doIBS 1 -makematrix 1 -doCov 1 \
-P 12 \
-minQ 20 -minMapQ 20 \
-setMinDepth 16 \
-setMaxDepth 64 \
-minInd 16 \
-remove_bads 1 \
-uniqueOnly 1 \
-only_proper_pairs 1 \
-minMaf 0.05 \
-SNP_pval 1e-6 \
-rf 'chr'$i.list
done
```
- .cov.npy file generation
```sh
module load pcangsd/0.98

for i in {1..10}; do
python3 /tools/pcangsd-0.98/pcangsd.py \
-beagle 'ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_'$i.beagle.gz \
-minMaf 0.05 \
-threads 16 \
-o 'ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_'$i \
> ch_ref_check_inversions_by_chr.log
done
```
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

- chr 1
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_1.pc1-2-1.jpg" alt="img" width="800"/>

- chr 2
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_1.pc1-2-1.jpg" alt="img" width="800"/>

- chr 3
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_1.pc1-2-1.jpg" alt="img" width="800"/>

- chr 4
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_1.pc1-2-1.jpg" alt="img" width="800"/>

- chr 5
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_1.pc1-2-1.jpg" alt="img" width="800"/>

- chr 6
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_1.pc1-2-1.jpg" alt="img" width="800"/>

- chr 7
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_1.pc1-2-1.jpg" alt="img" width="800"/>

- chr 8
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_1.pc1-2-1.jpg" alt="img" width="800"/>

- chr 9
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_1.pc1-2-1.jpg" alt="img" width="800"/>

- chr 10
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_1.pc1-2-1.jpg" alt="img" width="800"/>

---

Conclusion

1. Based on the PCA plots alone, it seems the inversions only exist in the chromosome 5 and 6

2. 
 





