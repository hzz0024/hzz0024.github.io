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

The major goal of this post is to check whether the inversions exist in other chromosomes except the chromosome 5 and 6. This was done by checking the PCA results and the Fst plots. The challenge dataset with 16 challenge (ch) and 16 reference (ref) was used for this test. 

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

- chr 1

Total number of sites analyzed: 52552080  
Number of sites retained after filtering: 444589  
Time used =  1512.00 sec

- chr 2

Total number of sites analyzed: 49791991  
Number of sites retained after filtering: 401766  
Time used =  1314.00 sec

- chr 3

Total number of sites analyzed: 59904312  
Number of sites retained after filtering: 423224  
Time used =  1708.00 sec

- chr 4

Total number of sites analyzed: 47970660  
Number of sites retained after filtering: 399808  
Time used =  1512.00 sec

- chr 5

Total number of sites analyzed: 83171346  
Number of sites retained after filtering: 783469  
Time used =  1958.00 sec

- chr 6

Total number of sites analyzed: 36791742  
Number of sites retained after filtering: 169199  
Time used =  747.00 sec

- chr 7

Total number of sites analyzed: 42934065  
Number of sites retained after filtering: 198263  
Time used =  843.00 sec

- chr 8

Total number of sites analyzed: 56545474  
Number of sites retained after filtering: 315429  
Time used =  1200.00 sec

- chr 9

Total number of sites analyzed: 76515152  
Number of sites retained after filtering: 386866  
Time used =  1589.00 sec

- chr 10

Total number of sites analyzed: 23311866  
Number of sites retained after filtering: 105533  
Time used =  431.00 sec

---
#### PCA plots

- chr 1
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_1.pc1-2-1.jpg" alt="img" width="800"/>

- chr 2
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_2.pc1-2-1.jpg" alt="img" width="800"/>

- chr 3
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_3.pc1-2-1.jpg" alt="img" width="800"/>

- chr 4
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_4.pc1-2-1.jpg" alt="img" width="800"/>

- chr 5
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_5.pc1-2-1.jpg" alt="img" width="800"/>

- chr 6
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_6.pc1-2-1.jpg" alt="img" width="800"/>

- chr 7
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_7.pc1-2-1.jpg" alt="img" width="800"/>

- chr 8
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_8.pc1-2-1.jpg" alt="img" width="800"/>

- chr 9
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_9.pc1-2-1.jpg" alt="img" width="800"/>

- chr 10
<img src="https://hzz0024.github.io/images/ch_ref_no32inv_minI16D16maxD64_MQ20_minMAF05_SNPe6_10.pc1-2-1.jpg" alt="img" width="800"/>

---

Conclusion

1. Based on the PCA plots alone, it seems the inversions only exist in the chromosome 5 and 6

2. 
 





