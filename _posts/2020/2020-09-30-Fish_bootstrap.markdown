---
comments: true
title: Bootstrap analyses around outliers
date: '2020-09-30 12:00'
tags:
  - DelBay19
  - ouliter
  - bootstrap
  - WGS
  - Fisher's exact
  - 
categories:
  - WGS data analysis
---

After identifing some potential outliers using the Fisher's approach, I performed the outlier analyses at the window scale. 

The detailed method is described in [Dixon et al. (2015) Genomic determinants of coral heat tolerance across latitudes](https://science.sciencemag.org/content/348/6242/1460). 

Below are the initial results

Manhattan plot for CH vs. REF + HC vs. SR. Note that the Red dots indicates outliers < FDR 5% threshold (No. of outliers = 10). Light green bars identify
regions with significant clustering of peaks (FDR < 0.1), dark green bars identify regions with significant clustering of peaks (FDR < 0.05, according to 100,000 bootstrapped replicates).

<img src="https://hzz0024.github.io/images/Fish_boot/hard_REF-CH-SR-HC_bootstrap.jpg" alt="img" width="800"/>

Manhattan plot for SR vs. REF + COH vs. ARN. No outliers < FDR 5% in this case.  

<img src="https://hzz0024.github.io/images/Fish_boot/hard_SR-REF-COH-ARN_bootstrap.jpg" alt="img" width="800"/>

Looks like the bootstrap analysis need further optimization. At least those green bars should not show in the control group (SR vs. REF + COH vs. ARN).

However, I noticed that in Dixon et al, the combined p values were used for bootstrap analyses but not the adjust ones, while their Fisher's outliers were uncovered after fdr correction. This sounds an unfair comparison between two approaches.

Instead, I tried to use the p-values AFTER fdr adjustment as the target values for bootstrap analysis, which produce more meaningful results shown below,

- REF-CH-SR-HC bootstrap results with 15 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/REF-CH-SR-HC_15_bootstrap.jpg" alt="img" width="800"/>

- REF-CH-SR-HC bootstrap results with 25 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/REF-CH-SR-HC_25_bootstrap.jpg" alt="img" width="800"/>

- REF-CH-SR-HC bootstrap results with 50 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/REF-CH-SR-HC_50_bootstrap.jpg" alt="img" width="800"/>

- REF-CH-SR-HC bootstrap results with 150 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/REF-CH-SR-HC_150_bootstrap.jpg" alt="img" width="800"/>

- REF-CH-SR-HC bootstrap results with 1500 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/REF-CH-SR-HC_1500_bootstrap.jpg" alt="img" width="800"/>

- REF-CH-SR-HC bootstrap results with 10000 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/REF-CH-SR-HC_10000_bootstrap.jpg" alt="img" width="800"/>

- REF-CH-NB-HC bootstrap results with 15 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/REF-CH-NB-HC_15_bootstrap.jpg" alt="img" width="800"/>

- REF-CH-NB-HC bootstrap results with 25 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/REF-CH-NB-HC_25_bootstrap.jpg" alt="img" width="800"/>

- REF-CH-NB-HC bootstrap results with 50 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/REF-CH-NB-HC_50_bootstrap.jpg" alt="img" width="800"/>

- REF-CH-NB-HC bootstrap results with 150 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/REF-CH-NB-HC_150_bootstrap.jpg" alt="img" width="800"/>

- REF-CH-NB-HC bootstrap results with 1500 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/REF-CH-NB-HC_1500_bootstrap.jpg" alt="img" width="800"/>

- REF-CH-NB-HC bootstrap results with 10000 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/REF-CH-NB-HC_10000_bootstrap.jpg" alt="img" width="800"/>

- control group (SR-REF-COH-ARN) bootstrap results with 15 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/SR-REF-COH-ARN_15_bootstrap.jpg" alt="img" width="800"/>

- control group (SR-REF-COH-ARN) bootstrap results with 25 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/SR-REF-COH-ARN_25_bootstrap.jpg" alt="img" width="800"/>

- control group (SR-REF-COH-ARN) bootstrap results with 50 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/SR-REF-COH-ARN_50_bootstrap.jpg" alt="img" width="800"/>

- control group (SR-REF-COH-ARN) bootstrap results with 150 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/SR-REF-COH-ARN_150_bootstrap.jpg" alt="img" width="800"/>

- control group (SR-REF-COH-ARN) bootstrap results with 1500 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/SR-REF-COH-ARN_1500_bootstrap.jpg" alt="img" width="800"/>

- control group (SR-REF-COH-ARN) bootstrap results with 10000 SNP/window

<img src="https://hzz0024.github.io/images/Fish_boot/SR-REF-COH-ARN_10000_bootstrap.jpg" alt="img" width="800"/>

---

In order to pick up the best window size, I calculate the average length for multiple window size (i.e. how many SNP per window but not the genomic span length),

```sh
for i in {0..9}; do
echo "chr $i"
cat ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 |grep "NC_03578$i.1" | awk -F ' ' '{print $2}' | head -n 1
cat ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 |grep "NC_03578$i.1" | awk -F ' ' '{print $2}' | tail -n 1
done
```

| Chromosome   |  start      |    end    | total bp  |    
|--------------|-------------|-----------|-----------|
|1             | 1466        | 65640332  | 65638866  |
|2             | 23617       | 61743264  | 61719647  |
|3             | 44088       | 76219443  | 76175355  |
|4             | 146974      | 59618701  | 59471727  |
|5             | 212281      | 98684645  | 78872364  |
|6             | 67761       | 51092571  | 31424810  |
|7             | 10954       | 57827871  | 57816917  |
|8             | 50363       | 75634506  | 75584143  |
|9             | 230804      | 104102130 | 103871326 |
|10            | 196012      | 32597966  | 32401954  |


The total SNPs spans 65638866+61719647+76175355+59471727+(98472364-19600000)+(51024810-14600000)+57816917+75584143+103871326+32401954=647977109 bp = 647.98 Mb

Average length per SNP is 647977109/1934038 = 335.038 bp

According to this,

15000 SNP ~ 5M bp
10000 SNP ~ 3M bp
1500 SNP ~ 500k bp
50 SNP ~ 17k bp
25 SNP ~ 8k bp
15 SNP ~ 5k bp

The coral paper used 15 SNP/window, corresponding (475 Mb/1448) * 15 ~ 4.92M bp/window, which equals ~ 15000 SNPs/window in my dataset (note genome size of coral *Acropora millepora* is recently published by [Fuller et al. 2020. Population genetics of the coral Acropora millepora: Toward genomic prediction of bleaching](https://science.sciencemag.org/content/369/6501/eaba4674?rss=1))

However, this is just rough estimate the distance among SNPs. I wrote an R script to calculate the average distance (bp) for each SNP window (e.g. 15 SNPs/window). The script can be found in the *DelBay_project/R_script/Calculate_window_size/cal_window_size.R*

- Average distance (bp) for each SNP window

| Chromosome   | No. SNP     |15 SNP/window|25 SNP/window |50 SNP/window|150 SNP/window|1500 SNP/window|10000 SNP/window |
|--------------|-------------|-------------|--------------|-------------|--------------|---------------|-----------------|
|1             | 273241      | 3603        | 6013         | 12041       |  36442       |  361938       |  2370021        |
|2             | 294774      | 3155        | 5234         | 10487       |  31483       |  313958       |  2159153        |
|3             | 315491      | 3621        | 6057         | 12072       |  36247       |  361699       |  2423132        |
|4             | 283496      | 3147        | 5246         | 10488       |  31466       |  314667       |  2076772        |
|5             | 408592      | 3625        | 6024         | 12067       |  36150       |  362800       |  2432094        |
|6             | 14784       | 51783       | 86249        | 172387      | 521769       | 5153410       | 30395754        |
|7             | 75924       | 11421       | 19037        |  38134      | 115216       | 1134255       | 7558580         |
|8             | 107832      |10515        |  17542       | 35121       | 105259       | 1050089       |  6964755        |
|9             | 132474      | 11760       | 19602        | 39351       | 119879       | 1172390       | 8060879         |
|10            | 27430       | 17715       | 29511        | 59023       | 177082       | 1715049       | 11492948        |    

From this table we can see that the first 5 chromosomes have relatively higher SNP density and shorter distance per window. 

So far, I'm not sure what might the best window size. Perhaps LD decay pattern is a good indicator for the window size choice. Given that this is a new dataset from the reproduced Angsd global SNP list. I conducted a LD decay analysis for chromosome 5, with a window size of 25 and 50 SNPs.

```sh
Rscript --vanilla --slave fit_LDdecay.R --ld_files LD.list --out ALL_chr5_k25_max10k.jpg --fit_level 0 --max_kb_dist 50 --fit_level 50 --plot_size 1.5,2
Rscript --vanilla --slave fit_LDdecay.R --ld_files LD.list --out ALL_chr5_k50_max10k.jpg --fit_level 0 --max_kb_dist 50 --fit_level 50 --plot_size 1.5,2
```

<img src="https://hzz0024.github.io/images/ngsLD/ALL_chr5_k25_max10k.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/ngsLD/ALL_chr5_k50_max10k.jpg" alt="img" width="800"/>

Another way to determine the window size may be the average protein coding "gene" size. As mentioned in the recent paper by Stern and Lee [Evolutionary origins of genomic adaptations in an invasive copepod](https://static-content.springer.com/esm/art%3A10.1038%2Fs41559-020-1201-y/MediaObjects/41559_2020_1201_MOESM1_ESM.pdf), "The window size of 10 kb was chosen capture selection targets on the same ‘gene’ and surrounding genomic region, given an average protein coding ‘gene’ size of ~8.5 kb in the *E. affinis* genome."

For eastern oyster, the average "gene" size is 10,828 bp [https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Crassostrea_virginica/100/](https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Crassostrea_virginica/100/), which roughly equals 25 SNP/window.




