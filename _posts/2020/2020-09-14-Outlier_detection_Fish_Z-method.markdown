---
comments: true
title: DelBay19 outlier Fisher’s approach vs weighted Z-method
date: '2020-09-14 12:00'
tags:
  - DelBay19
  - ouliter
  - WGS
  - Fisher's exact
  - 
categories:
  - WGS data analysis
---

In this post I made some modifications in Fisher's exact tests. I also compared the results from Fisher's approach vs weighted Z-method.  

1) Modification on Fisher's approach. 

Previously I performed Fisher's approach on DelBay datasets. The replicates are set as challenge vs reference plus wild contrast: either Hope Creek (HC) vs Shell Rock (SR) or Hope Creek (HC) vs New Bed (NB). I used REF vs SR + COH vs ARN as the control group.  

However, I found an important error when I double check the codes in [Dixon et al. (2015) Genomic determinants of coral heat tolerance across latitudes](https://science.sciencemag.org/content/348/6242/1460). This happens when they tried to combine the p-values from replicate tests:

```sh
#FUNCTION fisher.method():
#performs fisher's method to combine a set of p values
#ARGUMENTS: ps = a vector of p values you want combined
fisher.method=function(ps) {
  ftest=-2*sum(log(ps))
  df=length(2*ps)
  pv=1-pchisq(q=ftest,df=df)
} 
```

Here the *df* sould be *df=2 x length(ps)* but not *df=length(2 x ps)*. The detailed explanation can be found in this post: [https://brainder.org/2012/05/11/the-logic-of-the-fisher-method-to-combine-p-values/]. I am not sure how this might impact the results in coral heat selection paper, but this does impact the number of outliers in my datasets.

2) According to [Whitlock (2005)](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1420-9101.2005.00917.x), the weighted Z-test (also called ‘Stouffer’s method) is more power and more precision than does Fisher’s test. It favours symmetric rejection and is less sensitive to a single low p-value, requiring more consistently low p-values to yield a low combined p-value.

Here I used the *combinePValues* function in the [scran: Methods for Single-Cell RNA-Seq Data Analysis](https://rdrr.io/bioc/scran/man/combinePValues.html) package to obtain the combined p-values.

---

### Fisher's exact tests

The R scripts of Fisher's exact test for each repeated study are located in GitHub/DelBay_project/R_scripts/Fisher_exact/Fish_repeat

The R scripts of Fisher’s approach and Z-method are located in GitHub/DelBay_project/R_scripts/Fisher_exact/Fish_combine_p

---

### SNPs and Allele Frequency Files (mafs.gz)

Again I used the regenerated dataset (doMajorMinor 3 + doMaf 1) that includes 1,934,038 SNPs from global Angsd calling. This SNP list is used for indivudal population SNP dataset creation.

---

### Fisher’s approach

|   Test       | fdr < 0.1 | fdr < 0.05| fdr < 0.01|
|--------------|-----------|-----------|-----------|
|REF-CH-SR-HC  |     41    | 10        |      1    |
|REF-CH-NB-HC  |     32    | 16        |      6    |
|SR-REF-COH-ARN|     20    | 0         |      0    | 

|Group compared           | fdr < 0.1  | 
|-------------------------|-------------|
|CH vs. REF & HC vs. SR   |      41     |
|CH vs. REF & HC vs. NB   |      32     |
|Shared                   |      10     |

Manhattan plot for CH vs. REF + HC vs. SR. Red dashed line indicates 10% FDR threshold (No. of outliers = 41)

<img src="https://hzz0024.github.io/images/Fish/manhattan_REF_CH_SR_HC_fish.jpeg" alt="img" width="800"/>

Manhattan plot for CH vs. REF + HC vs. NB. Red dashed line indicates 10% FDR threshold (No. of outliers = 32)

<img src="https://hzz0024.github.io/images/Fish/manhattan_REF_CH_NB_HC_fish.jpeg" alt="img" width="800"/>

Manhattan plot for SR vs. REF + COH vs. ARN. Red dashed line indicates 10% FDR threshold (No. of outliers = 20) 

<img src="https://hzz0024.github.io/images/Fish/manhattan_SR-REF-COH-ARN_fish.jpeg" alt="img" width="800"/>

---

### Z or Stouffer’s approach

|   Test       | fdr < 0.1 | fdr < 0.05| fdr < 0.01|
|--------------|-----------|-----------|-----------|
|REF-CH-SR-HC  |     11    |  4        |      0    |
|REF-CH-NB-HC  |     31    |  8        |      1    |
|SR-REF-COH-ARN|      0    |  0        |      0    | 

|Group compared           | fdr < 0.1  | 
|-------------------------|-------------|
|CH vs. REF & HC vs. SR   |      11     |
|CH vs. REF & HC vs. NB   |      31     |
|Shared                   |       2     |

Manhattan plot for CH vs. REF + HC vs. SR. Red dashed line indicates 10% FDR threshold (No. of outliers = 11)

<img src="https://hzz0024.github.io/images/Fish/manhattan_REF-CH-SR-HC_z.jpeg" alt="img" width="800"/>

Manhattan plot for CH vs. REF + HC vs. NB. Red dashed line indicates 10% FDR threshold (No. of outliers = 31)

<img src="https://hzz0024.github.io/images/Fish/manhattan_REF-CH-NB-HC_z.jpeg" alt="img" width="800"/>

Manhattan plot for SR vs. REF + COH vs. ARN. Red dashed line indicates 10% FDR threshold (No. of outliers = 0) 

<img src="https://hzz0024.github.io/images/Fish/manhattan_SR-REF-COH-ARN_z.jpeg" alt="img" width="800"/>

### Shared SNPs (under FDR < 0.1)

I then checked the SNP lists in each group the tried to find if any outlier enriched in particular genome regions.

- Fisher’s approach

```sh
f_a = make_ID('REF-CH-SR-HC_out_all_fish.txt','REF-CH-SR-HC_')
f_b = make_ID('REF-CH-NB-HC_out_all_fish.txt','REF-CH-NB-HC_')
f_c = make_ID('SR-REF-COH-ARN_out_all_fish.txt','SR-REF-COH-ARN_')

f_a - REF-CH-SR-HC outlier using Fisher’s approach (n=41)
[1] "1_32280434" "1_58674134" "2_856732"   "2_20934830" "2_20940319" "2_26318648" "2_30708505" "2_59868440"
[9] "3_24092285" "3_35945753" "4_16409729" "4_18041038" "4_20681942" "4_34660703" "4_34660756" "5_3932376" 
[17] "5_13147035" "5_14331879" "5_14331978" "5_16551755" "5_16551904" "5_16552716" "5_21404457" "5_21405068"
[25] "5_21405082" "5_28997935" "5_37379685" "5_39887057" "5_43506040" "5_46534680" "5_52873868" "5_53318977"
[33] "7_5060336"  "7_33713495" "7_42579646" "8_28557322" "8_37282952" "8_46071132" "8_52952352" "9_51929934"
[41] "9_52184780"

f_b - REF-CH-NB-HC outlier using Fisher’s approach (n=32)
[1] "1_1595742"   "1_2845646"   "1_29942348"  "1_32280271"  "1_32280434"  "1_34970650"  "1_60393516"  "1_64431121" 
[9] "2_41329026"  "2_59868440"  "2_59868476"  "3_64556328"  "3_64973955"  "4_16409729"  "4_18041038"  "4_18361759" 
[17] "5_12533591"  "5_12625464"  "5_12699252"  "5_13147035"  "5_16551904"  "5_16552716"  "5_18212375"  "5_39887057" 
[25] "5_59795885"  "6_990216"    "7_43057673"  "7_52629606"  "8_37282952"  "8_52952352"  "8_56354987"  "10_11389354"

f_c - SR-REF-COH-ARN outlier using Fisher’s approach (n=20)
[1] "1_21189753" "1_33970557" "2_9884637"  "2_52944243" "3_13794806" "3_61366056" "4_5577858"  "4_15986526"
[9] "4_34572201" "4_58016127" "5_11059836" "5_25993454" "5_25993714" "5_40806736" "5_41202785" "5_45823024"
[17] "5_89367996" "8_14238691" "9_65948298" "9_81015150"

f_d <- intersect(f_a, f_b) # REF-CH-SR-HC & REF-CH-NB-HC shared outlier

f_d
[1] "1_32280434" "2_59868440" "4_16409729" "4_18041038" "5_13147035" "5_16551904" "5_16552716" "5_39887057"
[9] "8_37282952" "8_52952352"

No shared outliers between f_c & f_a or f_c & f_b
```

- Z method

```sh
z_a = make_ID('REF-CH-SR-HC_out_all_z.txt','REF-CH-SR-HC_')
z_b = make_ID('REF-CH-NB-HC_out_all_z.txt','REF-CH-NB-HC_')
z_c = make_ID('SR-REF-COH-ARN_out_all_z.txt','SR-REF-COH-ARN_')

z_a - REF-CH-SR-HC outlier using Z method (n=11)
[1] "1_58674134" "2_26318648" "2_30708505" "3_24092285" "4_16409729" "4_34660756" "5_16552716" "5_21404457"
[9] "5_52873868" "7_42579646" "8_46071132"

z_b - REF-CH-NB-HC outlier using Z method (n=31)
[1] "1_1595742"  "1_1634613"  "1_25834012" "1_29942348" "1_32280271" "1_32280434" "1_34970650" "1_59197311"
[9] "1_60393516" "2_41329026" "3_64556328" "3_64973955" "4_2804616"  "4_18361759" "4_25385293" "4_28716656"
[17] "5_8259087"  "5_12625464" "5_13147035" "5_14653344" "5_16551904" "5_16552716" "5_18028707" "5_18212375"
[25] "5_39887057" "5_59795885" "6_990216"   "8_37282952" "8_46071132" "8_56354987" "8_57708473"

z_c - SR-REF-COH-ARN outlier using Z method (n=0)

z_d <- intersect(z_a, z_b) # REF-CH-SR-HC & REF-CH-NB-HC shared outlier

z_d
[1] "5_16552716" "8_46071132"

No shared outliers between z_c & z_a or z_c & z_b
```

- shared by methods

```sh

m_a - REF-CH-SR-HC shared outlier between Fisher's and Z approaches (n = 11)

m_a <- intersect(z_a, f_a) 
length(m_a)
[1] 11

m_a - REF-CH-SR-HC shared outlier between Fisher's and Z approaches (n = 20)

m_b <- intersect(z_b, f_b) 
length(z_b)
[1] 31
length(f_b)
[1] 32
length(m_b)
[1] 20

m_b
[1] "1_1595742"  "1_29942348" "1_32280271" "1_32280434" "1_34970650" "1_60393516" "2_41329026" "3_64556328"
[9] "3_64973955" "4_18361759" "5_12625464" "5_13147035" "5_16551904" "5_16552716" "5_18212375" "5_39887057"
[17] "5_59795885" "6_990216"   "8_37282952" "8_56354987"
```

- Enrichment analysis

I performed a bootstrap analysis to demonstate how rare peaks are following the method in [Dixon et al. (2015) Genomic determinants of coral heat tolerance across latitudes](https://science.sciencemag.org/content/348/6242/1460). 

An example output from Fisher's test looks like this,

Manhattan plot for CH vs. REF + HC vs. SR. Note that the Red dots indicates outliers < FDR 5% threshold (No. of outliers = 10). Light green bars identify
regions with significant clustering of peaks (FDR < 0.1), dark green bars identify regions with significant clustering of peaks (FDR < 0.05, according to 100,000 bootstrapped replicates).

<img src="https://hzz0024.github.io/images/Fish_boot/hard_REF-CH-SR-HC_bootstrap.jpg" alt="img" width="800"/>

Manhattan plot for SR vs. REF + COH vs. ARN. No outliers < FDR 5% in this case.  

<img src="https://hzz0024.github.io/images/Fish_boot/hard_SR-REF-COH-ARN_bootstrap.jpg" alt="img" width="800"/>

Looks like the bootstrap analysis need further optimization. At least those green bars should not show in the control group (SR vs. REF + COH vs. ARN).

---

### Check SNP related SNPs

Some of the outliers were clustered in particular regions in the genomes. For example,

SNPs 5_16551904 and 5_16552716 are located closely and consistently shown in the shared outliers, and gene associated with these two SNPs are actin-depolymerizing factor 1-like (location = NC_035784.1: 16551175..16553910, https://www.ncbi.nlm.nih.gov/gene/111134891). This gene has been reported as  an actin-binding protein that controls actin assembly and has also been reported in Pacific oyster and Manila clam hemocytes (belong to cell structure and motility, GO:0003779). See [Roberts et al. 2009. Analysis of Genes Isolated from Plated Hemocytes of the Pacific Oyster, Crassostreas gigas](https://link.springer.com/content/pdf/10.1007/s10126-008-9117-6.pdf) 

SNP 3_64556328 belongs to the [gene transmembrane protein 184C-like, or tmem184C](https://www.ncbi.nlm.nih.gov/gene/111125354 NC_035782.1), which spans chromosome 3 at 64547915 to 64558382. 




