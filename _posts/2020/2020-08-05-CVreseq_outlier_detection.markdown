---
comments: true
title: CVreseq outlier detection
date: '2020-08-05 12:00'
tags:
  - CVreseq
  - Fst
  - ouliter
  - WGS
categories:
  - WGS data analysis
---

Inspired by Sutherland et al. 2020, I'd like to identify the poential outliers under domestication selection using the thinned vcf data. Four methods will be used for data analysis. They are summarized below

| 			 |Key parameters |   Status     |  
| -----------|---------------|--------------|
|   pcadapt  | k             |     done     |
|   BayeScan | pr_odds       |     done     |
| Percentile | percentile    |     done     |
| outFlank   | k             |     done     |


---
### Pcadapt

To evaluate if LD might be an issue for the dataset, I made some plots to show the loadings (contributions of each SNP to the PC) and to evaluate if the loadings are clustered in a single or several genomic regions. 

Population: LA-LSSL (SL) vs. LA-OBOY (OBOYS2) Louisiana wild vs. selected line

SNP number after LD prunning: 14255; r^2: 0.1

<img src="https://hzz0024.github.io/images/pcadapt/LA_r0.1_14255.jpeg" alt="img" width="800"/>

SNP number after LD prunning: 27767; r^2: 0.2

<img src="https://hzz0024.github.io/images/pcadapt/LA_r0.2_27767.jpeg" alt="img" width="800"/>

SNP number after LD prunning: 56893; r^2: 0.3

<img src="https://hzz0024.github.io/images/pcadapt/LA_r0.3_56893.jpeg" alt="img" width="800"/>

SNP number after LD prunning: 96993; r^2: 0.4

<img src="https://hzz0024.github.io/images/pcadapt/LA_r0.4_98693.jpeg" alt="img" width="800"/>

SNP number after LD prunning: 159890; r^2: 0.5

<img src="https://hzz0024.github.io/images/pcadapt/LA_r0.5_159890.jpeg" alt="img" width="800"/>

SNP number after LD prunning: 213149; r^2: 0.6

<img src="https://hzz0024.github.io/images/pcadapt/LA_r0.6_213149.jpeg" alt="img" width="800"/>

SNP number after LD prunning: 251625; r^2: 0.7

<img src="https://hzz0024.github.io/images/pcadapt/LA_r0.7_251625.jpeg" alt="img" width="800"/>

SNP number after LD prunning: 264918; r^2: 0.8

<img src="https://hzz0024.github.io/images/pcadapt/LA_r0.8_264918.jpeg" alt="img" width="800"/>

The distribution of the loadings is evenly distributed in all plots, the loading pattern in r^2=0.2 looks good to me because it 1) largely captures the potential outliers based on PC loading values; 2) does not show the "dashed line" patter (posibilly due the linked SNPs); 3) output an reasonalable number of outliers

I then use the r^2=0.2 as a parameter for other data LD pruning. 

we can have a look at the genome scan, which correctly identifies regions involved in adaptation.

> intersect(DB_1_list,DB_2_list)
[1] "1_55485508" "2_2405242"  "2_59641818" "5_49795690" "7_53403"    "8_36715103" "9_50117923"
> intersect(DB_1_list,LA_list)
[1] "2_656131"   "3_54596543" "8_63266137" "9_12057031"
> intersect(DB_2_list,LA_list)
[1] "3_69917769" "5_25337621" "5_79382471" "8_10601236" "9_52331392"

---
### Bayescan

- Parameter settings:  
-n 5000 Number of outputted iterations   
-nbp 20 Number of pilot runs   
-pilot 5000 Length of pilot runs   
-burn 50000 Burn-in length   
-pr_odds 10 or 100 Prior odds for the neutral model   

The key parameter of Bayescan running is the pr_odds, which indicates how much more likely we think the neutral model is compared to the model with selection. For example a pr_odds = 10 indicates that we think the neutral model is 10 times more likely than the model with selection. With low pr_odds settings, the number of false positives is likely to be large unless we use a very stringent FDR
threshold (see detailed parameter explaination and tutorial [here](http://evomicsorg.wpengine.netdna-cdn.com/wp-content/uploads/2016/01/BayeScan_BayeScEnv_exercises.pdf))

- Data processing and running

```sh
# still using the LD prunned snp list _tmp.prune.in
vcftools --vcf DB_1.sort.id.vcf --snps DB_1_tmp.prune.in --recode --out DB_1_prune
# kept 24112 out of a possible 334011 Sites
vcftools --vcf DB_2.sort.id.vcf --snps DB_2_tmp.prune.in --recode --out DB_2_prune
# kept 28320 out of a possible 334011 Sites
vcftools --vcf LA.sort.id.vcf --snps LA_tmp.prune.in --recode --out LA_prune
# kept 27767 out of a possible 334011 Sites

./vcf2genepop_hg.pl vcf=LA_prune.recode.vcf pops=LA,OBOYS2 > LA.gen
./vcf2genepop_hg.pl vcf=DB_1_prune.recode.vcf pops=CS,NEH > DB_1.gen
./vcf2genepop_hg.pl vcf=DB_2_prune.recode.vcf pops=CS,DEBY > DB_2.gen

# next using PGSspider to convert the genepop to bayescan format

# load the bayescan format files to the cluster and run 

./BayeScan2.1_linux64bits LA_odds1 -o LA_odds1_output -pr_odds 1
./BayeScan2.1_linux64bits DB_1_odds1 -o DB_1_odds1_output -pr_odds 1
./BayeScan2.1_linux64bits DB_2_odds1 -o DB_2_odds1_output -pr_odds 1

- Results

In order to decide if a locus is a good candidate for being under the influence of
selection one need to check the q_value in the *_fst.txt* file

Here I generated six Bayescan outputs for LA, DB_1, and DB_2 comparsions, with pr_odds set as 10 and 100

A python script was developed to filter out snps based on qvalues (FDR alpha value), the usage is shown below,

```sh
python3 qvalue_filter.py -i LA_odds10_output_fst.txt -r LA_snplist.txt -o LA_bayescan_outlier.txt -q 0.1
>[22] records has q value < 0.1.

python3 qvalue_filter.py -i DB_1_odds10_output_fst.txt -r DB_1_snplist.txt -o DB_1_bayescan_outlier.txt -q 0.1
>[58] records has q value < 0.1.

python3 qvalue_filter.py -i DB_2_odds10_output_fst.txt -r DB_2_snplist.txt -o DB_2_bayescan_outlier.txt -q 0.1
> [27] records has q value < 0.1.
```

| Group	     |Populations|   pr_odds    |  qval threshold| No. outlier|
| -----------|-----------|--------------|----------------|------------|
|   LA       | SL-OBOYS2 |    10        |      0.1       |    22      |
|   DB_1     | CS-NEH    |    10        |      0.1       |    58      |
|   DB_2     | CS-DEBY   |    10        |      0.1       |    27      |

SL - Louisiana wild line   
OBOYS2 - Louisiana selected line   
NEH - Delaware Bay selected NEH line   
DEBY - Chesapeake Bay selected line (initially from DB)   
CS - Cape Shore (Delaware Bay) wild line 

Here I'd like to compare the common shared outliers between two methods above. 

```sh
# DB_1_list is the result from pcadapt, Bas_DB_1$V1 is the result from bayescan. Same for the other populations.
> intersect(DB_1_list,Bas_DB_1$V1)
 [1] "1_30849333" "1_32752499" "1_58534084" "1_60259825" "2_656131"   "2_5071524"  "2_21500266" "2_45981261" "2_52157647" "3_40979403" "3_41730028" "3_70055986"
[13] "4_15417342" "4_18695731" "4_18722099" "4_31415808" "4_34685619" "4_55313638" "4_57208049" "5_16630602" "5_17325109" "5_21286358" "5_27887499" "5_36361325"
[25] "5_41952828" "5_43217405" "5_53541677" "5_62075320" "5_65420222" "6_4960897"  "6_18975650" "6_35317580" "6_35394384" "6_42853582" "7_17539134" "7_26413160"
[37] "7_49024371" "7_56713582" "8_28227013" "8_36715103" "8_57605961" "8_63266137"
length(intersect(DB_1_list,Bas_DB_1$V1))
[1] 42
> intersect(DB_2_list,Bas_DB_2$V1)
 [1] "2_9316508"  "2_21494062" "2_59711343" "3_20546749" "4_8027221"  "4_9429353"  "4_34547560" "5_19572926" "5_32632886" "5_48031603" "5_54049060" "5_65124487"
[13] "6_7992530"  "6_14666712" "7_7204094"  "7_15410073" "8_15737746" "8_19629404" "8_34017649" "8_43295276" "9_22343202" "9_39026257" "9_81096281" "10_4996789"
> length(intersect(DB_2_list,Bas_DB_2$V1))
[1] 24
> intersect(LA_list,Bas_LA$V1)
 [1] "1_15977907" "2_31084149" "4_10967020" "5_32604685" "5_40191744" "5_93475688" "8_42853491" "9_22044992" "10_1206118" "10_7829874"
> length(intersect(LA_list,Bas_LA$V1))
[1] 10
# compare the results from pcadapt
> intersect(DB_1_list,DB_2_list)
[1] "1_55485508" "2_2405242"  "2_59641818" "5_49795690" "7_53403"    "8_36715103" "9_50117923"
> length(intersect(DB_1_list,DB_2_list))
[1] 7
> intersect(DB_1_list,LA_list)
[1] "2_656131"   "3_54596543" "8_63266137" "9_12057031"
> length(intersect(DB_1_list,LA_list))
[1] 4
> intersect(DB_2_list,LA_list)
[1] "3_69917769" "5_25337621" "5_79382471" "8_10601236" "9_52331392"
> length(intersect(DB_2_list,LA_list))
[1] 5
```

See if there is any SNP outliers shared by two populations (from pcadapt results as none are shared by Bayescan methods) and between the methods,

```sh
> intersect(intersect(DB_1_list,DB_2_list), intersect(DB_1_list,Bas_DB_1$V1))
[1] "8_36715103"

> intersect(intersect(DB_1_list,LA_list), intersect(DB_1_list,Bas_DB_1$V1))
[1] "2_656131"   "8_63266137"
```

see vcftools [manual](https://vcftools.github.io/man_latest.html) and biostars [post](https://www.biostars.org/p/46858/) for detailed parameter explanation. Here the Fst values are estimated using Weir and Cockerham’s (1984) method

- calculate the percentile and extract the SNP information

A python script is made for this purpose, for example

```python
python3 percentile.py -i input -o output -p target_percentile
```

| Group	     |Populations|  Percentile  |  Fst threshold | No. outlier|
| -----------|-----------|--------------|----------------|------------|
|   LA       | SL-OBOYS2 |    99.9      |      0.70      |    269     |
|   DB_1     | CS-NEH    |    99.9      |      0.73      |    260     |
|   DB_2     | CS-DEBY   |    99.9      |      0.70      |    239     |

Nubmber of shared outliers

| Group	     |Fst ranges    | No. of shared outlier|
| -----------|--------------|----------------------|
|   LA-DB_1  |    0.8       |    2                 |
|   LA-DB_2  |    0.8       |    1                 |
| DB_1-DB_2  |    0.73-1    |    7                 |
|DB_1-DB_2-LA|    0.8       |    1                 |

Again, SNP 5_11774642 is shared among three groups.

---
### Conclusion and questions

1. In general, outliers were found throughout the genome. Top outliers with the strongest evidence and consistent observations included SNP 5_11774642. 
 
2. The genomic location including containing or nearby these outliers (within 10 kb) need further investigation.

3. Due to the limited sample size (6 for each populations), the Fst may not reflect the real genomic differentiation.




