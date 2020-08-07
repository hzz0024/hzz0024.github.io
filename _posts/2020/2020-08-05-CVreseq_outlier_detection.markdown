---
comments: true
title: Examine the potential CVreseq outliers
date: '2020-08-05 12:00'
tags:
  - CVreseq
  - Fst
  - ouliter
  - WGS
categories:
  - WGS data analysis
---

| 			 |Key parameters |
| -----------|---------------|
|   pcadapt  | k             |
|   BayeScan | pr_odds       |   

---
### Pcadapt

To evaluate the LD effects for the outlier identification, I made some plots to show the loadings (contributions of each SNP to the PC) and to evaluate if the loadings are clustered in a single or several genomic regions. 

Using populations from LA as an example: LA-LSSL (SL) vs. LA-OBOY (OBOYS2) Louisiana wild vs. selected line

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

It seems the distribution of the loadings is evenly distributed among these plots, the loading pattern in r^2=0.2 looks good because it 1) captures most of the potential outliers with heavy PC loading values; 2) does not show the "dashed line" patterns. The snps that contribute to this "dashed line" are also outliers to me, but due to the small sample size (6) they may share the same loading values. The main purpose of this post is to examine some of the most significant outlier across the genome, therefore I excluded those snps now; 3) output an reasonalable number of outliers

I then use the r^2=0.2 as a parameter for other data LD pruning. 

| Group	     |Populations| 
| -----------|-----------|
|   LA       | SL-OBOYS2 |
|   DB_1     | CS-NEH    |
|   DB_2     | CS-DEBY   |

SL - Louisiana wild line   
OBOYS2 - Louisiana selected line   
NEH - Delaware Bay selected NEH line    
DEBY - Chesapeake Bay selected line (initially from DB)    
CS - Cape Shore (Delaware Bay) wild line 

LA pcadapt result (666 outliers): 

<img src="https://hzz0024.github.io/images/outlier/LA.jpg" alt="img" width="800"/>

It is interesting that an potential individual outlier is identified from the PCA result, need to remove this one from further analysis (not done yet)

DB_1 pcadapt result (637 outliers): 

<img src="https://hzz0024.github.io/images/outlier/DB_1.jpg" alt="img" width="800"/>

DB_2 pcadapt result (593 outliers): 

<img src="https://hzz0024.github.io/images/outlier/DB_2.jpg" alt="img" width="800"/>

below are results for some common shared outliers (DB_1_list, DB_2_list and LA_list are outlier list in the results)

```sh
> intersect(DB_1_list,DB_2_list)
[1] "1_55485508" "2_2405242"  "2_59641818" "5_49795690" "7_53403"    "8_36715103" "9_50117923"
> intersect(DB_1_list,LA_list)
[1] "2_656131"   "3_54596543" "8_63266137" "9_12057031"
> intersect(DB_2_list,LA_list)
[1] "3_69917769" "5_25337621" "5_79382471" "8_10601236" "9_52331392"
```
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
# intersect is the function in r to check common shared values in two or more vectors

# first compare the methods but the same population pair
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
```

Compare the results between populations from pcadapt results,

```sh

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

See if there is any SNP outliers shared by population pairs (from pcadapt results none are shared by Bayescan methods) and among the methods,

```sh
> intersect(intersect(DB_1_list,DB_2_list), intersect(DB_1_list,Bas_DB_1$V1))
[1] "8_36715103"

> intersect(intersect(DB_1_list,LA_list), intersect(DB_1_list,Bas_DB_1$V1))
[1] "2_656131"   "8_63266137"
```

---
### Examine the pi patterns using potential outliers

For pi estimate in each population, I use vcftools with code below,

```sh
# for single snp diversity
for pop in SL OBOYS2 CS NEH DEBY; do
        vcftools --vcf $pop'_sort.vcf' --site-pi --out $pop
done
# diversity in windows
for pop in SL OBOYS2 CS NEH DEBY; do
    for win in 100 200 500 1000 5000; do
        vcftools --vcf $pop'_sort.vcf' --window-pi $win --out $pop'_'$win
    done
done
```

| Group	     |Populations| 
| -----------|-----------|
|   LA       | SL-OBOYS2 |
|   DB_1     | CS-NEH    |
|   DB_2     | CS-DEBY   |

SL - Louisiana wild line   
OBOYS2 - Louisiana selected line   
NEH - Delaware Bay selected NEH line    
DEBY - Chesapeake Bay selected line (initially from DB)    
CS - Cape Shore (Delaware Bay) wild line

- pi distribution

SNP 8_36715103 is shared by two methods (pcadatp and Bayescan) in CS-NEH comparsion, and shared by CS-NEH and CS-DEBY using pcadatp. The dashed line indicates the location of this snp. From the figure CS-NEH (DB_1), a lot of pi differences are shown around this outlier position. The wild CS population has increased pi at 36.5-36.8M genome region. However, in CS-DEBY I did not see this pattern, although this snp has also been identified as an outlier.

<img src="https://hzz0024.github.io/images/outlier/8_36715103_pi_DB_1.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/outlier/8_36715103_pi_DB_2.jpg" alt="img" width="800"/>

SNP 8_63266137 is shared by two methods (pcadatp and Bayescan) in CS-NEH comparsion, and shared by CS-NEH and SL-OBOYS2 using pcadatp. Here I zoom in (200k bp) to take a closer look at this loci due to the lack of values around this region. This loci looks like a false-positive as the observed pi values in SL-OBOYS2 (LA) is opposite from the expection, potentially due the individual outlier in the LA populaiton pair.

<img src="https://hzz0024.github.io/images/outlier/8_63266137_DB_1.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/outlier/8_63266137_LA.jpg" alt="img" width="800"/>

SNP 5_36361325 is only detected in CS-NEH comparsion, but SNPs in chromosome 5 are constantly identified as outliers. I'd like to check the pi patterns for regions around this snp. It looks like this outlier is population pair specific. I only observed elevated pi pattern in DB_1 population pair.

<img src="https://hzz0024.github.io/images/outlier/5_36361325_DB_1.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/outlier/5_36361325_DB_2.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/outlier/5_36361325_LA.jpg" alt="img" width="800"/>

- Tajima'd

Currently I am not able to plot the Tajima'd and pi in the same plot because the length of data are not identical due to some "nan" values (particular for Tajima's D). I need to figure out the way of making dual y-axis plots for both values. 

Here I simplyly use SNP 8_36715103 as an example (values below are esimated using 1000 bp as window size) to show the results in a population pair (DB_1 here),

Pi comparison

> CS population 1000.windowed.pi

|CHR  | Window_start | Window_end |NO. SNPs | pi      | 
|-----|--------------|------------|---------|---------|  
|8	  |36714001	     |36715000	  | 11	    |0.0040303|   
|8	  |36715001	     |36716000	  |11	    |0.0041818|     
|8	  |36717001      |36718000	  |4	    |0.0014090| 

> NEH population 1000.windowed.pi

|CHR  | Window_start | Window_end |NO. SNPs | pi      | 
|-----|--------------|------------|---------|---------|  
|8	  |36714001	     |36715000	  | 4 	    |0.00131818|   
|8	  |36715001	     |36716000	  | 1	    |0.00030303|     
|8	  |36717001      |36718000	  | 1	    |0.00030303| 

Tajima.D comparison

> CS population 1000.windowed.tajima's D

|CHR  | Window_start | Window_end |NO. SNPs | Tajima'D| 
|-----|--------------|------------|---------|---------|  
|8	  |36714001	     |36715000	  | 11	    |0.442939 |   
|8	  |36715001	     |36716000	  |11	    |0.61601  |     
|8	  |36717001      |36718000	  |4	    |na       |


> NEH population 1000.windowed.tajima's D

|CHR  | Window_start | Window_end |NO. SNPs | Tajima'D| 
|-----|--------------|------------|---------|---------|  
|8	  |36714001	     |36715000	  | 4	    |-0.016928|   
|8	  |36715001	     |36716000	  |1	    |-0.194921|     
|8	  |36717001      |36718000	  |0	    |na       |


 Here I observed Tajima's D in NEH domesticated population to be negative due to positive selection (or selective sweep), and positive in wild CS population. Need to check some other outilers.

