---
comments: true
title: CVreseq data analyses plan
date: '2020-07-27 12:00'
tags:
  - CVreseq
  - theta
  - ouliter
  - WGS
categories:
  - WGS data analysis
---

Inspired by Sutherland et al. 2020, I'd like to identify the poential outliers under domestication selection using the thinned vcf data. Four methods will be used for data analysis. They are summarized below

| 			 |Key parameters |   Status     |  
| -----------|----------|--------------|
|   pcadapt  | K        |              |
|   BayeScan | pr_odds      |     done     |
|    RDA     | standard deviations|              |
| Percentile | percentile|     done     |
 

---
### Pcadapt

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

- Results

In order to decide if a locus is a good candidate for being under the influence of
selection one need to check the q_value in the *_fst.txt* file

Here I generated six Bayescan outputs for LA, DB_1, and DB_2 comparsions, with pr_odds set as 10 and 100

A python script was developed to filter out snps based on qvalues (FDR alpha value), the usage is shown below,

```python
python3 qvalue_filter.py -i LA_odd100_output_fst.txt -r LA_snplist.txt -o LA_bayescan_outlier.txt -q 0.1

>LA_odd100
>[14] records has q value < 0.1.
>LA_odd10
>[96] records has q value < 0.1.
>DB_1_odd100
>[0] records has q value < 0.1.
>DB_1_odd10
>[49] records has q value < 0.1.
>DB_2_odd100
>[17] records has q value < 0.1.
>DB_2_odd10
>[82] records has q value < 0.1.
```

| Group	     |Populations|   pr_odds    |  qval threshold| No. outlier|
| -----------|-----------|--------------|----------------|------------|
|   LA       | SL-OBOYS2 |    100       |      0.1       |    14      |
|   LA       | SL-OBOYS2 |    10        |      0.1       |    96      |
|   DB_1     | CS-NEH    |    100       |      0.1       |     0      |
|   DB_1     | CS-NEH    |    10        |      0.1       |    49      |
|   DB_2     | CS-DEBY   |    100       |      0.1       |    17      |
|   DB_2     | CS-DEBY   |    10        |      0.1       |    82      |

SL - Louisiana wild line   
OBOYS2 - Louisiana selected line   
NEH - Delaware Bay selected NEH line   
DEBY - Chesapeake Bay selected line (initially from DB)   
CS - Cape Shore (Delaware Bay) wild line 

Among these outliers, only one SNP (5_11774642) was shared between DB_2 (pr_odd 100) and LA (pr_odd 100) results.

- Outlier plots

Mahattan plot for SL-OBOYS2 (pr_odds 10)

<img src="https://hzz0024.github.io/images/CVreseq_bayescan/odd10_LA.jpg" alt="img" width="800"/>

Mahattan plot for SL-OBOYS2 (pr_odds 100)

<img src="https://hzz0024.github.io/images/CVreseq_bayescan/odd100_LA.jpg" alt="img" width="800"/>

Mahattan plot for CS-NEH (pr_odds 10)

<img src="https://hzz0024.github.io/images/CVreseq_bayescan/odd10_DB_1.jpg" alt="img" width="800"/>

Mahattan plot for CS-NEH (pr_odds 100)

<img src="https://hzz0024.github.io/images/CVreseq_bayescan/odd100_DB_1.jpg" alt="img" width="800"/>

Mahattan plot for CS-DEBY (pr_odds 10)

<img src="https://hzz0024.github.io/images/CVreseq_bayescan/odd10_DB_2.jpg" alt="img" width="800"/>

Mahattan plot for CS-DEBY (pr_odds 100)

<img src="https://hzz0024.github.io/images/CVreseq_bayescan/odd100_DB_2.jpg" alt="img" width="800"/>

---
### RDA

---
### Percentile

- get polymorphic loci for each comparsion

```sh
vcftools --vcf LA.recode.vcf --maf 0.05 --recode --recode-INFO-all --out LA_maf
>After filtering, kept 284587 out of a possible 334011 Sites

vcftools --vcf DB_1.recode.vcf --maf 0.05 --recode --recode-INFO-all --out DB_1_maf
>After filtering, kept 282067 out of a possible 334011 Sites

vcftools --vcf DB_2.recode.vcf --maf 0.05 --recode --recode-INFO-all --out DB_2_maf
>After filtering, kept 293946 out of a possible 334011 Sites
```

- estimate the Fst values using VCFtools 

```sh
./vcftools --vcf vcf_file1.vcf --weir-fst-pop individual_list_1.txt --weir-fst-pop individual_list_2.txt

vcftools --vcf LA_maf.recode.vcf --weir-fst-pop SL --weir-fst-pop OBOYS2
Weir and Cockerham mean Fst estimate: 0.024711
Weir and Cockerham weighted Fst estimate: 0.036727
After filtering, kept 284587 out of a possible 284587 Sites

vcftools --vcf DB_1_maf.recode.vcf --weir-fst-pop CS --weir-fst-pop NEH
Weir and Cockerham mean Fst estimate: 0.048347
Weir and Cockerham weighted Fst estimate: 0.063669
After filtering, kept 282067 out of a possible 282067 Sites

vcftools --vcf DB_2_maf.recode.vcf --weir-fst-pop CS --weir-fst-pop DEBY
Weir and Cockerham mean Fst estimate: 0.02112
Weir and Cockerham weighted Fst estimate: 0.033042
After filtering, kept 293946 out of a possible 293946 Sites
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

- Detecting runs of homozygosity (ROH)

I am still working on the ROH analysis using the vcf file and follow the instructure from bcftools site [https://samtools.github.io/bcftools/howtos/roh-calling.html](https://samtools.github.io/bcftools/howtos/roh-calling.html)




