---
comments: true
title: Percentile outlier Mantel test 
date: '2020-07-13 12:00'
tags:
  - DelBay19
  - ouliter
  - Mantel test
  - WGS
categories:
  - WGS data analysis
---

This post shows the steps for Mantel test. The Fst datasets were created from this post [Rebuild DelBay19 data](https://hzz0024.github.io/2020/07/08/DelBay_data_redo.html)

Given that LD is minimal beyond 200 bp, I used 200bp as the window size to produce the Fst values. The wild transect has 5 population, leading to 10 pairwise comparsions. Below are the steps for data prepariation and Mantel test:

1) Replace the chromosome name (string) with numbers

```sh
for i in *.txt; do
sed -i.bak 's/NC_035780.1/1/g;s/NC_035781.1/2/g;s/NC_035782.1/3/g;s/NC_035783.1/4/g;s/NC_035784.1/5/g;s/NC_035785.1/6/g;s/NC_035786.1/7/g;s/NC_035787.1/8/g;s/NC_035788.1/9/g;s/NC_035789.1/10/g;s/NC_007175.2/11/g' $i
done

```

2) extract the chr, window_start, window_end, and Fst values for each SNP

```sh
awk -F"\t" '{print $2, $3, $3, $5}' CH_REF_maf0.05_pctind0.7_cv30_nochr56invers_fold.200_win_200_fst.txt | awk '$2-=100' | awk '$3+=100' > ch_ref_200_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' ARN_HC_maf0.05_pctind0.7_cv30_no56invers_fold.200_win_200_fst.txt | awk '$2-=100' | awk '$3+=100' > HC_ARN_200_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' COH_HC_maf0.05_pctind0.7_cv30_no56invers_fold.200_win_200_fst.txt | awk '$2-=100' | awk '$3+=100' > HC_COH_200_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' HC_SR_maf0.05_pctind0.7_cv30_no56invers_fold.200_win_200_fst.txt  | awk '$2-=100' | awk '$3+=100' > HC_SR_200_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' HC_NB_maf0.05_pctind0.7_cv30_no56invers_fold.200_win_200_fst.txt  | awk '$2-=100' | awk '$3+=100' > HC_NB_200_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' ARN_COH_maf0.05_pctind0.7_cv30_no56invers_fold.200_win_200_fst.txt| awk '$2-=100' | awk '$3+=100' > ARN_COH_200_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' ARN_SR_maf0.05_pctind0.7_cv30_no56invers_fold.200_win_200_fst.txt | awk '$2-=100' | awk '$3+=100' > ARN_SR_200_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' ARN_NB_maf0.05_pctind0.7_cv30_no56invers_fold.200_win_200_fst.txt | awk '$2-=100' | awk '$3+=100' > ARN_NB_200_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' COH_SR_maf0.05_pctind0.7_cv30_no56invers_fold.200_win_200_fst.txt | awk '$2-=100' | awk '$3+=100' > COH_SR_200_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' COH_NB_maf0.05_pctind0.7_cv30_no56invers_fold.200_win_200_fst.txt | awk '$2-=100' | awk '$3+=100' > COH_NB_200_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' NB_SR_maf0.05_pctind0.7_cv30_no56invers_fold.200_win_200_fst.txt  | awk '$2-=100' | awk '$3+=100' > SR_NB_200_fold.fst

```

3) add the header for each fst file

```sh
for i in *.fst;do
    python3 3_edit_line.py $i $i
done 
```

4) extract the percentile data (99.9% percentile)

```sh
python3 4_percentile.py -i ch_ref_200_fold.fst -o ch_ref_200.csv -p 99.9 > ch_ref_200.log
python3 4_percentile.py -i HC_ARN_200_fold.fst -o HC_ARN_200.csv -p 99.9 > HC_ARN_200.log
python3 4_percentile.py -i HC_COH_200_fold.fst -o HC_COH_200.csv -p 99.9 > HC_COH_200.log
python3 4_percentile.py -i HC_SR_200_fold.fst -o HC_SR_200.csv -p 99.9 > HC_SR_200.log
python3 4_percentile.py -i HC_NB_200_fold.fst -o HC_NB_200.csv -p 99.9 > HC_NB_200.log
python3 4_percentile.py -i ARN_COH_200_fold.fst -o ARN_COH_200.csv -p 99.9 > ARN_COH_200.log
python3 4_percentile.py -i ARN_SR_200_fold.fst -o ARN_SR_200.csv -p 99.9 > ARN_SR_200.log
python3 4_percentile.py -i ARN_NB_200_fold.fst -o ARN_NB_200.csv -p 99.9 > ARN_NB_200.log
python3 4_percentile.py -i COH_SR_200_fold.fst -o COH_SR_200.csv -p 99.9 > COH_SR_200.log
python3 4_percentile.py -i COH_NB_200_fold.fst -o COH_NB_200.csv -p 99.9 > COH_NB_200.log
python3 4_percentile.py -i SR_NB_200_fold.fst -o SR_NB_200.csv -p 99.9 > SR_NB_200.log

```


5) rouguly calculate the geographic distance between two sites. The approximate distance was calculated on the basis of a spherical earch and used the website here [https://www.movable-type.co.uk/scripts/latlong.html](https://www.movable-type.co.uk/scripts/latlong.html).

|ID  |pop |coordinate 1|coordinate 2|
|----|----|-----------|-------------|
|1	 |HC  |39.45084	  |-75.518749   |
|2	 |ARN |	39.387956 |	-75.453399  |
|3	 |COH |	39.325903 |	-75.374292  |
|4	 |SR  |	39.285095 |	-75.327286  |
|5	 |NB  |	39.26557  |	-75.277987  |


6) estimate average Fst across all 200bp windowed bins for each pairwise comparison, see results below

|pop1|pop2|average Fst|Fst/(1-Fst)|Distance(km)|
|----|----|-----------|-----------|------------|
|HC	 |ARN |0.058069561|0.061649522| 8.967      |
|HC	 |COH |0.061133346|0.065113982|18.63       |
|HC	 |SR  |0.05861558 |0.062265297|24.71       | 
|HC	 |NB  |0.058487175|0.062120423|29.2        | 
|ARN |COH |0.060670008|0.064588599|9.689       | 
|ARN |SR  |0.056563392|0.05995463 |15.76       | 
|ARN |NB  |0.057787892|0.061332148|20.32       |  
|COH |SR  |0.063892914|0.068253851|6.078       | 
|COH |NB  |0.064080804|0.068468308|10.66       | 
|SR	 |NB  |0.059037481|0.062741586|4.767       | 


7) create the matrix input for IBD

```sh
GENETIC_DISTANCE
1 2 0.061649522
1 3 0.065113982
1 4 0.062265297
1 5 0.062120423
2 3 0.064588599
2 4 0.05995463
2 5 0.061332148
3 4 0.068253851
3 5 0.068468308
4 5 0.062741586
GEOGRAPHIC_DISTANCE
1 2 8.967
1 3 18.63
1 4 24.71
1 5 29.2
2 3 9.689
2 4 15.76
2 5 20.32
3 4 6.078
3 5 10.66
4 5 4.767
```

8) run the IBD (only windows version is working)

the outcome and log file are stored [here](https://docs.google.com/spreadsheets/d/1xcbEXtuEYNLG2BUW8Ivh0o79iArKyGI7mixq_wzaAlc/edit?usp=sharing)

9) plot the genetic vs. geographic distance

<img src="https://hzz0024.github.io/images/Mantel/mantel1.jpeg" alt="img" width="800"/>

