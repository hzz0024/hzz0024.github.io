---
comments: true
title: CVreseq theta & Tajima's D estimates 
date: '2020-03-29 12:00'
tags:
  - angsd
  - SFS
  - theta
  - pi
  - WGS
categories:
  - WGS data analysis
---

This post lists the script details of theta estimation using the angsd pipelines. The origianl command lines of theta calculation is [here](http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests). Given that we do not have ancestral reference genome, the folded sfs will be used for the calculation. 

We have three datasets, each with 12 individuals. They are:

- LA-LSSL vs. LA-OBOY

Louisiana wild vs. selected line, source of oysters from selection from approx same environment (Oyster Bayou)

- DB-HSCS vs. DB-NEHD

Wild (Cape Shore, Delaware Bay wild high salinity) vs. selected (core line NEH)

- DB-HSCS vs. CB-DEBY

Wild (Cape Shore, Delaware Bay wild high salinity) vs. Chesapeake Bay selected (initially from DB)

Advice from angsd authors and example command:

> If you don't have the ancestral states, you can still calculate the Watterson and Tajima theta, which means you can perform the Tajima's D neutrality test statistic. But this requires you to use the folded sfs. The output files will have the same format, but only the thetaW and thetaD, and tajimas D is meaningful. Below is an example based on the earlier example where we now base our analysis on the folded spectrum. Notice the -fold 1

> First estimate the folded site allele frequency likelihood
```sh
./angsd -bam bam.filelist -doSaf 1 -anc hg19.fa -GL 1 -P 24 -out outFold -fold 1
```

> Obtain the maximum likelihood estimate of the SFS
```sh
misc/realSFS outFold.saf.idx -P 24 > outFold.sfs
```

> Calculate the thetas (remember to fold)
```sh
./angsd -bam bam.filelist -out outFold -doThetas 1 -doSaf 1 -pest outFold.sfs -anc hg19.fa -GL 1 -fold 1
```

> Estimate Tajimas D
```sh
thetaStat do_stat outFold.thetas.idx -nChr 10 
```

### CVreseq dataset

- saf file generation

```sh
module load angsd/0.931
angsd -b LA.bamlist -anc CVtot_3.0_newchr.fa.bgz -out 1_saf/LA_no12inv -dosaf 1 -doDepth 1 -doCounts 1 -GL 1 -P 16 -fold 1
angsd -b DB.bamlist -anc CVtot_3.0_newchr.fa.bgz -out 1_saf/DB_no12inv -dosaf 1 -doDepth 1 -doCounts 1 -GL 1 -P 16 -fold 1
angsd -b CB.bamlist -anc CVtot_3.0_newchr.fa.bgz -out 1_saf/DB_no12inv -dosaf 1 -doDepth 1 -doCounts 1 -GL 1 -P 16 -fold 1
```
- sfs, theta and Tajima's D

```sh
## Get SFS from saf 
module load angsd/0.931

BASE_DIR='/scratch/hzz0024/CVreseq/2_theta/'
for POP in {DB,CB,LA}; do
    echo $POP' saf starts'
    /tools/angsd-0.931/misc/realSFS \
      $BASE_DIR$POP'_no12inv.saf.idx' \
      -P 12 \
      > $BASE_DIR$POP'_no12inv.sfs'
done
## Estimate theta
BASE_DIR='/scratch/hzz0024/CVreseq/2_theta/'
for POP in {DB,CB,LA}; do
    echo $POP' saf theta estimate'
    /tools/angsd-0.931/angsd \
    -bam $BASE_DIR$POP'.bamlist' \
    -out $BASE_DIR$POP'_no12inv' \
    -doThetas 1 \
    -doSaf 1 \
    -fold 1 \
    -pest $BASE_DIR$POP'_no12inv.sfs' \
    -anc CVtot_3.0_newchr.fa.bgz \
    -GL 1 \
    -P 12
done

## Print per-SNP theta
BASE_DIR='/scratch/hzz0024/CVreseq/2_theta/'
for POP in {DB,CB,LA}; do
    /tools/angsd-0.931/misc/thetaStat print \
    $BASE_DIR$POP'_no12inv.thetas.idx' \
    > $BASE_DIR$POP'_no12inv.thetas.tsv'
done

## Do fixed window theta
BASE_DIR='/scratch/hzz0024/CVreseq/2_theta/'
for POP in {DB,CB,LA}; do
    /tools/angsd-0.931/misc/thetaStat do_stat \
    $BASE_DIR$POP'_no12inv.thetas.idx' \
    -win 10000 -step 10000 \
    -outnames $BASE_DIR$POP'_no12inv.thetas.idx' \
    > $BASE_DIR$POP'_no12inv.windowed_thetas.log'
done

## Do per-chromosome average theta
BASE_DIR='/scratch/hzz0024/CVreseq/2_theta/'
for POP in {DB,CB,LA}; do
    /tools/angsd-0.931/misc/thetaStat do_stat \
    $BASE_DIR$POP'_no12inv.thetas.idx' \
    -outnames $BASE_DIR$POP'_no12inv.average_thetas.idx' \
    > BASE_DIR$POP'_no12inv.average_thetas.log'
done
```
---
#### Results

- Output in the ./thetaStat print thetas.idx are the log scaled per site estimates of the thetas     

```sh
# script used to check .theta.idx
/tools/angsd-0.931/misc/thetaStat print LA_no12inv.thetas.idx
```
- Output in the pestPG file are the sum of the per site estimates for a region   

> The .pestPG file is a 14 column file (tab seperated). The first column contains information about the region. The second and third column is the reference name and the center of the window. We then have 5 different estimators of theta, these are: Watterson, pairwise, FuLi, fayH, L. And we have 5 different neutrality test statistics: Tajima's D, Fu&Li F's, Fu&Li's D, Fay's H, Zeng's E. The final column is the effetive number of sites with data in the window. Most likely you are just interest in the wincenter (column 3) and the column 9 which is the Tajima's D statistic.



 





