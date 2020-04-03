---
comments: true
title: CVreseq theta update 
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

This post lists the script updates of theta estimation using the angsd pipelines. The origianl command lines of theta calculation is [here](http://www.popgen.dk/angsd/index.php/Thetas,Tajima,Neutrality_tests). Given that we do not have ancestral reference genome, the folded sfs will be used for the calculation. 

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

SURPRISINGLY, the angsd author changed the webpage of "Thetas,Tajima,Neutrality_tests" this week, and removed the command lines above. The current code is shown below,

> First estimate the folded site allele frequency likelihood
```sh
./angsd -bam bam.filelist -doSaf 1 -anc chimpHg19.fa -GL 1 -P 24 -out out
```
> Obtain the maximum likelihood estimate of the SFS (remember to fold)
```sh
./misc/realSFS out.saf.idx -P 24 > out.sfs
# use -fold 1 in the above command if you dont have ancestral state.
```
> Calculate the thetas (remember to fold), this is the major step that has been changed!
```sh
./misc/realSFS saf2theta out.saf.idx -outname out -sfs out.sfs
```
> Estimate thetas for every Chromosome/scaffold
```sh
./misc/thetaStat do_stat out.thetas.idx
```
> Estimate Tajimas D and other statistics
```sh
#calculate Tajimas D
./misc/thetaStat do_stat out.thetas.idx
```
> Do a sliding window analysis based on the output from the make_bed command.
```sh
./misc/thetaStat do_stat out.thetas.idx -win 50000 -step 10000  -outnames theta.thetasWindow.gz
```

The output formats and files are still the same. Output in the ./thetaStat print thetas.idx are the log scaled per site estimates of the thetas. Column from left to right:    
1. chromosome   
2. position   
3. ThetaWatterson   
4. ThetaD (nucleotide diversity)   
5. Theta? (singleton category)   
6. ThetaH   
7. ThetaL 

The .pestPG file is a 14 column file (tab seperated). The first column contains information about the region. The second and third column is the reference name and the center of the window. We then have 5 different estimators of theta, these are: Watterson, pairwise, FuLi, fayH, L. And we have 5 different neutrality test statistics: Tajima's D, Fu&Li F's, Fu&Li's D, Fay's H, Zeng's E. The final column is the effetive number of sites with data in the window. Most likely you are just interest in the wincenter (column 3) and the column 9 which is the Tajima's D statistic.

### Data running

I ran the updated theta command lines using the five seperate populations (one population from Cape Shore, Delaware Bay is shared for DB and CB comparsion) in CVreseq, they are

SL - Louisiana wild line
OBOYS2 - Louisiana selected line
NEH - Delaware Bay selected NEH line
DEBY - Chesapeake Bay selected line (initially from DB)
CS - Cape Shore (Delaware Bay) wild line

- saf file generation

```sh
module load angsd/0.931
angsd -b SL.bamlist -anc CVtot_3.0_newchr.fa.bgz -out 1_saf/SL_no6inv -dosaf 1 -doDepth 1 -doCounts 1 -GL 1 -P 16
angsd -b OBOYS2.bamlist -anc CVtot_3.0_newchr.fa.bgz -out 1_saf/OBOYS2_no6inv -dosaf 1 -doDepth 1 -doCounts 1 -GL 1 -P 16
angsd -b NEH.bamlist -anc CVtot_3.0_newchr.fa.bgz -out 1_saf/NEH_no6inv -dosaf 1 -doDepth 1 -doCounts 1 -GL 1 -P 16
angsd -b DEBY.bamlist -anc CVtot_3.0_newchr.fa.bgz -out 1_saf/DEBY_no6inv -dosaf 1 -doDepth 1 -doCounts 1 -GL 1 -P 16
angsd -b CS.bamlist -anc CVtot_3.0_newchr.fa.bgz -out 1_saf/CS_no6inv -dosaf 1 -doDepth 1 -doCounts 1 -GL 1 -P 16
```
- sfs, theta and Tajima's D

```sh
module load angsd/0.931

## Step 1 Get SFS from saf
BASE_DIR='/scratch/hzz0024/CVreseq/2_theta/'
for POP in {OBOYS2,DEBY,CS,SL,NEH}; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 1 '$POP' sfs starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/realSFS \
    $BASE_DIR$POP'_no6inv.saf.idx' \
    -P 12 \
    -fold 1 \
    > $BASE_DIR$POP'_no6inv.sfs'
done

## Step 2a Estimate theta
for POP in {OBOYS2,DEBY,CS,SL,NEH}; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 2a '$POP' theta estimate starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/realSFS saf2theta \
    $BASE_DIR$POP'_no6inv.saf.idx' \
    -sfs $BASE_DIR$POP'_no6inv.sfs' \
    -fold 1 \
    -outname $BASE_DIR$POP'_no6inv'
done

## Step 2b Print per-SNP theta
for POP in {OBOYS2,DEBY,CS,SL,NEH}; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 2a '$POP' per-SNP theta print starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/thetaStat print \
    $BASE_DIR$POP'_no6inv.thetas.idx' \
    > $BASE_DIR$POP'_no6inv.thetas.tsv'
done

## Step 3a do fixed window theta
for POP in {OBOYS2,DEBY,CS,SL,NEH}; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 3a '$POP' fixed window theta starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/thetaStat do_stat \
    $BASE_DIR$POP'_no6inv.thetas.idx' \
    -win 10000 -step 10000 \
    -outnames $BASE_DIR$POP'_no6inv.thetas.window.idx'
done

## Step 3b do per-chromosome average theta
for POP in {OBOYS2,DEBY,CS,SL,NEH}; do
    echo '~~~~~~~~~~~~~~~~~~~~~~~step 3b '$POP' per-chromosome average theta starts~~~~~~~~~~~~~~~~~~~~~~~'
    /tools/angsd-0.931/misc/thetaStat do_stat \
    $BASE_DIR$POP'_no6inv.thetas.idx' \
    -outnames $BASE_DIR$POP'_no6inv.thetas.average.idx'
done
```

I also used angsd command to estimate the loci depth distribution of the 30 CVreseq samples
```sh
module load angsd/0.931
angsd -b all.bamlist -anc CVtot_3.0_newchr.fa.bgz -out QC_CVreseq_maxD2000 -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -minMapQ 1 -minQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 -maxDepth 2000 -P 16 >& QC.log
#plot the qc
Rscript plotQC.R QC_CVreseq_maxD2000 2> /dev/null
```
---
#### Results

Note: for folded results, only thetaW and thetaD, and tajimas D are meaningful.

- Output in the ./thetaStat print thetas.idx are the log scaled per site estimates of the thetas. Column from left to right:    
1. chromosome   
2. position   
3. ThetaWatterson   
4. ThetaD (nucleotide diversity)   
5. Theta? (singleton category)   
6. ThetaH   
7. ThetaL        

Due to the huge idx file size, I only plot the thetas for chromosome 10 in NEH population.

Overall site depth distribution
<img src="https://hzz0024.github.io/images/QC_CVreseq_maxD2000-page-001.jpg" alt="img" width="800"/>

NEH chr 10 Watterson Theta
<img src="https://hzz0024.github.io/images/Mahattan_NEH_no6inv_fold_chr10_Wt.jpg" alt="img" width="800"/>

NEH chr 10 Tajima Theta
<img src="https://hzz0024.github.io/images/Mahattan_NEH_no6inv_fold_chr10.jpg" alt="img" width="800"/>

I tried to manually calculate the maximum value of Watterson Theta (R code),
```R
watterson_estimator <- function(Sn, n) {
  a_n = 0
  for (i in 1:(n-1)) { 
    a_n = a_n + 1/i
  }
  theta = Sn/a_n
  return(theta)
}

theta = watterson_estimator(Sn=1,n=12)
print(theta)
[1] 0.3311393
```

and compared that to the maximum value in the dataset,
```R
max(exp(DT$Watterson), na.rm = TRUE)
[1] 0.3311393
```

#### Questions

1. what if we have differnt population size, does that impact the theta results? 

2. what is the optimal depth filtering parameters?







