---
comments: true
title: DelBay19 fst permutation and outlier detection
date: '2020-06-22 12:00'
tags:
  - DelBay19
  - Fst
  - ouliter
  - null model
  - WGS
categories:
  - WGS data analysis
---

In this post I will perform the significance test for Fst outliers and examine the deltap among these outliters. The first step is to "label" the samples from two populations and shuffle them to recreate two new population. Next I would generate 1,000 neutral dataset of Fst (i.e. permutation test), and use that dataset for p-value calculation (see below). 

1) For each SNP k: have the observed Fst_k

2) p-value = proportion of resamples that are larger than Fst_k - WINNER! i.e. the proportion of resamples where the Ha is true. Assume that p-value =0.01 -> that means that only 1% of resamples are greater than Fst_k

3) Apply FDR to p-values

Note: the alternative way to calculate the p-value is based on simple t-test. That is, for each neutral Fst_k, we have 1,000 resampled Fst values - this will have a mean and SE. In this case, we could run a one sample t-test (could be a one-side test or two-sided)    
   Ho: mean resampled Fst >= Fst_k    
   Ha: mean resampled Fst < Fst_k 

### Data generation

-minQ 20
-minInd 70% 
-minMapQ 25
-maf 0.20

Step 1 Using all samples of two contrast groups (ch and ref here) for maf file generation. The maf file will be used for minor allele frequency filter. Here I used three filters, 0.05, 0.10, and 0.20 to generate the SNP list for angsd running.

```sh

files=$(ls *.mafs.gz)
for file in $files; do
    file=${file/_doMAF.mafs.gz/}
    zcat $file'_doMAF.mafs.gz' | tail -n +2 | awk '$6>0.05 {print ;}' | awk '{print $1,$2,$3,$4}' > $file'_maf05.snplist'
    zcat $file'_doMAF.mafs.gz' | tail -n +2 | awk '$6>0.10 {print ;}' | awk '{print $1,$2,$3,$4}' > $file'_maf10.snplist'
    zcat $file'_doMAF.mafs.gz' | tail -n +2 | awk '$6>0.20 {print ;}' | awk '{print $1,$2,$3,$4}' > $file'_maf20.snplist'
done

walltime used =  17648.00 sec
```

|    MAF     | No. SNPs |  
| -----------|----------|
|  no filter | 34313641 | 
|    0.05    | 872160   | 
|    0.10    | 539901   |
|    0.20    | 291145   |

Step 2 Generate the saf file for downstream Fst calling. The SNP list generated from step 1 will be used here. Fst dataset from this step is used as observed Fst.

```sh

#ch
module load angsd/0.931
angsd -b /scratch/hzz0024/fst_pt/sample/ch.list -anc /scratch/hzz0024/fst_pt/genome/cv30.fa -out /scratch/hzz0024/fst_pt/output/ch_ref_98_ch_doMAF_filter -dosaf 1 -doMaf 1 -doGlf 2 -GL 1 -doPost 1 -doMajorMinor 1 -doDepth 1 -doCounts 1 -P 20 -minQ 20 -minMapQ 25 -setMinDepth 25 -setMaxDepth 100 -minInd 35 -sites ch_ref_98_maf20.snplist >& /scratch/hzz0024/fst_pt/output/ch_ref_98_ch_doMAF.log
#ref
angsd -b /scratch/hzz0024/fst_pt/sample/ref.list -anc /scratch/hzz0024/fst_pt/genome/cv30.fa -out /scratch/hzz0024/fst_pt/output/ch_ref_98_ref_doMAF_filter -dosaf 1 -doMaf 1 -doGlf 2 -GL 1 -doPost 1 -doMajorMinor 1 -doDepth 1 -doCounts 1 -P 20 -minQ 20 -minMapQ 25 -setMinDepth 25 -setMaxDepth 100 -minInd 35 -sites ch_ref_98_maf20.snplist >& /scratch/hzz0024/fst_pt/output/ch_ref_98_ref_doMAF.log
# realSFS for fst estimates
./get_fst_ch.sh /scratch/hzz0024/fst_pt/output pop_ch.txt 1 _doMAF_filter >& get_fst_ch_06042020.log

```

|            |No. SNPs  | SNPs retained|  
| -----------|----------|--------------|
|    CH      | 520698396| 291145       |
|    REF     | 545117898| 291145       |

Step 3 Perform permutation tests. A python script is developed for this purpose. Each permutation is randomly generated independently of previous permutations or of the original data configuration. It also means that random permutations may include repetitions of the same permutations. This is the most common way that permutation tests are used in practice, especially in the genomic literature. 

Sometimes the p-value is calibrated by adding a 1 to the numerator and denominator to account for misestimation of the p-value (see https://genomicsclass.github.io/book/pages/permutation_tests.html). For more details see [Permutation P-values should never be zero](https://pubmed.ncbi.nlm.nih.gov/21044043/). 

```sh
# write script for angsd run
for i in {0..999}; do
    echo -e 'module load angsd/0.931\nangsd -b /scratch/hzz0024/fst_pt/sample/ch-ref_'$i'a.list -anc /scratch/hzz0024/fst_pt/genome/cv30.fa -out /scratch/hzz0024/fst_pt/output/ch-ref_'$i'a_doMAF_filter -dosaf 1 -doMaf 1 -doGlf 2 -GL 1 -doPost 1 -doMajorMinor 1 -doDepth 1 -doCounts 1 -P 20 -minQ 20 -minMapQ 25 -setMinDepth 25 -setMaxDepth 100 -minInd 35 -sites ch_ref_98_maf20.snplist >& /scratch/hzz0024/fst_pt/output/ch-ref_'$i'a_doMAF.log' >> 'ch-ref_'$i'a.sh'
    echo -e 'module load angsd/0.931\nangsd -b /scratch/hzz0024/fst_pt/sample/ch-ref_'$i'b.list -anc /scratch/hzz0024/fst_pt/genome/cv30.fa -out /scratch/hzz0024/fst_pt/output/ch-ref_'$i'b_doMAF_filter -dosaf 1 -doMaf 1 -doGlf 2 -GL 1 -doPost 1 -doMajorMinor 1 -doDepth 1 -doCounts 1 -P 20 -minQ 20 -minMapQ 25 -setMinDepth 25 -setMaxDepth 100 -minInd 35 -sites ch_ref_98_maf20.snplist >& /scratch/hzz0024/fst_pt/output/ch-ref_'$i'b_doMAF.log' >> 'ch-ref_'$i'b.sh'
done
# write script for pop txt and get fst
for i in {0..999}; do
    echo -e 'Pop\nch-ref_'$i'a\nch-ref_'$i'b' >> 'pop_ch_'$i'.txt'
    echo -e './get_fst_ch.sh /scratch/hzz0024/fst_pt/output_random pop_ch_'$i'.txt 1 _doMAF_filter >& get_fst_ch_'$i'_06042020.log' >> 'get_fst_ch_'$i'_run.sh'
done
# write script for job submission
files=$(ls *.sh)
for file in $files; do
    file=${file/.sh/}
    echo -e 'EMAIL=`whoami`"@XXXXX.edu";\nCWD=`pwd`;\nqsub -q general -N '$file'_06032020 -j oe -e '$file'.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V '$file'.sh\nexit 0;' >> 'run_'$file'_06042020.sh'
done
```

 - for each run, the saf step will cost > 3 hours, while the realSFS step costs ~ 15 mins.

Step 4 Compare the Fst values between observed and neutral datasets. Perform statistic analysis for each SNP using the Fst_cnt.py.

- Average weighted Fst comparsion 

|            | Observed | Neutral (averaged from 1,000 runs)     |  
| -----------|----------|--------------|
|Average Fst | 0.001167 | 0.000778     |

Note: the mean weighted fst ranges from 0.000328 to 0.001694

Step 5 For each shared SNP (291,145 in total), I selected the snp outliers by calculating the proportion of resamples that are larger than Fst_obs. That proportion would be the p-value.

|     P-value    |  No. of outliers  |  
| ---------------|-------------------|
|      0.0100    |      3254         |
|      0.0075    |      2626         |
|      0.0050    |      1673         |
|      0.0025    |      1054         |
|      0.0010    |      386          |

Step 6 Examine the deltap change. 

1) extract the snp outliers identifed from permutation test. This is done by creating a file named "outlier_list.txt", followed by *extract.py* script. The mafs files for both populations (e.g. challenge and reference) are also needed for script running.

```sh
python extract.py

# add header to the outputs
chromo  position  major minor anc knownEM nInd
```

2) calculate the deltap using R script

Two versions of deltap were available, one for absolute value, another for actual deltap. 

example usage: 

```sh
Rscript get_deltaP.R -d /workdir/xxx/ -p 1.mafs.gz -q 2.mafs.gz -t 291145 -o deltap.output

#for obs data
Rscript deltaP_abs.R -d /Volumes/cornell/DelBay19_Hopper/permutation/3_deltap -p ch_ref_98_ref_doMAF_filter.mafs.extracted -q ch_ref_98_ch_doMAF_filter.mafs.extracted -t 291145 -o obs_deltap.output

#for neutral data 

```

- outliers with p-value < 0.001 (386 snps)

1) create the outlier list

2) extract the outlier maf information

```py
python extract_neu.py
```

3) add the header to the extracted file

```sh
./addname.sh
```

4) calculate deltap (absolute value) for each neutral data using deltaP_abs.R

```sh
./repeat_run.sh
```

5) compare the observed deltap and neutral deltap (from 1000 permutations) for 386 SNPs

```py
python deltap_cnt.py

# need to change the No. of SNPs at line 21
```

Note that Fst outliers with p-value < 0.001 are 386, while deltap (based on absolute values) outliers with p-value < 0.001 are 259 out of 386 SNPs. I will use another post to show the deltap outlier identification from permutation test, and compare that results to single-generation selection tests (i.e. SGS used by sea bream paper)

6) deltap plotting

- absolute deltap plotting for all 386 SNPs, here the red dots means observation deltap (absolute value for each SNP), the black dots means medium deltap (absolute, 50% quantile) value calculated from the 1000 permutation test. The yellow part covers the 25% to 99% quantile range for the neutral deltap (absolute).

<img src="https://hzz0024.github.io/images/outlier/386_deltap.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/outlier/obs_deltap_p.jpeg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/outlier/neu_deltap_p.jpeg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/outlier/386_mac.jpeg" alt="img" width="800"/>




   

