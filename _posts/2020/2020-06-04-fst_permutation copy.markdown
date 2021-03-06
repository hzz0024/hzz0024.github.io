---
comments: true
title: DelBay19 fst outlier detection using permutation method
date: '2020-06-04 12:00'
tags:
  - DelBay19
  - Fst
  - ouliter
  - WGS
categories:
  - WGS data analysis
---

In this post I tried to perform the significance test for Fst outliers. This could be done by randomly drawing individauls from two populations and obtain a neutral dataset of Fst (or permutation test), and using that dataset for p-value calculation. The workframe is shown below,

1) For each SNP k: have the observed Fst_k

2) Repeated scrambled the subpopulation membership to get Fst_k_hat, say M = 100 times. For each Fst_k, we have 100 resampled Fst_k_hat values - this will have a mean and SE

3) Or p-value = proportion of resamples that are larger than Fst_k - WINNER! i.e. the proportion of resamples where the Ha is true. Assume that p-value =0.01 -> that means that only 1% of resamples are greater than Fst_k

4) we could run a one sample t-test (could be a one-side test or two-sided)    
   Ho: mean resampled Fst >= Fst_k    
   Ha: mean resampled Fst < Fst_k   

5) Apply FDR to p-values

### Data generation

-minQ 20
-minInd 70% 
-minMapQ 25

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

| 	 MAF	   | No. SNPs |  
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
|    CH      | 520698396| 243239       |
|    REF     | 545117898| 188681       |

Here the SNP number should be same, in my next round test I will remove the minInd settings to keep the SNP number same as the snplist

Step 3 Perform the LD pruning using ngsLD tool

Step 4 Perform permutation test by randomly drawing individauls from two populations and obtain a neutral dataset of Fst. A python script is developed for individual shuffle. The script is designed to ensure that each redraw has a equal number of samples from each population. 

```sh
# write script for angsd run
for i in {0..99}; do
    echo -e 'module load angsd/0.931\nangsd -b /scratch/hzz0024/fst_pt/sample/ch-ref_'$i'a.list -anc /scratch/hzz0024/fst_pt/genome/cv30.fa -out /scratch/hzz0024/fst_pt/output/ch-ref_'$i'a_doMAF_filter -dosaf 1 -doMaf 1 -doGlf 2 -GL 1 -doPost 1 -doMajorMinor 1 -doDepth 1 -doCounts 1 -P 20 -minQ 20 -minMapQ 25 -setMinDepth 25 -setMaxDepth 100 -minInd 35 -sites ch_ref_98_maf20.snplist >& /scratch/hzz0024/fst_pt/output/ch-ref_'$i'a_doMAF.log' >> 'ch-ref_'$i'a.sh'
    echo -e 'module load angsd/0.931\nangsd -b /scratch/hzz0024/fst_pt/sample/ch-ref_'$i'b.list -anc /scratch/hzz0024/fst_pt/genome/cv30.fa -out /scratch/hzz0024/fst_pt/output/ch-ref_'$i'b_doMAF_filter -dosaf 1 -doMaf 1 -doGlf 2 -GL 1 -doPost 1 -doMajorMinor 1 -doDepth 1 -doCounts 1 -P 20 -minQ 20 -minMapQ 25 -setMinDepth 25 -setMaxDepth 100 -minInd 35 -sites ch_ref_98_maf20.snplist >& /scratch/hzz0024/fst_pt/output/ch-ref_'$i'b_doMAF.log' >> 'ch-ref_'$i'b.sh'
done
# write script for pop txt and get fst
for i in {0..99}; do
    echo -e 'Pop\nch-ref_'$i'a\nch-ref_'$i'b' >> 'pop_ch_'$i'.txt'
    echo -e './get_fst_ch.sh /scratch/hzz0024/fst_pt/output_random pop_ch_'$i'.txt 1 _doMAF_filter >& get_fst_ch_'$i'_06042020.log' >> 'get_fst_ch_'$i'_run.sh'
done
# write script for job submission
files=$(ls *.sh)
for file in $files; do
    file=${file/.sh/}
    echo -e 'EMAIL=`whoami`"@auburn.edu";\nCWD=`pwd`;\nqsub -q general -N '$file'_06032020 -j oe -e '$file'.error -l nodes=1:ppn=20,mem=120GB,walltime=240:00:00 -m be -M $EMAIL -d $CWD -V '$file'.sh\nexit 0;' >> 'run_'$file'_06042020.sh'
done
```

 - for each run, the saf step will cost > 3 hours, while the realSFS step costs ~ 15 mins.


Step 4 Compare the Fst outputs between observed and neutral datasets. Perform statistic analyses on shared SNPs.

- Average weighted Fst comparsion 

|            | Observed | Neutral      |  
| -----------|----------|--------------|
|Average Fst | 0.001167 | 0.000446     |

- Fst distribution 

<img src="https://hzz0024.github.io/images/slim/fst_density_all.jpeg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/slim/fst_density.jpeg" alt="img" width="800"/>

- qqplot

<img src="https://hzz0024.github.io/images/slim/qqplot.jpeg" alt="img" width="800"/>


Some literatures showed alternative ways to achieve this goal,

[Experimental evidence for ecological selection on genome variation in the wild](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12238)

<img src="https://hzz0024.github.io/images/outlier/index.jpg" alt="img" width="800"/>

> Figure S2. Schematic of the null models of the absence of selection (null models 1 and 2). Each box or circle contains the genetic data for a hypothetical locus. Black numbers or letters denote individuals that survived and gray characters denote individuals that lived. A) In null model 1 the squares denote the true survivors and circles denote simulated sets of survivors. Density plots show the simulated distribution of allele frequency change for each locus and the vertical blue lines denote the observed allele frequency change. We reject the null hypothesis in model 1 for loci where the observed change (blue line) is an extreme tail of this null distribution. B) In null model 2 the real and replicate simulated data sets are show in boxes and numbers rather than genotypes are used to highlight the fact that the survival data is applied to all loci. Each of these data sets (the true data set and 100 replicate simulated data sets) are subjected to the null model 1 analysis (i.e., pane A) and the number of exceptional allele frequency change loci designated from the true survival data (vertical blue line) and replicate simulated data sets (gray histogram) are tabulated to determine whether the null hypothesis of no selection across the genome (null model 2) can be rejected.

[Within-Generation Polygenic Selection Shapes Fitness-Related Traits across Environments in Juvenile Sea Bream](https://www.mdpi.com/2073-4425/11/4/398/htm#app1-genes-11-00398)





   

