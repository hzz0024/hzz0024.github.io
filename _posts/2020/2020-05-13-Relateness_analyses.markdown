---
comments: true
title: Relateness analyese
date: '2020-05-13 12:00'
tags:
  - CVreseq
  - DelBay19
  - angsd
  - ADMIXTURE
  - relateness
categories:
  - WGS data analysis
---

This post will summarize the measurements the relatedness in DelBay19 dataset.

Some scripts are adopted from [here](https://github.com/grovesdixon/caveRAD/blob/master/cave_RAD_processing_walkthrough.txt)

--- 

### Details of data analyses

```sh

module load angsd/0.931

# angsd run
# Here I target for 
# -minIND=69 #70% of 98 challenge samples and MININD=165 #70% of 235 wild samples
# -minMapQ 25
angsd -b challenge_98.list -anc cv30.fa -ref cv30.fa -out mds_output/ch_ref_98_pca_minI69D49maxD347_minQ30_minMAF05_SNPe6_nochr56invers_70 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doCov 1 -makeMatrix 1 -doIBS 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 20 -minQ 20 -minMapQ 25 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 49 -setMaxDepth 347 -minInd 69 -doMajorMinor 3 -sites ch_ref_98_pca_minI49D49maxD347_minQ20_minMAF05_SNPe6_nochr56invers.snplist -rf chr.list >& mds_output/ch_ref_98_pca_70_nochr56invers.log

angsd -b wild_235.list -anc cv30.fa -ref cv30.fa -out mds_output/wild_235_pca_minI165D118maxD1003_minQ30_minMAF05_SNPe6_nochr56invers_70 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doPost 1 -doCov 1 -makeMatrix 1 -doIBS 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 20 -minQ 20 -minMapQ 25 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 118 -setMaxDepth 1003 -minInd 165 -doMajorMinor 3 -sites wild_235_pca_minI118D118maxD1003_minQ20_minMAF05_SNPe6_nochr56invers.snplist -rf chr.list >& mds_output/wild_235_pca_70_nochr56invers.lo

# Number of SNPs for each output

zcat ch_ref_98_pca_minI69D49maxD347_minQ30_minMAF05_SNPe6_nochr56invers_70.mafs.gz | wc -l
# 1552833
zcat wild_235_pca_minI165D118maxD1003_minQ30_minMAF05_SNPe6_nochr56invers_70.mafs.gz | wc -l
# 1982310

#--------ADMIXTURE---------#
#USE NGSadmix ON GENOTYPE LIKELIHOODS TO INCLUDE UNCERTAINTY

for K in `seq 2 5` ; 
do 
NGSadmix -likes ibsResults.beagle.gz -K $K -P 10 -o mydata_k${K};
done


#outputs (one for each k):
  mydata_k2.filter
  mydata_k2.fopt.gz   -- estimation of the allele frequencies in each of the ancestral populations (each column is a population, each row a SNP)
  mydata_k2.log       -- log of the run
  mydata_k2.qopt      -- Infered admixture proportions. Each line is an individual and each column is a population.
  

#send the following files to PC for plotting with plot_ngsAdmix.R 
*.qopt
```
still under running

```sh

#--------RELATNEDNESS----------#

#First have to re-run angsd with '-doGlf 3' argument:
# -doGlf  0
# 1: binary glf (10 log likes)  .glf.gz
# 2: beagle likelihood file .beagle.gz
# 3: binary 3 times likelihood  .glf.gz
# 4: text version (10 log likes)  .glf.gz

module load angsd/0.931

# the only change from 1_mds.sh is the change of -doGlf 2 to -doGlf 3

angsd -b challenge_98.list -anc cv30.fa -ref cv30.fa -out relat_output/ch_ref_98_pca_minI69D49maxD347_minQ30_minMAF05_SNPe6_nochr56invers_70 -dosaf 1 -GL 1 -doGlf 3 -doMaf 1 -doPost 1 -doGeno 8 -doCov 1 -makeMatrix 1 -doIBS 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 20 -minQ 20 -minMapQ 25 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 49 -setMaxDepth 347 -minInd 69 -doMajorMinor 3 -sites ch_ref_98_pca_minI49D49maxD347_minQ20_minMAF05_SNPe6_nochr56invers.snplist -rf chr.list >& relat_output/ch_ref_98_pca_70_nochr56invers.log

angsd -b wild_235.list -anc cv30.fa -ref cv30.fa -out relat_output/wild_235_pca_minI165D118maxD1003_minQ30_minMAF05_SNPe6_nochr56invers_70 -dosaf 1 -GL 1 -doGlf 3 -doMaf 1 -doPost 1 -doGeno 8 -doCov 1 -makeMatrix 1 -doIBS 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 20 -minQ 20 -minMapQ 25 -minMaf 0.05 -remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -SNP_pval 1e-6 -setMinDepth 118 -setMaxDepth 1003 -minInd 165 -doMajorMinor 3 -sites wild_235_pca_minI118D118maxD1003_minQ20_minMAF05_SNPe6_nochr56invers.snplist -rf chr.list >& relat_output/wild_235_pca_70_nochr56invers.log

#now get relatedness

# extract the snp frequency data from mafs zip file
# note here use cut -f7, which is different from the cave_RAD post
zcat ch_ref_98_pca_minI69D49maxD347_minQ30_minMAF05_SNPe6_nochr56invers_70.mafs.gz | cut -f7 |sed 1d > ch_ref_freq
zcat wild_235_pca_minI165D118maxD1003_minQ30_minMAF05_SNPe6_nochr56invers_70.mafs.gz | cut -f7 |sed 1d >wild_freq

# run through ngsrelate
module load ngsrelate

ngsRelate -f ch_ref_freq -g ch_ref_98_pca_minI69D49maxD347_minQ30_minMAF05_SNPe6_nochr56invers_70.glf.gz -n 98 -z challenge_98.list -O ch_output.res

ngsRelate -f wild_freq -g wild_235_pca_minI165D118maxD1003_minQ30_minMAF05_SNPe6_nochr56invers_70.glf.gz -n 235 -z wild_235.list -O wild_output.res

# example output for challenge group

head ch_output.res 

a b ida idb nSites  J9  J8  J7  J6  J5  J4  J3  J2  J1  rab Fa  Fb  theta inbred_relatedness_1_2  inbred_relatedness_2_1  fraternity  identity  zygosity  2of3_IDB  F_diff_a_b  loglh nIter bestoptimll coverage  2dsfs R0  R1  KING  2dsfs_loglike 2dsfsf_niter
0 1 SR0719ch_002.bam  SR0719ch_019.bam  666348  0.927751  0.000000  0.000000  0.000552  0.000000  0.049241  0.000000  0.022456  0.000000  0.000000  0.071697  0.023008  0.000000  0.000000  0.000000  0.022456  0.000000  0.022456  0.047353  0.024344  -709136.469929  70  -1  0.429118  5.280951e-01,1.551008e-01,2.511597e-02,1.435009e-01,7.577034e-02,1.671133e-02,2.723838e-02,1.973095e-02,8.736163e-03  0.690961  0.195588  -0.059472 -790567.826697  43
0 2 SR0719ch_002.bam  SR0719ch_064.bam  596550  0.870193  0.000000  0.000000  0.058360  0.000000  0.034342  0.000000  0.037105  0.000000  0.000000  0.071447  0.095465  0.000000  0.000000  0.000000  0.037105  0.000000  0.037105  0.083456  -0.012009 -616565.041839  170 -1  0.384169  5.349475e-01,1.433755e-01,3.026484e-02,1.482602e-01,6.737879e-02,2.015998e-02,2.788225e-02,1.899071e-02,8.740264e-03  0.862988  0.173240  -0.105072-686932.450194 39

# that is a lot of stats summary
```

By checking the ngsrelate website (https://github.com/ANGSD/NgsRelate), some description about output format is shown below

> The first two columns contain indices of the two individuals used for the analysis. The third column is the number of genomic sites considered. The following nine columns are the maximum likelihood (ML) estimates of the nine jacquard coefficients, where K0==J9; K1==J8; K2==J7 in absence of inbreeding. Based on these Jacquard coefficients, NgsRelate calculates 11 summary statistics:

> 13 rab is the pairwise relatedness (J1+J7+0.75*(J3+J5)+.5*J8) Hedrick et al   
> 14 Fa is the inbreeding coefficient of individual a J1+J2+J3+J4 Jacquard   
> 15 Fb is the inbreeding coefficient of individual b J1+J2+J5+J6 Jacquard   
> 16 theta is the coefficient of kinship J1 + 0.5*(J3+J5+J7) + 0.25*J8) Jacquard   
> 17 inbred_relatedness_1_2 J1+0.5*J3 Ackerman et al   
> 18 inbred_relatedness_2_1 J1+0.5*J5 Ackerman et al   
> 19 fraternity J2+J7 Ackerman et al   
> 20 identity J1 Ackerman et al   
> 21 zygosity J1+J2+J7 Ackerman et al   
> 22 Two-out-of-three IBD J1+J2+J3+J5+J7+0.5*(J4+J6+J8) Miklos csuros   
> 23 Inbreeding difference 0.5*(J4-J6) Miklos csuros   
> 24 the log-likelihood of the ML estimate.   
> 25 number of EM iterations. If a -1 is displayed. A boundary estimate had a higher likelihood.   
> 26 fraction of sites used for the ML estimate   

The remaining columns relate to statistics based on a 2D-SFS.

> 27 2dsfs estimates using the same methodology as implemented in realSFS, see ANGSD   
> 28 R0 Waples et al   
> 29 R1 Waples et al   
> 30 KING Waples et al   
> 31 the log-likelihood of the 2dsfs estimate.   
> 32 number of iterations for 2dsfs estimate   

