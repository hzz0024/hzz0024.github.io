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

This post will measure the relatedness for DelBay19 dataset.

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
  

#send the following files to PC for plotting with plot_ngsAdmix.R and angsd_ibs_pca.R:
*.qopt, *Mat, bams

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
zcat ch_ref_98_pca_minI69D49maxD347_minQ30_minMAF05_SNPe6_nochr56invers_70.mafs.gz | cut -f5 |sed 1d > ch_ref_freq
zcat wild_235_pca_minI165D118maxD1003_minQ30_minMAF05_SNPe6_nochr56invers_70.mafs.gz | cut -f5 |sed 1d >wild_freq



