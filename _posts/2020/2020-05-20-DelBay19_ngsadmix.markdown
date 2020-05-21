---
comments: true
title: ngsadmix analyese
date: '2020-05-13 12:00'
tags:
  - DelBay19
  - angsd
  - ADMIXTURE
  - relateness
categories:
  - DelBay19 data analysis
---

This post will summarize the admixture results in DelBay19 dataset.

Some scripts are adopted from [here](https://github.com/alexkrohn/AmargosaVoleTutorials/blob/master/ngsAdmix_tutorial.md)

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
module load ngsadmix/1
# loop over ch-ref populations
for i in {1..8}
do
    for j in {1..10}
    do
        NGSadmix -likes ch_ref_98_pca_minI69D49maxD347_minQ30_minMAF05_SNPe6_nochr56invers_70.beagle.gz -K $i -o chr_k${i}_run${j}
    done
done
#loop over wild populations
module load ngsadmix/1

for i in {1..8}
do
    for j in {1..10}
    do
        NGSadmix -likes wild_235_pca_minI165D118maxD1003_minQ30_minMAF05_SNPe6_nochr56invers_70.beagle.gz -K $i -o wild_k${i}_run${j}
    done
done

#outputs (one for each k):
  mydata_k2.filter
  mydata_k2.fopt.gz   -- estimation of the allele frequencies in each of the ancestral populations (each column is a population, each row a SNP)
  mydata_k2.log       -- log of the run
  mydata_k2.qopt      -- Infered admixture proportions. Each line is an individual and each column is a population.

```
> Create the CLUMPAK Input File to Calculate the Best K and make barplot
```sh
#only k1-4 now, I will include other run results later
logs <- as.data.frame(read.table("logfile"))
logs$K <- c(rep("1", 7), rep("2", 7), rep("3", 7), rep("4", 7))
write.table(logs[, c(2, 1)], "logfile_formatted", row.names = F, 
            col.names = F, quote = F)

# start plot
#send the following files to PC for plotting with plot_ngsAdmix.R 
# start plot
admix <- t(as.matrix(read.table("chr_k2_run3.qopt")))
head(admix[, 1:98])
ref <- read.table('challenge_98.list', head=TRUE)
populations <- as.vector(ref[,3])
# sort by population
orders <- order(populations, decreasing=FALSE)
admix <- admix[,orders]
populations <- populations[orders]
# draw abline
pops = c("CH1", "CH2", "CH3", "CH4", "REF")
idx = c(rep(0,5))
for (i in c(1:length(populations))){
  for (j in c(1:length(pops))) {
    if (populations[i] == pops[j]){
      idx[j] = i
    }
    
  }
}

orders1 <- order(idx, decreasing=FALSE)
pops <- pops[orders1]
idx <- idx[orders1]
text_loc = c()
runner = 0
for (i in idx){
    text_loc <- c(text_loc, (runner+i)/2)
    runner = i
}

barplot(admix, col = c("blue", "green",'black','red'), space = 0, border = NA, 
        ylab = "Admixture", main = "ch-ref group (K=2)")
text(text_loc, -0.05, unique(populations), 
     xpd = T)
abline(v = idx, lty = 5, lwd = 2, col = "white")

```

### plot the admixture results

ch_ref_k2_run1
<img src="https://hzz0024.github.io/images/ngsadmix/chr_k2_run1.jepg" alt="img" width="800"/>

ch_ref_k2_run2
<img src="https://hzz0024.github.io/images/ngsadmix/chr_k2_run2.jepg" alt="img" width="800"/>

ch_ref_k2_run3
<img src="https://hzz0024.github.io/images/ngsadmix/chr_k2_run3.jepg" alt="img" width="800"/>

ch_ref_k2_run4
<img src="https://hzz0024.github.io/images/ngsadmix/chr_k2_run4.jepg" alt="img" width="800"/>

By looking at the admixture plots with equal K, the q-value for the same individual varies a lot. Need double check the ngsadmix again to see why this happens.

DeltaK can only calculate best K with K > 1, but based on likelihood estimate itself, the best K is 4. Not very informative based on admixature plot with K=4. 

best k from likelihood estimates
<img src="https://hzz0024.github.io/images/ngsadmix/bestk.jpg" alt="img" width="800"/>

ch_ref_k1_run1
<img src="https://hzz0024.github.io/images/ngsadmix/chr_k4_run1.ejpg" alt="img" width="800"/>



