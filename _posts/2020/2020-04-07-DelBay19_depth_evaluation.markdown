---
comments: true
title: DelBay19 depth evaluation
date: '2020-04-07 23:30'
tags:
  - DelBay19
  - depth
  - WGS
categories:
  - WGS data analysis
---

Finally I have the bamfiles ready for angsd analyses. The current dataset is composed by two challenge (ch vs. ref) populations and five wild groups (ARN, COH, HC, NB, SR). The first step is to evaluate the global depth of each datasets. 

---

- scripts for global depth distribution

```sh
module load angsd/0.931

angsd -b all_333.list -anc cv30.fa -out QC_all_maxD2000 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
-minMapQ 20 -minQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 \
-maxDepth 2000 -P 16 >& QC_all.log

angsd -b SR_48.list -anc cv30.fa -out QC_SR_48_maxD250 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
-minMapQ 20 -minQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 \
-maxDepth 250 -P 16 >& QC_SR_48.log

angsd -b NB_48.list -anc cv30.fa -out QC_NB_48_maxD250 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
-minMapQ 20 -minQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 \
-maxDepth 250 -P 16 >& QC_NB_48.log

angsd -b HC_48.list -anc cv30.fa -out QC_HC_48_maxD250 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
-minMapQ 20 -minQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 \
-maxDepth 250 -P 16 >& QC_HC_48.log

angsd -b ARN_47.list -anc cv30.fa -out QC_ARN_47_maxD250 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
-minMapQ 20 -minQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 \
-maxDepth 250 -P 16 >& QC_ARN_47.log

angsd -b COH_44.list -anc cv30.fa -out QC_COH_44_maxD250 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
-minMapQ 20 -minQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 \
-maxDepth 250 -P 16 >& QC_COH_44.log

angsd -b ref_48.list -anc cv30.fa -out QC_ref_48_maxD250 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
-minMapQ 20 -minQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 \
-maxDepth 250 -P 16 >& QC_ref_48.log

angsd -b ch_50.list -anc cv30.fa -out QC_ch_50_maxD250 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
-minMapQ 20 -minQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 \
-maxDepth 250 -P 16 >& QC_ch_50.log

angsd -b ch1_12.list -anc cv30.fa -out QC_ch1_12_maxD250 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
-minMapQ 20 -minQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 \
-maxDepth 250 -P 16 >& QC_ch1_12.log

angsd -b ch2_12.list -anc cv30.fa -out QC_ch2_12_maxD250 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
-minMapQ 20 -minQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 \
-maxDepth 250 -P 16 >& QC_ch2_12.log

angsd -b ch3_12.list -anc cv30.fa -out QC_ch3_12_maxD250 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 \
-minMapQ 20 -minQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 \
-maxDepth 250 -P 16 >& QC_ch3_12.log

angsd -b ch4_14.list -anc cv30.fa -out QC_ch4_14_maxD250 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 
-minMapQ 20 -minQ 20 -doQsDist 1 -doDepth 1 -doCounts 1 \
-maxDepth 250 -P 16 >& QC_ch4_14.log

```

- Estimate the mean and standard deviation from depthGlobal files

I made an R script that could be used for mean and standard deviation calculation using the depthGlobal files

The script can be downloaded [here](https://hzz0024.github.io/scripts/mean_sd.R), or just copy the R code below,

```R
# usage: Rscript mean_sd.R prefix_of_the_depthGlobal 2> /dev/null

fin <- commandArgs(T)
cat("", file=paste(fin,"_mean_sd.log", sep="", collapse=""))

dep <- as.numeric(scan(paste(fin,".depthGlobal", sep = "", collapse=""), what = "char", quiet = T))
cnt = seq(0,length(dep)-1)
Mean = sum((dep * cnt))/(sum(dep))
Deviation = sum((cnt-Mean)^2*dep)/(sum(dep)-1)
Sd = sqrt(Deviation)

cat("\nMean\tDeviation\tSD\n", file=paste(fin,"_mean_sd.log", sep="", collapse=""),append=T)
write.table(cbind(Mean, Deviation, Sd), row.names = F, col.names = F, quote = F, sep = "\t", file=paste (fin,"_mean_sd.log", sep="", collapse=""), append = T)
  
```
Running script

```sh
Rscript mean_sd.R QC_all_maxD2000 2> /dev/null
Rscript mean_sd.R QC_SR_48_maxD250 2> /dev/null
Rscript mean_sd.R QC_NB_48_maxD250 2> /dev/null
Rscript mean_sd.R QC_HC_48_maxD250 2> /dev/null
Rscript mean_sd.R QC_ARN_47_maxD250 2> /dev/null
Rscript mean_sd.R QC_COH_44_maxD250 2> /dev/null
Rscript mean_sd.R QC_ref_48_maxD250 2> /dev/null
Rscript mean_sd.R QC_ch_50_maxD250 2> /dev/null
Rscript mean_sd.R QC_ch1_12_maxD250 2> /dev/null
Rscript mean_sd.R QC_ch2_12_maxD250 2> /dev/null
Rscript mean_sd.R QC_ch3_12_maxD250 2> /dev/null
Rscript mean_sd.R QC_ch4_14_maxD250 2> /dev/null
```
---

### Results

- QC_all_maxD2000
<img src="https://hzz0024.github.io/images/QC_all_maxD2000.jpg" alt="img" width="800"/>

- QC_SR_48_maxD250
<img src="https://hzz0024.github.io/images/QC_SR_48_maxD250.jpg" alt="img" width="800"/>

- QC_NB_48_maxD250
<img src="https://hzz0024.github.io/images/QC_NB_48_maxD250.jpg" alt="img" width="800"/>

- QC_HC_48_maxD250
<img src="https://hzz0024.github.io/images/QC_HC_48_maxD250.jpg" alt="img" width="800"/>

- QC_ARN_47_maxD250
<img src="https://hzz0024.github.io/images/QC_ARN_47_maxD250.jpg" alt="img" width="800"/>

- QC_COH_44_maxD250 
<img src="https://hzz0024.github.io/images/QC_COH_44_maxD250.jpg" alt="img" width="800"/>

- QC_ref_48_maxD250
<img src="https://hzz0024.github.io/images/QC_ref_48_maxD250.jpg" alt="img" width="800"/>

- QC_ch_50_maxD250
<img src="https://hzz0024.github.io/images/QC_ch_50_maxD250.jpg" alt="img" width="800"/>

- QC_ch1_12_maxD250
<img src="https://hzz0024.github.io/images/QC_ch1_12_maxD250.jpg" alt="img" width="800"/>

- QC_ch2_12_maxD250
<img src="https://hzz0024.github.io/images/QC_ch2_12_maxD250.jpg" alt="img" width="800"/>

- QC_ch3_12_maxD250
<img src="https://hzz0024.github.io/images/QC_ch3_12_maxD250.jpg" alt="img" width="800"/>

- QC_ch4_14_maxD250
<img src="https://hzz0024.github.io/images/QC_ch4_14_maxD250.jpg" alt="img" width="800"/>

---

### Conclusion 

The average depth coverage in these data is 1X per individual, and 3 standard deviations from the mean includes all but the thinnest part of the tail. 

The depth results of individual population here (including the ch1-4) is helpful to determine the MaxDepth cutoff for Fst estimation. For population with ~ 50 samples, I would evaluate how the maxDepth impact the Fst patterns by running angsd on two example populations (ch_50 and ref_48). The 2, 3 and 4 standard deviations above the mean will be set as the maxDepth and used for Fst plot generation.
