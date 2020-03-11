---
comments: true
title: DelBay1_angsd_SFS_FST_2
date: '2020-03-10 12:00'
tags:
  - angsd
  - SFS
  - fst
  - paralogs
  - WGS
categories:
  - WGS data analysis
---

THis is an update about the DelBay1_angsd_SFS_FST post,

Here are posts I found that are useful for this goal:

1. [Calculating FST: when to fold SFS](https://github.com/ANGSD/angsd/issues/259)

2. [NaN in fst values](https://github.com/ANGSD/angsd/issues/274)

3. [thetaStat fails with folded saf](https://github.com/ANGSD/angsd/issues/286)

4. [fuzzy calculator of probability that heterozygotes constitute](https://github.com/z0on/2bRAD_denovo/blob/master/HetMajorityProb.py)

5. [2d SFS Estimation](http://popgen.dk/angsd/index.php/2d_SFS_Estimation)

and of course, Matt's notes!

First of all, I followed the steps of notes and got the individual pop beagle files (adjust minInd and minDepth and maxDepth) - needed for SFS and pairwise Fst calcs

### Get beagle file
```
# Notes from Matt

/programs/angsd20191002/angsd/angsd -b FIS.bamlist -anc genome/6565_refCV_10fixchr.fa.gz -out angsd_mq20/FIS_minI30D30maxD100_MQ20_minMAF05_SNPe6_no56inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 60 -minInd 13 -sites angsd_mq20/FISHHSV_minMAF05_MQ20_4col_nochr56invers.snplist >& angsd_mq20/FIS_minMAF05.log &
#	[1] 12671, 2 hrs
#	    -> Total number of sites analyzed: 397860158
#       -> Number of sites retained after filtering: 1225591 
/programs/angsd20191002/angsd/angsd -b SV.bamlist -anc genome/6565_refCV_10fixchr.fa.gz -out angsd_mq20/SV_minI30D30maxD100_MQ20_minMAF05_SNPe6_no56inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 60 -minInd 13 -sites angsd_mq20/FISHHSV_minMAF05_MQ20_4col_nochr56invers.snplist >& angsd_mq20/SV_minMAF05.log &
#	[2] 12908, 2 hrs
#        -> Total number of sites analyzed: 457600436
#        -> Number of sites retained after filtering: 1863128
/programs/angsd20191002/angsd/angsd -b HH.bamlist -anc genome/6565_refCV_10fixchr.fa.gz -out angsd_mq20/HH_minI30D30maxD100_MQ20_minMAF05_SNPe6_no56inv -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 3 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 2 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 60 -minInd 13 -sites angsd_mq20/FISHHSV_minMAF05_MQ20_4col_nochr56invers.snplist >& angsd_mq20/HH_minMAF05.log &
#	[3] 13194, 2 hrs
#       -> Total number of sites analyzed: 453391700
#        -> Number of sites retained after filtering: 1785526
# HALF THE SITES THAT I WAS GETTING BEFORE MQ20 FILTERING

# below are my script contents:

# for challenge populations

# load the most recent angsd github release, angsd version: 0.931-18-g37be4d5 (htslib: 1.5-20-ga159aa4) build(Mar  4 2020 08:26:29) 

module load angsd/0.931

angsd -b ch_48.bamlist -anc cv30.fa -out Fst_ch/ch_no48inv_minI2D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12
# 6 hrs
# Total number of sites analyzed: 530191757
# Number of sites retained after filtering: 113396147

angsd -b r_47.bamlist -anc cv30.fa -out Fst_ch/ref_no47inv_minI2D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12
# 6 hrs
# Total number of sites analyzed: 519303947
# Number of sites retained after filtering: 124171639

# for wild populations

module load angsd/0.931

angsd -b ARN_44.bamlist -anc cv30.fa -out Fst_wild_mph/ARN_no44inv_minI12D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12 
# 6 hrs
# Total number of sites analyzed: 536691844
# Number of sites retained after filtering: 94713203

angsd -b COH_42.bamlist -anc cv30.fa -out Fst_wild_mph/COH_no42inv_minI12D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12 
# 6 hrs
# Total number of sites analyzed: 539000040
# Number of sites retained after filtering: 94915508

angsd -b HC_44.bamlist -anc cv30.fa -out Fst_wild_mph/HC_no44inv_minI12D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12 
# 6 hrs
# Total number of sites analyzed: 539427341
# Number of sites retained after filtering: 96294679

angsd -b NB_47.bamlist -anc cv30.fa -out Fst_wild_mph/NB_no47inv_minI12D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12 
# 6 hrs
# Total number of sites analyzed: 534922807
# Number of sites retained after filtering: 100848073

angsd -b SR_48.bamlist -anc cv30.fa -out Fst_wild_mph/SR_no48inv_minI12D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12 
# 6 hrs
# Total number of sites analyzed: 545667117
# Number of sites retained after filtering: 94714818

```	
### Get FST
```	
# Notes from Matt

# modify get_fst script for window sizes you want

./get_fst_windows.sh /workdir/mph75/angsd_mq20 /workdir/mph75/3pops_names.txt 1 _minI30D30maxD100_MQ20_minMAF05_SNPe6_no56inv >& Fst_MQ20_log &

# this step is memory-cost process, may need require more memory than usual

# using most recent angsd: /programs/angsd20191002/angsd/misc/realSFS

# [1] 10515 started 16:55 am, finished 17:45 (50 min), includes individual SNP and windowed output
#		  Fst		weighted Fst (from log files)
#FIS x HH: 0.103902	0.113671
#FIS x SV: 0.094523	0.102044
#HH x SV:  0.006597	0.007601
#Need to convert windowed output to Manhattan plot format
#region	chr	midPos	Nsites
#(19,62)(17786,29325)(15000,30000)	NC_035780.1	22500	45	0.002800
#(63,151)(30473,44941)(30000,45000)	NC_035780.1	37500	90	0.010688
#(152,216)(45027,53309)(45000,60000)	NC_035780.1	52500	66	0.001701
#to
#chr pos_start pos_end Angsd_HHxSV_Fst
#NC_035780.1 1535 1536 0.00843683	#except now I want the full 15K interval reflected, not single SNP
#NC_035780.1 1862 1863 0.0216682
#NC_035780.1 2129 2130 0.0509635
```
I used get_fst script to estimate the SFS for challenge and wild groups

Here I compared the fold vc unfold option and examined how it impact the SFS and Fst outputs, below are results,

### wild groups

|  | Fst | weighted Fst |
| -----| ----| --------|
| ARN x COH (fold) | 0.003315 | 0.003711 |
| ARN x COH (unfold) | 0.003314 | 0.003736 |
| ch x ref (fold) | 0.005207 | 0.002855 |
| ch x ref (unfold) | 0.005270 | 0.002872 |

```shell
awk -F"\t" '{print $2, $3, $3, $5}' ARN_COH_MQ20_minMAF05_SNPe6_fold.1kb_win_1kb_fst.txt | awk '$2-=500' | awk '$3+=500' > ARN_COH_1kb_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' ARN_COH_MQ20_minMAF05_SNPe6_fold.5kb_win_5kb_fst.txt | awk '$2-=2500' | awk '$3+=2500' > ARN_COH_5kb_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' ARN_COH_MQ20_minMAF05_SNPe6_fold.15kb_win_15kb_fst.txt | awk '$2-=7500' | awk '$3+=7500' > ARN_COH_15kb_fold.fst

awk -F"\t" '{print $2, $3, $3, $5}' ARN_COH_MQ20_minMAF05_SNPe6.1kb_win_1kb_fst.txt | awk '$2-=500' | awk '$3+=500' > ARN_COH_1kb_unfold.fst
awk -F"\t" '{print $2, $3, $3, $5}' ARN_COH_MQ20_minMAF05_SNPe6.5kb_win_5kb_fst.txt | awk '$2-=2500' | awk '$3+=2500' > ARN_COH_5kb_unfold.fst
awk -F"\t" '{print $2, $3, $3, $5}' ARN_COH_MQ20_minMAF05_SNPe6.15kb_win_15kb_fst.txt | awk '$2-=7500' | awk '$3+=7500' > ARN_COH_15kb_unfold.fst

awk -F"\t" '{print $2, $3, $3, $5}' ch_ref_MQ20_minMAF05_SNPe6.1kb_win_1kb_fst.txt | awk '$2-=500' | awk '$3+=500' > ch_ref_1kb_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' ch_ref_MQ20_minMAF05_SNPe6.5kb_win_5kb_fst.txt | awk '$2-=2500' | awk '$3+=2500' > ch_ref_5kb_fold.fst
awk -F"\t" '{print $2, $3, $3, $5}' ch_ref_MQ20_minMAF05_SNPe6.15kb_win_15kb_fst.txt | awk '$2-=7500' | awk '$3+=7500' > ch_ref_15kb_fold.fst

```
then edit headers manually: chr window_start window_end angsd_Fst
then edit headers manually (single SNP): chr	pos	c1	c2	angsd_Fst

### qqman Manhattan plotting

NOTE - for qqman Manhattan plotting the chr need to be integer

```shell 
# using sed to replace the chr name with integer
for i in *.fst; do
sed -i .bak 's/NC_035780.1/1/g;s/NC_035781.1/2/g;s/NC_035782.1/3/g;s/NC_035783.1/4/g;s/NC_035784.1/5/g;s/NC_035785.1/6/g;s/NC_035786.1/7/g;s/NC_035787.1/8/g;s/NC_035788.1/9/g;s/NC_035789.1/10/g' $i;
done
```
Plot with Manhattan_loop.R (several loop options in script)

```R
library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)

require(data.table)

##### script for single-SNP plot (due to difficulty in opening the pdf, I export the jepg plot here)
setwd("~/Documents/Ryan_workplace/oyster/WGS/WGS_fst/Fst_challenge/fold/plot/single_snp")
DT <- fread("ch_ref_MQ20_minMAF05_SNPe6.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#pdf("Mahattan_ch_ref_single_SNP_fold.pdf",width=15,height=10)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_single_SNP_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="pos",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref angsd Fst", cex.lab=1.4) #main = "Chromosome",
dev.off()

setwd("~/Documents/Ryan_workplace/oyster/WGS/WGS_fst/Fst_challenge/unfold/plot/single_snp")
DT <- fread("ch_ref_MQ20_minMAF05_SNPe6.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#pdf("Mahattan_ch_ref_single_SNP_unfold.pdf",width=15,height=10)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_single_SNP_unfold.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="pos",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref angsd Fst", cex.lab=1.4) #main = "Chromosome",
dev.off()

##### script for chromosome-wide plots at 1kb, 5kb, and 15kb (in one window)
setwd("~/Documents/Ryan_workplace/oyster/WGS/WGS_fst/Fst_challenge/fold/plot")

DT1 <- fread("ch_ref_1kb_fold.fst")
DT2 <- fread("ch_ref_5kb_fold.fst")
DT3 <- fread("ch_ref_15kb_fold.fst")
DT1$chr <- as.numeric(DT1$chr)
DT2$chr <- as.numeric(DT2$chr)
DT3$chr <- as.numeric(DT3$chr)
#pdf("Mahattan_ch_ref_fold.pdf",width=15,height=10)
jpeg("Mahattan_ch_ref_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
#png("Mahattan_ch_ref_fold.png", width = 6, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1)) 
manhattan(DT1,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref 1kb Fst", cex.lab=1.4) #main = "Chromosome",

manhattan(DT2,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref 5kb Fst", cex.lab=1.4) #main = "Chromosome",

manhattan(DT3,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref 15kb Fst", cex.lab=1.4) #main = "Chromosome",
dev.off()

setwd("~/Documents/Ryan_workplace/oyster/WGS/WGS_fst/Fst_challenge/unfold/plot")
DT1 <- fread("ch_ref_1kb_unfold.fst")
DT2 <- fread("ch_ref_5kb_unfold.fst")
DT3 <- fread("ch_ref_15kb_unfold.fst")
DT1$chr <- as.numeric(DT1$chr)
DT2$chr <- as.numeric(DT2$chr)
DT3$chr <- as.numeric(DT3$chr)
#pdf("Mahattan_ch_ref_unfold.pdf",width=15,height=10)
jpeg("Mahattan_ch_ref_unfold.jpg", width = 16, height = 9, units = 'in', res = 300)
#png("Mahattan_ch_ref_unfold.png", width = 6, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1)) 
manhattan(DT1,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref 1kb Fst", cex.lab=1.4) #main = "Chromosome",

manhattan(DT2,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref 5kb Fst", cex.lab=1.4) #main = "Chromosome",

manhattan(DT3,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref 15kb Fst", cex.lab=1.4) #main = "Chromosome",
dev.off()

```

|  | Fst | weighted Fst |
| -----| ----| --------|
| ARN x COH (fold)| 0.003315 | 0.003711 |



|  | Fst | weighted Fst |
| -----| ----| --------|
| ARN x COH (unfold)| 0.003314 | 0.003736 |



|  | Fst | weighted Fst |
| -----| ----| --------|
| ch x ref (fold)| 0.005207 | 0.002855 |

1k (top), 5k (middle), 15k (bottom) windowed Fst plot 

<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_fold.jpg" alt="img" width="800"/>

single SNP based on Fst plot

<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_single_SNP_fold.jpg" alt="img" width="800"/>

|  | Fst | weighted Fst |
| -----| ----| --------|
| ch x ref (unfold)| 0.005270 | 0.002872 |

1k (top), 5k (middle), 15k (bottom) windowed Fst plot

<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_unfold.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/Mahattan_ch_ref_single_SNP_unfold.jpg" alt="img" width="800"/>

### Questions

1. Following the suggestion by [Mikhail V Matz](https://github.com/z0on), perhaps removing the sites that are potentially from lumped paralogs in the reference assembly, see angsd issue here [Calculating FST: when to fold SFS · Issue #259 · ANGSD/angsd · GitHub](https://github.com/ANGSD/angsd/issues/259)

2. What if using folded option for SFS calculation, given we do not have real ancester genome (mentioned by https://github.com/ANGSD/angsd/issues/259)?

3. What does -whichFst 1 do?
