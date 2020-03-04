---
comments: true
title: angsd_SFS_FST
date: '2020-03-03 12:00'
tags:
  - angsd
  - SFS
  - fst
  - paralogs
  - WGS
categories:
  - WGS data analysis
---

In the last few days I was trying to figure out the best ways of generating the fst Manhattan plots for DelBay19 samples and identify the potential inversion regions.

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

# for wild populations

module load angsd/0.931

angsd -b ARN_44.bamlist -anc cv30.fa -out Fst_wild_mph/ARN_no44inv_minI12D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12 

angsd -b COH_42.bamlist -anc cv30.fa -out Fst_wild_mph/COH_no42inv_minI12D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12 

angsd -b HC_44.bamlist -anc cv30.fa -out Fst_wild_mph/HC_no44inv_minI12D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12 

angsd -b NB_47.bamlist -anc cv30.fa -out Fst_wild_mph/NB_no47inv_minI12D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12 

angsd -b SR_48.bamlist -anc cv30.fa -out Fst_wild_mph/SR_no48inv_minI12D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12 

# for challenge populations

module load angsd/0.931

angsd -b ch1_11.bamlist -anc cv30.fa -out Fst_ch/ch1_no11inv_minI2D5maxD12_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 5 -setMaxDepth 12 -minInd 2 

angsd -b ch2_11.bamlist -anc cv30.fa -out Fst_ch/ch2_no11inv_minI2D5maxD12_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 5 -setMaxDepth 12 -minInd 2 

angsd -b ch3_12.bamlist -anc cv30.fa -out Fst_ch/ch3_no12inv_minI2D5maxD12_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 5 -setMaxDepth 12 -minInd 2 

angsd -b ch4_14.bamlist -anc cv30.fa -out Fst_ch/ch4_no14inv_minI2D5maxD12_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 5 -setMaxDepth 12 -minInd 2

angsd -b ch_48.bamlist -anc cv30.fa -out Fst_ch/ch_no48inv_minI2D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12

angsd -b r_47.bamlist -anc cv30.fa -out Fst_ch/ref_no47inv_minI2D20maxD50_MQ20_minMAF05_SNPe6 -dosaf 1 -GL 1 -doGlf 2 -doMaf 1 -doMajorMinor 1 -doPost 1 -doVcf 1 -doCounts 1 -doDepth 1 -dumpCounts 1 -P 10 -minQ 20 -minMapQ 20 -setMinDepth 20 -setMaxDepth 50 -minInd 12

```	
### Get FST
```	
# Notes from Matt

# modify get_fst script for window sizes you want

./get_fst_windows.sh /workdir/mph75/angsd_mq20 /workdir/mph75/3pops_names.txt 1 _minI30D30maxD100_MQ20_minMAF05_SNPe6_no56inv >& Fst_MQ20_log &

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

awk -F"\t" '{print $2, $3, $3, $5}' HH_SV_minI30D30maxD100_MQ20_minMAF05_SNPe6_no56inv.15kb_win_15kb_fst | awk '$2-=7500' | awk '$3+=7500' > HH_SV_angsd_15K.fst
then edit headers manually: chr window_start window_end angsd_Fst
awk -F"\t" '{print $2, $3, $3, $5}' FIS_SV_minI30D30maxD100_MQ20_minMAF05_SNPe6_no56inv.15kb_win_15kb_fst | awk '$2-=7500' | awk '$3+=7500' > FIS_SV_angsd_15K.fst
awk -F"\t" '{print $2, $3, $3, $5}' FIS_HH_minI30D30maxD100_MQ20_minMAF05_SNPe6_no56inv.15kb_win_15kb_fst | awk '$2-=7500' | awk '$3+=7500' > FIS_HH_angsd_15K.fst

```
### qqman Manhattan plotting

NOTE - for qqman Manhattan plotting the chr need to be integer 

Plot with Manhattan_loop.R (several loop options in script)

plotting by chromosome is helpful for initial screening for inversions and other interesting patterns

Not sure about the most efficient way to calculate an overall Fst for a set of populations

### Questions

1. What does -whichFst 1 do?

2. Following the suggestion by [Mikhail V Matz](https://github.com/z0on), perhaps removing the sites that are potentially from lumped paralogs in the reference assembly, see angsd issue here [Calculating FST: when to fold SFS · Issue #259 · ANGSD/angsd · GitHub](https://github.com/ANGSD/angsd/issues/259)

3. What if using folded option for SFS calculation, given we do not have real ancester genome (mentioned by https://github.com/ANGSD/angsd/issues/259)?
