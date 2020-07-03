---
comments: true
title: Rebuild DelBay19 data 
date: '2020-07-01 12:00'
tags:
  - DelBay19
  - deltap
  - ouliter
  - null model
  - maf
  - WGS
categories:
  - WGS data analysis
---

In this post I am trying to rebuild the DelBay19 data for downstream deltap and outlier analyses. This is inspired by that the major and minor alleles are sometimes reversed if I run ref and challenge population independently. Instead, I need to use a inferred snp list to set the major and minor alleles. That is, major and minor alleles were inferred from genotype likelihoods across all individuals (run ‐doMajorMinor 1 for challenge and wild populations separately), and these were then set for downstream analyses (‐doMajorMinor 3).

The angsd pipeline from [Claire Mérot's github](https://github.com/clairemerot) is well-organized and could be used as a good reference for this redo process. 

Besides, I'd like to compare the Fst, theta, and/or outlier identification results before/after genome masking. 

### 1) set up the 01_config.sh

```sh
#path to bam list
CH=/scratch/hzz0024/DelBay19_HG/02_info/ch_50.list
REF=/scratch/hzz0024/DelBay19_HG/02_info/ref_48.list
CHR=/scratch/hzz0024/DelBay19_HG/02_info/challenge_98.list
WILD=/scratch/hzz0024/DelBay19_HG/02_info/wild_235.list
HC=/scratch/hzz0024/DelBay19_HG/02_info/HC_48.list
ARN=/scratch/hzz0024/DelBay19_HG/02_info/ARN_47.list
COH=/scratch/hzz0024/DelBay19_HG/02_info/COH_44.list
SR=/scratch/hzz0024/DelBay19_HG/02_info/SR_48.list
NB=/scratch/hzz0024/DelBay19_HG/02_info/NB_48.list
HC_ARN=/scratch/hzz0024/DelBay19_HG/02_info/HC_ARN_95.list
ARN_COH=/scratch/hzz0024/DelBay19_HG/02_info/ARN_COH_91.list
COH_SR=/scratch/hzz0024/DelBay19_HG/02_info/COH_SR_92.list
SR_NB=/scratch/hzz0024/DelBay19_HG/02_info/SR_NB_96.list

#path to the anc genome
ANC=/scratch/hzz0024/DelBay19_HG/genome/cv30.fa
ANC_MASKED=/scratch/hzz0024/DelBay19_HG/genome/cv30_masked.fasta

#path to bam folder
BAM_PATH=../run

#filter : will keep SNP above this allele frequency (over all individuals)
MIN_MAF=0.05

#filter : will keep SNP with at least one read for this percentage of individuals (over all individuals in step 03, and within each pop at step 07)
PERCENT_IND=0.7

#window size for sliding window FST & Thetas
WINDOW=1000

#window step
WINDOW_STEP=1000

#min nb of pop to consider for NGS admix
K_MIN=2

#maximum nb of pop to consider for NGS admix
K_MAX=5
```

### 2) set the maxdepth for challenge and wild datasets (02_depth_CHR.sh and 02_depth_WILD.sh)

```sh
#!/bin/bash
module load angsd/0.931

###this script will work on all bamfiles and calculate saf, maf & genotype likelihood
#maybe edit
# need to change $CHR at first echo part (L19) and angsd command (L29)
# also change the target below
target=CHR
NB_CPU=20 #change accordingly
REGIONS="-rf chr_list.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh
N_IND=$(wc -l $CHR | cut -d " " -f 1)
MIN_IND_FLOAT=$(echo "($N_IND * 0.7)"| bc -l)
MIN_IND=${MIN_IND_FLOAT%.*}

echo "Ouput of -doQsDist can be used for depth evaluation with  all individuals listed in $CHR"
echo "keep loci with at leat one read for n individuals = $MIN_IND, which is 70% of total $N_IND individuals"
echo "filter on allele frequency = $MIN_MAF"

####prepare the output for depth evaluation
angsd -P $NB_CPU \
-doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 \
-doDepth 1 -dumpCounts 1 \
-anc $ANC_MASKED -remove_bads 1 -doQsDist 1 -minMapQ 30 -minQ 20 \
-minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth 2000 \
-b $CHR -out /scratch/hzz0024/DelBay19_HG/03_depth/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"

#main features
# -P nb of threads
# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 1 samtools method - export GL in beagle format  -doGLF2)
# -doMajorMinor 1 use the most frequent allele as major
# -anc provide a ancestral sequence = reference in our case
# -rf (file with the region written) work on a defined region : OPTIONAL
# -b (bamlist) input file
# -out  output file

#main filters
#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ (minimum quality of reads?)
#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 70%
#filter on allele frequency -minMaf, set to 0.05
```


