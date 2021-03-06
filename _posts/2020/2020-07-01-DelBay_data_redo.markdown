---
comments: true
title: Rebuild DelBay19 data 
date: '2020-07-08 12:00'
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

In this post I am trying to rebuild the DelBay19 data for downstream deltap and outlier analyses. This is inspired by that the major and minor alleles are sometimes reversed if I run ref and challenge population independently. Instead, I need to use a inferred snp list to set the major and minor alleles. That is, the major and minor alleles were inferred from genotype likelihoods across all individuals, and these were then set for downstream analyses (‐doMajorMinor 3).

The angsd pipeline from [Claire Mérot's github](https://github.com/clairemerot) is well-organized and could be used as a good reference for this redo process. 

Besides, I'd like to compare the Fst, theta, and/or outlier identification results before/after genome masking. Note that steps below used both original and masked genome for data rebuilding. 

### 1) set up the 01_config.sh

```sh
#path to bam list
CH="/scratch/hzz0024/DelBay19_July/02_info/ch_50.list"
REF="/scratch/hzz0024/DelBay19_July/02_info/ref_48.list"
CHR="/scratch/hzz0024/DelBay19_July/02_info/challenge_98.list"
WILD="/scratch/hzz0024/DelBay19_July/02_info/wild_235.list"
ALL="/scratch/hzz0024/DelBay19_July/02_info/ALL_333.list"
HC="/scratch/hzz0024/DelBay19_July/02_info/HC_48.list"
ARN="/scratch/hzz0024/DelBay19_July/02_info/ARN_47.list"
COH="/scratch/hzz0024/DelBay19_July/02_info/COH_44.list"
SR="/scratch/hzz0024/DelBay19_July/02_info/SR_48.list"
NB="/scratch/hzz0024/DelBay19_July/02_info/NB_48.list"
HC_ARN="/scratch/hzz0024/DelBay19_July/02_info/HC_ARN_95.list"
ARN_COH="/scratch/hzz0024/DelBay19_July/02_info/ARN_COH_91.list"
COH_SR="/scratch/hzz0024/DelBay19_July/02_info/COH_SR_92.list"
SR_NB="/scratch/hzz0024/DelBay19_July/02_info/SR_NB_96.list"
test="/scratch/hzz0024/DelBay19_July/02_info/test.list"

#path to the anc genome
ANC="/scratch/hzz0024/DelBay19_July/genome/cv30.fa"
ANC_MASKED="/scratch/hzz0024/DelBay19_July/genome/cv30_masked.fasta"

#path to bam folder
#BAM_PATH=../02_info

#path to pcaangsd
#PCA_ANGSD_PATH="/project/lbernatchez/programs/pcangsd "

#filter : will keep SNP above this allele frequency (over all individuals)
MIN_MAF=0.05

#filter : will keep SNP with at least one read for this percentage of individuals (over all individuals in step 03, and within each pop at step 07)
PERCENT_IND=0.7

#filter: will keep SNP with at least a coverage of this factor multiplied by the number of ind - across all ind. usually set 2-4
#times the coverage to remove repeated regions
#MAX_DEPTH_FACTOR=3

#window size for sliding window FST & Thetas
WINDOW=200

#window step
WINDOW_STEP=200

#min nb of pop to consider for NGS admix
K_MIN=2

#maximum nb of pop to consider for NGS admix
K_MAX=5
```

### 2) set up the bam list file in 02_info

create the bam file list for each population

### 3) set the maxdepth for challenge and wild datasets

```sh
#!/bin/bash
# write up the script for maxDepth evaluation
# first echo is for results from masked genome, second one is from original genome
for pop in CHR WILD ALL
do
    echo -e 'module load angsd/0.931\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_list.txt" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_July/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\necho "Ouput of -doDepth can be used for depth evaluation with all individuals listed in "$'$pop'\necho "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"\necho "filter on allele frequency = "$MIN_MAF\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -anc $ANC_MASKED -remove_bads 1 -doQsDist 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-6 -minInd $MIN_IND -minMaf $MIN_MAF -maxDepth 2000 -b $'$pop' -out "/scratch/hzz0024/DelBay19_July/03_depth/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mask"' >> formal/'02_depth_'$pop'_mask.sh'
    echo -e 'module load angsd/0.931\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_list.txt" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_July/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\necho "Ouput of -doDepth can be used for depth evaluation with all individuals listed in "$'$pop'\necho "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"\necho "filter on allele frequency = "$MIN_MAF\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -anc $ANC -remove_bads 1 -doQsDist 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-6 -minInd $MIN_IND -minMaf $MIN_MAF -maxDepth 2000 -b $'$pop' -out "/scratch/hzz0024/DelBay19_July/03_depth/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"' >> formal/'02_depth_'$pop'_cv30.sh'
done

# An example script is shown below,
module load angsd/0.931
# maybe edit
target="ALL"
NB_CPU=20 #change accordingly
REGIONS="-rf chr_list.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /scratch/hzz0024/DelBay19_July/01_scripts/01_config.sh
N_IND=$(wc -l $ALL | cut -d " " -f 1)
MIN_IND=$(($N_IND*7/10))

echo "Ouput of -doDepth can be used for depth evaluation with all individuals listed in "$ALL
echo "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"
echo "filter on allele frequency = "$MIN_MAF

angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -anc $ANC -remove_bads 1 -doQsDist 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-6 -minInd $MIN_IND -minMaf $MIN_MAF -maxDepth 2000 -b $ALL -out "/scratch/hzz0024/DelBay19_July/03_depth/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"
```

- Output 

CHR walltime: 17420 s       
WILD walltime: 38991 s       

|Group|Mean| SD| Mean+3SD|
|-----|----|---|---------|
|CHR  |205 |60 | 384     |
|WILD |353 |203| 962     |
|ALL  |777 |210| 1407    |

### 4) create maf, saf, and beagle files for challenge (n=98), wild (n=235) and ALL (n=333) datasets

```sh
#!/bin/bash
# Read a string with spaces using for loop
for pop in ALL
do
    echo -e 'module load angsd/0.931\n###this script will work on bamfiles by population and calculate saf & maf\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_list.txt" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_July/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\necho "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"\necho "filter on allele frequency = "$MIN_MAF\n\nangsd -P $NB_CPU -doMaf 1 -doVcf 1 -dosaf 1 -doPost 1 -doIBS 1 -doCov 1 -makematrix 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 4 -anc $ANC -remove_bads 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-6 -minInd $MIN_IND -minMaf $MIN_MAF $REGIONS -setMaxDepth 962 -b $'$pop' -out "/scratch/hzz0024/DelBay19_July/04_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 2 -anc $ANC -remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND $REGIONS -setMaxDepth 962 -b $'$pop' -out "/scratch/hzz0024/DelBay19_July/04_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30_allvar"\n#main features\n# -P nb of threads\n# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 1 samtools method - export GL in beagle format  -doGLF2)\n# -doMajorMinor 1 use the most frequent allele as major\n# -anc provide a ancestral sequence = reference in our case\n# -rf (file with the region written) work on a defined region : OPTIONAL\n# -b (bamlist) input file\n# -out  output file\n\n#main filters\n#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ minimum quality of reads\n#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 70%\n#filter on allele frequency -minMaf, set to 0.05\n\n#extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles\n#order the sites file by chromosome names\n#makes a region file matching the sites files and with same order\n#index sites file\necho "from the maf file, extract a list of SNP chr, positoin, major all, minor all"\ncd /scratch/hzz0024/DelBay19_July/04_saf_maf_gl_all\nzcat "$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30.mafs.gz | tail -n +2 > FILE_cv30.tmp && mv FILE_cv30.tmp "$target"_snplist_cv30\nawk '\''{print $1,$2,$3,$4}'\'' "$target"_snplist_cv30 > "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_cv30\nangsd sites index "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_cv30' >> formal/'04_saf_maf_al_all_'$pop'_cv30_noinvers.sh'

    echo -e 'module load angsd/0.931\n###this script will work on bamfiles by population and calculate saf & maf\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_list.txt" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_July/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\necho "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"\necho "filter on allele frequency = "$MIN_MAF\n\nangsd -P $NB_CPU -doMaf 1 -doVcf 1 -dosaf 1 -doPost 1 -doIBS 1 -doCov 1 -makematrix 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 4 -anc $ANC_MASKED -remove_bads 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-6 -minInd $MIN_IND -minMaf $MIN_MAF $REGIONS -setMaxDepth 962 -b $'$pop' -out "/scratch/hzz0024/DelBay19_July/04_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mask"\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 2 -anc $ANC_MASKED -remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND $REGIONS -setMaxDepth 962 -b $'$pop' -out "/scratch/hzz0024/DelBay19_July/04_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mask_allvar"\n#main features\n# -P nb of threads\n# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 1 samtools method - export GL in beagle format  -doGLF2)\n# -doMajorMinor 1 use the most frequent allele as major\n# -anc provide a ancestral sequence = reference in our case\n# -rf (file with the region written) work on a defined region : OPTIONAL\n# -b (bamlist) input file\n# -out  output file\n\n#main filters\n#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ minimum quality of reads\n#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 70%\n#filter on allele frequency -minMaf, set to 0.05\n\n#extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles\n#order the sites file by chromosome names\n#makes a region file matching the sites files and with same order\n#index sites file\necho "from the maf file, extract a list of SNP chr, positoin, major all, minor all"\ncd /scratch/hzz0024/DelBay19_July/04_saf_maf_gl_all\nzcat "$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mask.mafs.gz | tail -n +2 > FILE_mask.tmp && mv FILE_mask.tmp "$target"_snplist_mask\nawk '\''{print $1,$2,$3,$4}'\'' "$target"_snplist_mask > "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_mask\nangsd sites index "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_mask' >> formal/'04_saf_maf_al_all_'$pop'_mask_noinvers.sh'
done

# An example script is shown below,

module load angsd/0.931
###this script will work on bamfiles by population and calculate saf & maf
# maybe edit
target="ALL"
NB_CPU=20 #change accordingly
REGIONS="-rf chr_list.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /scratch/hzz0024/DelBay19_July/01_scripts/01_config.sh
N_IND=$(wc -l $ALL | cut -d " " -f 1)
MIN_IND=$(($N_IND*7/10))

echo "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"
echo "filter on allele frequency = "$MIN_MAF

angsd -P $NB_CPU -doMaf 1 -doVcf 1 -dosaf 1 -doPost 1 -doIBS 1 -doCov 1 -makematrix 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 4 -anc $ANC -remove_bads 1 -minMapQ 30 -minQ 20 -SNP_pval 1e-6 -minInd $MIN_IND -minMaf $MIN_MAF $REGIONS -setMaxDepth 962 -b $ALL -out "/scratch/hzz0024/DelBay19_July/04_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"
angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 2 -anc $ANC -remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND $REGIONS -setMaxDepth 962 -b $ALL -out "/scratch/hzz0024/DelBay19_July/04_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30_allvar"
#main features
# -P nb of threads
# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 1 samtools method - export GL in beagle format  -doGLF2)
# -doMajorMinor 1 use the most frequent allele as major
# -anc provide a ancestral sequence = reference in our case
# -rf (file with the region written) work on a defined region : OPTIONAL
# -b (bamlist) input file
# -out  output file

#main filters
#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ minimum quality of reads
#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 70%
#filter on allele frequency -minMaf, set to 0.05

#extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles
#order the sites file by chromosome names
#makes a region file matching the sites files and with same order
#index sites file
echo "from the maf file, extract a list of SNP chr, positoin, major all, minor all"
cd /scratch/hzz0024/DelBay19_July/04_saf_maf_gl_all
zcat "$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30.mafs.gz | tail -n +2 > FILE_cv30.tmp && mv FILE_cv30.tmp "$target"_snplist_cv30
awk '{print $1,$2,$3,$4}' "$target"_snplist_cv30 > "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_cv30
angsd sites index "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_cv30
```

- Output 

CHR (cv30/mask) walltime: 16370/16081 s     
WILD (cv30/mask) walltime:      

No. of SNPs for each group (before inversion removal)

|           |CHR_no_mask |CHR_mask|WILD_no_mask|WILD_mask|
|-----------|------------|--------|------------|---------|
|NC_035780.1|246460      |244926  |2293798     |291460   |
|NC_035781.1|263299      |260829  |315613      |312294   |
|NC_035782.1|282811      |280797  |337773      |334781   |
|NC_035783.1|252959      |250925  |304301      |301402   |
|NC_035784.1|463899      |462615  |555977      |553815   |
|NC_035785.1|68048       |67204   |80197       | 78960   |
|NC_035786.1|68401       |67166   |82698       |80886    |
|NC_035787.1|96466       |94616   |117951      |115320   |
|NC_035788.1|118142      |115916  |144945      |141895   |
|NC_035789.1|24316       |23751   |30442       |29684    |
|Total      |1884801     |1868745 |2263695     |2240497  |

No. of SNPs for each group (after inversion removal)

|           |CHR_no_mask |CHR_mask|WILD_no_mask|WILD_mask|
|-----------|------------|--------|------------|---------|
|NC_035780.1|246545      |244926  |293798      |291460   |
|NC_035781.1|263368      |260829  |315613      |312294   |
|NC_035782.1|282895      |280797  |337773      |334781   |
|NC_035783.1|253018      |250925  |304301      |301402   |
|NC_035784.1|366069      |364659  |439627      |437465   |
|NC_035785.1|13219       |12368   |16467       | 15258   |
|NC_035786.1|68417       |67166   |82698       |80886    |
|NC_035787.1|96490       |94616   |117951      |115320   |
|NC_035788.1|118175      |115916  |144945      |141895   |
|NC_035789.1|24318       |23751   |30442       |29684    |
|Total      |1732514     |1715953 |2083615     |2060445  |

### 5) create saf for Fst estimates

```sh
#!/bin/bash
# Change the -site ### for CHR and WILD group
for pop in CH REF
do
    echo -e 'module load angsd/0.931\n###this script will work on bamfiles by population and calculate saf  & maf\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_list.txt" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\necho "Ouput can be used for depth evaluation with all individuals listed in "$'$pop'\necho "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"\necho "filter on allele frequency = "$MIN_MAF\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doMajorMinor 3 -anc $ANC_MASKED -remove_bads 1 -minMapQ 30 -minQ 20 -b $'$pop' -sites CHR_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_mask -out "/scratch/hzz0024/DelBay19_HG/05_saf_maf_by_pop/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mask"' >> formal/'05_saf_maf_by_pop_'$pop'_mask.sh'
    echo -e 'module load angsd/0.931\n###this script will work on bamfiles by population and calculate saf  & maf\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_list.txt" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\necho "Ouput can be used for depth evaluation with all individuals listed in "$'$pop'\necho "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"\necho "filter on allele frequency = "$MIN_MAF\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doMajorMinor 3 -anc $ANC -remove_bads 1 -minMapQ 30 -minQ 20 -b $'$pop' -sites CHR_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 -out "/scratch/hzz0024/DelBay19_HG/05_saf_maf_by_pop/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"' >> formal/'05_saf_maf_by_pop_'$pop'_cv30.sh'
done

# An example script is shown below,

module load angsd/0.931
###this script will work on bamfiles by population and calculate saf  & maf
# maybe edit
target="CH"
NB_CPU=20 #change accordingly
REGIONS="-rf chr_list.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh
N_IND=$(wc -l $CH | cut -d " " -f 1)
MIN_IND=$(($N_IND*7/10))

echo "Ouput can be used for depth evaluation with all individuals listed in "$CH
echo "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"
echo "filter on allele frequency = "$MIN_MAF

angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doMajorMinor 3 -anc $ANC -remove_bads 1 -minMapQ 30 -minQ 20 -b $CH -sites CHR_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 -out "/scratch/hzz0024/DelBay19_HG/05_saf_maf_by_pop/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"
```

- Output 

CH (cv30/mask) walltime: 8856/8882 s     
REF (cv30/mask) walltime: 8168/8542 s

- Question here: is there any fixed allele frequency in the maf files?

```R
setwd("/DelBay19_HG/05_saf_maf_by_pop")
filename = 'CH_maf0.05_pctind0.7_cv30_nochr56invers.mafs'
dat <- read.delim(filename, header = TRUE, sep='\t')
> max_af=max(dat$knownEM)
> min_af=min(dat$knownEM)
> max_af
[1] 0.718242
> min_af
[1] 0
> filename = 'REF_maf0.05_pctind0.7_cv30_nochr56invers.mafs'
> dat <- read.delim(filename, header = TRUE, sep='\t')
> max_af=max(dat$knownEM)
> min_af=min(dat$knownEM)
> max_af
[1] 0.70511
> min_af
[1] 0
> filename = 'HC_maf0.05_pctind0.7_cv30_no56invers.mafs'
> dat <- read.delim(filename, header = TRUE, sep='\t')
> max_af=max(dat$knownEM)
> min_af=min(dat$knownEM)
> max_af
[1] 0.765922
> min_af
[1] 0
> filename = 'ARN_maf0.05_pctind0.7_cv30_no56invers.mafs'
> dat <- read.delim(filename, header = TRUE, sep='\t')
> max_af=max(dat$knownEM)
> min_af=min(dat$knownEM)
> max_af
[1] 0.766392
> min_af
[1] 0
> filename = 'COH_maf0.05_pctind0.7_cv30_no56invers.mafs'
> dat <- read.delim(filename, header = TRUE, sep='\t')
> max_af=max(dat$knownEM)
> min_af=min(dat$knownEM)
> max_af
[1] 0.832372
> min_af
[1] 0
> filename = 'SR_maf0.05_pctind0.7_cv30_no56invers.mafs'
> dat <- read.delim(filename, header = TRUE, sep='\t')
> max_af=max(dat$knownEM)
> min_af=min(dat$knownEM)
> max_af
[1] 0.774196
> min_af
[1] 0
> filename = 'NB_maf0.05_pctind0.7_cv30_no56invers.mafs'
> dat <- read.delim(filename, header = TRUE, sep='\t')
> max_af=max(dat$knownEM)
> min_af=min(dat$knownEM)
> max_af
[1] 0.754925
> min_af
[1] 0
> filename = 'maf_test.mafs'
> dat <- read.delim(filename, header = TRUE, sep='\t')
> max_af=max(dat$knownEM)
> min_af=min(dat$knownEM)
> max_af
[1] 1
> min_af
[1] 0
```

I performed a small test with 5 individual within the challenge population (CH) and examine the maf result. This small amount of sample sould produce some snp with allele frquency due to drift. 

```R
filename = 'maf_test.mafs'
dat <- read.delim(filename, header = TRUE, sep='\t')
max_maf=max(dat$knownEM)
max_maf
[1] 1
min_af=min(dat$knownEM)
min_af
[1] 0
```

I observed the fixed snp with allele frequency = 0, suggesting that Angsd did not show unexpected behanior for maf calculation.  

- Fst (weighted) pairwise comparsion using non-masked genome

|            |   HC   |  ARN   |  COH   |   SR   |   NB   |   CH   |   REF  | 
| -----------|--------|--------|--------|--------|--------|--------|--------|
|    HC      |   −    |0.000521|0.000475|0.000514|0.000732|    −   |    −   |
|    ARN     |0.000521|   −    |0.000433|0.000338|0.000530|    −   |    −   |
|    COH     |0.000475|0.000433|   −    |0.000539|0.000681|    −   |    −   |
|    SR      |0.000514|0.000338|0.000539|   −    |0.000533|    −   |    −   |
|    NB      |0.000732|0.000530|0.000681|0.000533|   −    |    −   |    −   |
|    CH      |   −    |   −    |   −    |   −    |   −    |    −   |0.000584|
|    REF     |   −    |   −    |   −    |   −    |   −    |0.000584|    −   |


- Fst (weighted) pairwise comparsion using masked genome

|            |   HC   |  ARN   |  COH   |   SR   |   NB   |   CH   |   REF  | 
| -----------|--------|--------|--------|--------|--------|--------|--------|
|    HC      |   −    |0.000521|0.000469|0.000513|0.000732|    −   |    −   |
|    ARN     |0.000521|   −    |0.000428|0.000334|0.000527|    −   |    −   |
|    COH     |0.000469|0.000428|   −    |0.000539|0.000682|    −   |    −   |
|    SR      |0.000513|0.000334|0.000539|   −    |0.000540|    −   |    −   |
|    NB      |0.000732|0.000527|0.000682|0.000540|   −    |    −   |    −   |
|    CH      |   −    |   −    |   −    |   −    |   −    |    −   |0.000585|
|    REF     |   −    |   −    |   −    |   −    |   −    |0.000585|    −   |

Compared to the non-masked genome results, Fst values from the masked genome have 5 decreases vs. 3 increases observations. May need further tests to see if the masked regions lead to such changes. 

### 6) conduct Fst estimates

```sh
#!/bin/bash
for genome in cv30 mask
do
    echo -e './get_fst.sh /scratch/hzz0024/DelBay19_HG/06_fst_by_pop_pair pop_ch.txt 1 _maf0.05_pctind0.7_'$genome'' >> formal/'06_fst_by_pop_pair_CHR_'$genome'.sh'
done

# An example script is shown below,

./get_fst.sh /scratch/hzz0024/DelBay19_HG/06_fst_by_pop_pair pop_ch.txt 1 _maf0.05_pctind0.7_cv30

```

- Output

walltime: ~ 1 hour for each run

### 7) Fst plots for challenge vs. reference

<img src="https://hzz0024.github.io/images/Fst/Mahattan_ch_ref_cv30.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/Fst/Mahattan_ch_ref_mask.jpg" alt="img" width="800"/>

I can only observe some subtle differences between the plots from non-masked (first plot) vs masked (second one) genomes, let us take a look at the single SNP plots between these two.

<img src="https://hzz0024.github.io/images/Fst/cv30_single.jpg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/Fst/mask_single.jpg" alt="img" width="800"/>

Note the 0.2 line does not mean anthing. It is just used for easier comparsion. By comparing plots below, it seems 

1) fst values for some potential outliers are increased with the masked genome but changes are subtle (for example, SNP marked with red arrow).     
2) SNPs may be removed due to the genome masking, see SNP marked with blue arrow.     

### 8) the SGS test is under running

```sh
# script for deltap calculation
Rscript deltaP_abs.R -d /Users/ryan/Documents/Ryan_workplace/DelBay19_HG/SGS -p CH_maf0.05_pctind0.7_cv30.mafs -q REF_maf0.05_pctind0.7_cv30.mafs -t 1885320 -o obs_deltap_cv30.output

Rscript deltaP_abs.R -d /Users/ryan/Documents/Ryan_workplace/DelBay19_HG/SGS -p CH_maf0.05_pctind0.7_mask.mafs -q REF_maf0.05_pctind0.7_mask.mafs -t 1868745 -o obs_deltap_mask.output
```

### 9) ngsLD

```sh
##### step 1 create the snp list for each group and chromosome 
for i in {0..9}
do
awk -v val=$i '(($1 == "NC_03578"val".1"))' CHR_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30_no56invers.snplist > "CHR_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30_no56invers_chr"$(($i+1))".list"
awk -v val=$i '(($1 == "NC_03578"val".1"))' WILD_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30_no56invers.snplist > "WILD_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30_no56invers_chr"$(($i+1))".list"
done
##### index the sitelist file using angsd
module load angsd/0.931
for i in {1..10}
do
angsd sites index "CHR_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_mask_no56invers_chr"$i".list"
done


# stpe 2 create beagle files for ngsLD calculation
# for challenge group
for i in {1..10}
do
    for pop in CH REF
    do
    echo -e 'module load angsd/0.931\n#this script is used to create beagle files for ngsLD calculation\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_'$i'.list" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 3 -anc $ANC_MASKED -remove_bads 1 -minMapQ 30 -minQ 20 -b $'$pop' -sites CHR_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30_no56invers_chr'$i'.list -out "/scratch/hzz0024/DelBay19_HG/07_ngsLD/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30_nochr56invers_ch'$i'"' >> formal/'07_ngsLD_'$pop'_cv30_no56invers_ch'$i'.sh'
    done
done
# for wild group
for i in {1..10}
do
    for pop in HC ARN COH SR NB
    do
    echo -e 'module load angsd/0.931\n#this script is used to create beagle files for ngsLD calculation\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_'$i'.list" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 3 -anc $ANC_MASKED -remove_bads 1 -minMapQ 30 -minQ 20 -b $'$pop' -sites WILD_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30_no56invers_chr'$i'.list -out "/scratch/hzz0024/DelBay19_HG/07_ngsLD/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30_nochr56invers_ch'$i'"' >> formal/'07_ngsLD_'$pop'_cv30_no56invers_ch'$i'.sh'
    done
done

# An example script is shown below,

module load angsd/0.931
#this script is used to create beagle files for ngsLD calculation
# maybe edit
target="CH"
NB_CPU=20 #change accordingly
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh
N_IND=$(wc -l $CH | cut -d " " -f 1)
MIN_IND=$(($N_IND*7/10))

angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 3 -anc $ANC_MASKED -remove_bads 1 -minMapQ 30 -minQ 20 -b $CH -sites CHR_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30_no56invers_chr1.list -out "/scratch/hzz0024/DelBay19_HG/07_ngsLD/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30_nochr56invers_ch1"

# step 3 run ngsLD

# create the position file
for i in {1..10};do
    zcat 'REF_maf0.05_pctind0.7_cv30_nochr56invers_ch'$i'.mafs.gz' | cut -f 1,2 | tail -n +2 > 'ref_chr'$i'_pos.txt'
done

# run ngsLD
module load perl/5.26.1
module load R/3.6.3
#module load gsl/intel/2.1
module load gcc/4.9.3
#module load zlib/1.2.8
module load curl/7.47.1
module load samtools/1.6
module load bzip2/1.6.0
module load intel/mpi/2017
module load zlib/1.2.8
module load gsl/intel/2.1

for i in {1..10};do

/tools/ngsld/ngsLD \
--geno 'REF_maf0.05_pctind0.7_cv30_nochr56invers_ch'$i'.beagle.gz' \
--pos 'ref_chr'$i'_pos.txt' \
--n_ind 48 \
--n_sites $(cat 'ref_chr'$i'_pos.txt' | wc -l) \
--out 'ref_chr'$i'.output' \
--probs \
--max_kb_dist  1 \
--min_maf 0.05 \
--n_threads 20

done
```

