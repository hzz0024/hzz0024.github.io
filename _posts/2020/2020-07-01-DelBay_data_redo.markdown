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

In this post I am trying to rebuild the DelBay19 data for downstream deltap and outlier analyses. This is inspired by that the major and minor alleles are sometimes reversed if I run ref and challenge population independently. Instead, I need to use a inferred snp list to set the major and minor alleles. That is, the major and minor alleles were inferred from genotype likelihoods across all individuals (run ‐doMajorMinor 1 for challenge and wild populations separately), and these were then set for downstream analyses (‐doMajorMinor 3).

The angsd pipeline from [Claire Mérot's github](https://github.com/clairemerot) is well-organized and could be used as a good reference for this redo process. 

Besides, I'd like to compare the Fst, theta, and/or outlier identification results before/after genome masking. Note that steps below used both original and masked genome for data rebuilding. 

### 1) set up the 01_config.sh

```sh
#path to bam list
CH="/scratch/hzz0024/DelBay19_HG/02_info/ch_50.list"
REF="/scratch/hzz0024/DelBay19_HG/02_info/ref_48.list"
CHR="/scratch/hzz0024/DelBay19_HG/02_info/challenge_98.list"
WILD="/scratch/hzz0024/DelBay19_HG/02_info/wild_235.list"
HC="/scratch/hzz0024/DelBay19_HG/02_info/HC_48.list"
ARN="/scratch/hzz0024/DelBay19_HG/02_info/ARN_47.list"
COH="/scratch/hzz0024/DelBay19_HG/02_info/COH_44.list"
SR="/scratch/hzz0024/DelBay19_HG/02_info/SR_48.list"
NB="/scratch/hzz0024/DelBay19_HG/02_info/NB_48.list"
HC_ARN="/scratch/hzz0024/DelBay19_HG/02_info/HC_ARN_95.list"
ARN_COH="/scratch/hzz0024/DelBay19_HG/02_info/ARN_COH_91.list"
COH_SR="/scratch/hzz0024/DelBay19_HG/02_info/COH_SR_92.list"
SR_NB="/scratch/hzz0024/DelBay19_HG/02_info/SR_NB_96.list"

#path to the anc genome
ANC="/scratch/hzz0024/DelBay19_HG/genome/cv30.fa"
ANC_MASKED="/scratch/hzz0024/DelBay19_HG/genome/cv30_masked.fasta"

#filter : will keep SNP above this allele frequency (over all individuals)
MIN_MAF=0.05

#filter : will keep SNP with at least one read for this percentage of individuals (over all individuals in step 03, and within each pop at step 07)
PERCENT_IND=0.7

#filter: will keep SNP with at least a coverage of this factor multiplied by the number of ind - across all ind. usually set 2-4
#times the coverage to remove repeated regions
#MAX_DEPTH_FACTOR=3

#window size for sliding window FST & Thetas
WINDOW=1000

#window step
WINDOW_STEP=1000

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
for pop in CHR WILD
do
    echo -e 'module load angsd/0.931\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_list.txt" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\necho "Ouput of -doDepth can be used for depth evaluation with all individuals listed in "$'$pop'\necho "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"\necho "filter on allele frequency = "$MIN_MAF\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -anc $ANC_MASKED -remove_bads 1 -doQsDist 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -minMaf $MIN_MAF -maxDepth 2000 -b $'$pop' -out "/scratch/hzz0024/DelBay19_HG/03_depth/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mask"' >> formal/'02_depth_'$pop'_mask.sh'
    echo -e 'module load angsd/0.931\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_list.txt" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\necho "Ouput of -doDepth can be used for depth evaluation with all individuals listed in "$'$pop'\necho "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"\necho "filter on allele frequency = "$MIN_MAF\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -anc $ANC -remove_bads 1 -doQsDist 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -minMaf $MIN_MAF -maxDepth 2000 -b $'$pop' -out "/scratch/hzz0024/DelBay19_HG/03_depth/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"' >> formal/'02_depth_'$pop'_cv30.sh'
done

# An example script is shown below,

module load angsd/0.931
# maybe edit
target="CHR"
NB_CPU=20 #change accordingly
REGIONS="-rf chr_list.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh
N_IND=$(wc -l $CHR | cut -d " " -f 1)
MIN_IND=$(($N_IND*7/10))

echo "Ouput of -doDepth can be used for depth evaluation with all individuals listed in "$CHR
echo "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"
echo "filter on allele frequency = "$MIN_MAF

angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -anc $ANC -remove_bads 1 -doQsDist 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -minMaf $MIN_MAF -maxDepth 2000 -b $CHR -out "/scratch/hzz0024/DelBay19_HG/03_depth/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"
```

- Output 

CHR walltime: 17420 s       
WILD walltime: 38991 s       

|Group|Mean| SD| Mean+3SD|
|-----|----|---|---------|
|CHR  |205 |60 | 384     |
|WILD |353 |203| 962     |

### 4) create maf, saf, and beagle files for challenge and wild datasets

```sh
#!/bin/bash
# write up the script for 04_saf_maf_gl_all
# change the WILD into CHR and -setMaxDepth 962 into -setMaxDepth 384 for CHR group, see example script
# note for WILD group, 200GB memory is needed due to large number of samples
for pop in WILD
do
    echo -e 'module load angsd/0.931\n###this script will work on bamfiles by population and calculate saf & maf\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_list.txt" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\necho "Ouput can be used for depth evaluation with all individuals listed in "$'$pop'\necho "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"\necho "filter on allele frequency = "$MIN_MAF\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 2 -anc $ANC_MASKED -remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth 962 -b $'$pop' -out "/scratch/hzz0024/DelBay19_HG/04_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mask"\n\n#main features\n# -P nb of threads\n# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 1 samtools method - export GL in beagle format  -doGLF2)\n# -doMajorMinor 1 use the most frequent allele as major\n# -anc provide a ancestral sequence = reference in our case\n# -rf (file with the region written) work on a defined region : OPTIONAL\n# -b (bamlist) input file\n# -out  output file\n\n#main filters\n#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ minimum quality of reads\n#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 70%\n#filter on allele frequency -minMaf, set to 0.05\n\n#extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles\n#order the sites file by chromosome names\n#makes a region file matching the sites files and with same order\n#index sites file\necho "from the maf file, extract a list of SNP chr, positoin, major all, minor all"\ncd /scratch/hzz0024/DelBay19_HG/04_saf_maf_gl_all\nzcat "$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_mask.mafs.gz | tail -n +2 > FILE_mask.tmp && mv FILE_mask.tmp "$target"_snplist_mask\nawk '\''{print $1,$2,$3,$4}'\'' "$target"_snplist_mask > "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_mask\nangsd sites index "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_mask' >> formal/'04_saf_maf_al_all_'$pop'_mask.sh'
    echo -e 'module load angsd/0.931\n###this script will work on bamfiles by population and calculate saf  & maf\n# maybe edit\ntarget="'$pop'"\nNB_CPU=20 #change accordingly\nREGIONS="-rf chr_list.txt" #optional\n#REGIONS="" # to remove the options to focus on a limited number of regions\n\n#prepare variables - avoid to modify\nsource /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh\nN_IND=$(wc -l $'$pop' | cut -d " " -f 1)\nMIN_IND=$(($N_IND*7/10))\n\necho "Ouput can be used for depth evaluation with all individuals listed in "$'$pop'\necho "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"\necho "filter on allele frequency = "$MIN_MAF\n\nangsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 2 -anc $ANC -remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth 962 -b $'$pop' -out "/scratch/hzz0024/DelBay19_HG/04_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"\n\n#main features\n# -P nb of threads\n# -doMaf 1 (allele frequencies)  -dosaf (prior for SFS) -GL (Genotype likelihood 1 samtools method - export GL in beagle format  -doGLF2)\n# -doMajorMinor 1 use the most frequent allele as major\n# -anc provide a ancestral sequence = reference in our case\n# -rf (file with the region written) work on a defined region : OPTIONAL\n# -b (bamlist) input file\n# -out  output file\n\n#main filters\n#filter on bam files -remove_bads (remove files with flag above 255) -minMapQ minimum mapquality -minQ minimum quality of reads\n#filter on frequency -minInd (minimum number of individuals with at least one read at this locus) we set it to 70%\n#filter on allele frequency -minMaf, set to 0.05\n\n#extract SNP which passed the MIN_MAF and PERCENT_IND filters & their Major-minor alleles\n#order the sites file by chromosome names\n#makes a region file matching the sites files and with same order\n#index sites file\necho "from the maf file, extract a list of SNP chr, positoin, major all, minor all"\ncd /scratch/hzz0024/DelBay19_HG/04_saf_maf_gl_all\nzcat "$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30.mafs.gz | tail -n +2 > FILE_cv30.tmp && mv FILE_cv30.tmp "$target"_snplist_cv30\nawk '\''{print $1,$2,$3,$4}'\'' "$target"_snplist_cv30 > "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_cv30\nangsd sites index "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_cv30' >> formal/'04_saf_maf_al_all_'$pop'_cv30.sh'
done

# An example script is shown below,

module load angsd/0.931
###this script will work on bamfiles by population and calculate saf & maf
# maybe edit
target="CHR"
NB_CPU=20 #change accordingly
REGIONS="-rf chr_list.txt" #optional
#REGIONS="" # to remove the options to focus on a limited number of regions

#prepare variables - avoid to modify
source /scratch/hzz0024/DelBay19_HG/01_scripts/01_config.sh
N_IND=$(wc -l $CHR | cut -d " " -f 1)
MIN_IND=$(($N_IND*7/10))

echo "Ouput can be used for depth evaluation with all individuals listed in "$CHR
echo "keep loci with at leat one read for n individuals = "$MIN_IND", which is 70% of total "$N_IND" individuals"
echo "filter on allele frequency = "$MIN_MAF

angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doGlf 2 -doMajorMinor 1 -doCounts 1 -doDepth 1 -maxDepth 1000 -dumpCounts 2 -anc $ANC -remove_bads 1 -minMapQ 30 -minQ 20 -minInd $MIN_IND -minMaf $MIN_MAF -setMaxDepth 384 -b $CHR -out "/scratch/hzz0024/DelBay19_HG/04_saf_maf_gl_all/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"

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
cd /scratch/hzz0024/DelBay19_HG/04_saf_maf_gl_all
zcat "$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30.mafs.gz | tail -n +2 > FILE_cv30.tmp && mv FILE_cv30.tmp "$target"_snplist_cv30
awk '{print $1,$2,$3,$4}' "$target"_snplist_cv30 > "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_cv30
angsd sites index "$target"_sites_all_maf"$MIN_MAF"_pctind"$PERCENT_IND"_maxdepth3dv_snplist_4col_cv30
```

- Output 

CHR (cv30/mask) walltime: 16370/16081 s     
WILD (cv30/mask) walltime:      

No. of SNPs for each group

|           |CHR_no_mask |CHR_mask|WILD_no_mask|WILD_mask|
|-----------|------------|--------|------------|---------|
|NC_035780.1|246545      |244926  |
|NC_035781.1|263368      |260829  |
|NC_035782.1|282895      |280797  |
|NC_035783.1|253018      |250925  |
|NC_035784.1|464025      |462615  |
|NC_035785.1|68069       |67204   |
|NC_035786.1|68417       |67166   |
|NC_035787.1|96490       |94616   |
|NC_035788.1|118175      |115916  |
|NC_035789.1|24318       |23751   |
|Total      | 1885320    |1868745 |

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



