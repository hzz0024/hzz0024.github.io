---
comments: true
title: DelBay19 theta estimate correction
date: '2020-07-16 12:00'
tags:
  - DelBay19
  - deltap
  - ouliter
  - null model
  - theta
  - drift
  - WGS
categories:
  - WGS data analysis
---

Link below stores the code for pi correction borrowed from Anna's github, I need to do some edits to make it useful for my DelBay19 data.

[https://github.com/atigano/Peromyscus_eremicus_genome/tree/master/variant_calling_ANGSD](https://github.com/atigano/Peromyscus_eremicus_genome/tree/master/variant_calling_ANGSD)


### Part I. Theta correction commands

- Files needed for this analyses:

```sh
angsd_peer_global_noout.sh - global calling

angsd -b $HOME/wgr_peromyscus/sample_lists/peer_goodbam_all_noout.txt \
-anc ~/genomes/Peer/Peer1.7.2.fasta \
-out $HOME/wgr_peromyscus/angsd_results_peer/angsd_peer_global_noout \
-P 24 -minMapQ 20 -minQ 20 -setMinDepth 13 -minInd 13 -GL 1 -minMaf 0.01 -SNP_pval 1e-6 \
-doMaf 1 -doMajorMinor 1 -doPost 1 -doCounts 1 -doDepth 1 -doGlf 3 -doGeno 32 -dumpCounts 1 -doSaf 1 -fold 1
```

output: angsd_peer_global_noout_sites.txt - variants list from the global calling with stringent quality filter

```sh
angsd_peer_global_noout_allvar.sh - variant and invariant calling

angsd -b $HOME/wgr_peromyscus/sample_lists/peer_goodbam_all_noout.txt \
-anc ~/genomes/Peer/Peer1.7.2.fasta \
-out $HOME/wgr_peromyscus/angsd_results_peer/angsd_peer_global_noout_allvar \
-P 24 -minMapQ 20 -minQ 20 -setMinDepth 13 -minInd 13 -GL 1 \
-doMaf 1 -doMajorMinor 1 -doPost 1 -doCounts 1 -doDepth 1
```

output: angsd_peer_global_noout_allvar.mafs - including all sites in the analysis (i.e. variant and invariant), need to to called again with angsd_peer_global_noout_allvar.sh

```sh
angsd -b ~/wgr_peromyscus/sample_lists/peer_goodbam_all_noout.txt -anc ~/genomes/Peer/Peer1.7.2.fasta -out ~/wgr_peromyscus/angsd_results_peer/angsd_peer_global_noout_folded.sfs -P 24 -doThetas 1 -doSaf 1 -pest angsd_peer_global_noout_folded.sfs -GL 1 -fold 1 -sites angsd_peer_global_noout_sites.txt
```

output: angsd_peer_global_noout_folded.sfs.thetas.site - file contailing pi for each loci (i.e. tsv file in theta estimation output)

1) step 1. split the reference genome into windows

```sh
### original code
WIN=$1 ###win size in kb

module load linuxbrew/colsa
###out of the loop
###make windows
bedtools makewindows -g ~/genomes/Peer/Peer1.7.2.fasta.fai -w ${WIN}000 | awk '$3 ~ /000$/' | sed 's/ /\t/g'> genome_windows_${WIN}k.bed

##################### for DelBay19 run ##################### 
module load bedtools/2.29.0
WIN=200
# make windows
bedtools makewindows -g cv30.fa.fai -w ${WIN} | awk '$3 ~ /000$/' | sed 's/ /\t/g'> genome_windows_${WIN}.bed
```

2) make bed file from the theta file containing pi

```sh
### original code
### make bed file from the theta file containing pi
tail -n +2 angsd_peer_global_noout_folded.sfs.thetas.site > angsd_peer_global_noout_folded.sfs.thetas.site_nohead
cut -f2 angsd_peer_global_noout_folded.sfs.thetas.site_nohead | awk '{$1 = $1 + 1; print}' | paste angsd_peer_global_noout_folded.sfs.thetas.site_nohead - | awk 'BEGIN {FS="\t"};{ print $1"\t"$2"\t"$8"\t"$4}' | sed 's/ //g' > pi_peer_global_noout.bed

##################### for DelBay19 run #####################  
tail -n +2 CHR_maf0.05_pctind0.7_cv30.thetas.tsv > CHR_maf0.05_pctind0.7_cv30.thetas.tsv.tsv_nohead
cut -f2 CHR_maf0.05_pctind0.7_cv30.thetas.tsv.tsv_nohead | awk '{$1 = $1 + 1; print}' | paste CHR_maf0.05_pctind0.7_cv30.thetas.tsv.tsv_nohead - | awk 'BEGIN {FS="\t"};{ print $1"\t"$2"\t"$8"\t"$4}' | sed 's/ //g' > CHR_pi_global.bed


head ref_pi_global.bed
NC_035780.1	1466	1467	-2.030720
NC_035780.1	1472	1473	-1.983084
NC_035780.1	4241	4242	-0.797251
NC_035780.1	4277	4278	-1.473429
NC_035780.1	4928	4929	-0.683224
NC_035780.1	4957	4958	-0.683153
NC_035780.1	5102	5103	-1.823825
NC_035780.1	5123	5124	-1.403415
NC_035780.1	5184	5185	-0.826394
NC_035780.1	5200	5201	-0.870824
```

3) make bed file for all variant and invariant sites and for each chromosome

```sh
###in the loop
for CHR in `cat chromosomes.txt`; do ###list of chromosomes in chromosomes.txt file; do
###make bed file for all variant and invariant sites for each chromosome
grep "$CHR" angsd_peer_global_noout_allvar.mafs > angsd_peer_global_noout_allvar_${CHR}.mafs
cut -f 1,2 angsd_peer_global_noout_allvar_${CHR}.mafs > peer_global_noout_allvar.sites_${CHR}.txt
cut -f2 peer_global_noout_allvar.sites_${CHR}.txt | awk '{$1 = $1 + 1; print}' | paste peer_global_noout_allvar.sites_${CHR}.txt - | sed 's/ //g'> peer_global_noout_allvar.sites_${CHR}.bed

###split the genome window file in chr
grep "$CHR" genome_windows_${WIN}k.bed > genome_windows_${WIN}k_${CHR}.bed

###calculate the number of sites in each window for each chromosome
bedtools coverage -a genome_windows_${WIN}k_${CHR}.bed -b peer_global_noout_allvar.sites_${CHR}.bed -counts > allsites_${WIN}kbwin_${CHR}.txt
cut -f4 allsites_${WIN}kbwin_${CHR}.txt | sed 's/\t0/NA/g' > allsites_${WIN}kwin_NA_${CHR}.txt

###split the per site theta file
grep "$CHR" pi_peer_global_noout.bed > pi_peer_global_noout_${CHR}.bed

awk '{print exp($4)}' pi_peer_global_noout_${CHR}.bed | paste pi_peer_global_noout_${CHR}.bed - > pi_peer_global_noout_log_${CHR}.bed

bedtools map -a genome_windows_${WIN}k_${CHR}.bed -b pi_peer_global_noout_log_${CHR}.bed -c 5 -o sum | sed 's/\t[.]/\tNA/g' - > pi_peer_global_noout_log_${WIN}kbwin_${CHR}.txt
###pi_peer_global_noout_log_50kbwin_chr2_pilon.txt 

paste pi_peer_global_noout_log_${WIN}kbwin_${CHR}.txt allsites_${WIN}kwin_NA_${CHR}.txt | sed 's/[.]\t/NA\t/g' - > pi_peer_global_noout_log_${WIN}kbwin_sites_${CHR}.txt

awk '{if(/NA/)var="NA";else var=$4/$5;print var}' pi_peer_global_noout_log_${WIN}kbwin_sites_${CHR}.txt | paste pi_peer_global_noout_log_${WIN}kbwin_sites_${CHR}.txt - > pi_peer_global_noout_log_${WIN}kbwin_sites_corrected_${CHR}.txt

done
##################### for DelBay19 run ##################### 

###in the loop
for CHR in `cat chromosomes.txt`; do ###list of chromosomes in chromosomes.txt file; do
###make bed file for all variant and invariant sites for each chromosome
grep "$CHR" CHR_maf0.05_pctind0.7_cv30_allvar.mafs > CHR_maf0.05_pctind0.7_cv30_allvar_${CHR}.mafs
cut -f 1,2 CHR_maf0.05_pctind0.7_cv30_allvar_${CHR}.mafs > CHR_maf0.05_pctind0.7_cv30_allvar_${CHR}.txt
cut -f2 CHR_maf0.05_pctind0.7_cv30_allvar_${CHR}.txt | awk '{$1 = $1 + 1; print}' | paste CHR_maf0.05_pctind0.7_cv30_allvar_${CHR}.txt - | sed 's/ //g'> CHR_maf0.05_pctind0.7_cv30_allvar_${CHR}.bed

###split the genome window file in chr
grep "$CHR" genome_windows_${WIN}.bed > genome_windows_${WIN}_${CHR}.bed


###calculate the number of sites in each window for each chromosome
bedtools coverage -a genome_windows_${WIN}_${CHR}.bed -b CHR_maf0.05_pctind0.7_cv30_allvar_${CHR}.bed -counts > allvar_${WIN}bwin_${CHR}.txt
## not sure why replace the \t0 here
cut -f4 allvar_${WIN}bwin_${CHR}.txt | sed 's/\t0/NA/g' > allvar_${WIN}win_NA_${CHR}.txt

grep "$CHR" CHR_pi_global.bed > CHR_pi_global_${CHR}.bed
## pate - will add the new column to the end of data
awk '{print exp($4)}' CHR_pi_global_${CHR}.bed | paste CHR_pi_global_${CHR}.bed - > CHR_pi_global_log_${CHR}.bed
# bedtools is used to sum up the theta value in each window
bedtools map -a genome_windows_${WIN}_${CHR}.bed -b CHR_pi_global_log_${CHR}.bed -c 5 -o sum | sed 's/\t[.]/\tNA/g' - > CHR_pi_global_log_${WIN}bwin_${CHR}.txt
###pi_peer_global_noout_log_50kbwin_chr2_pilon.txt

paste CHR_pi_global_log_${WIN}bwin_${CHR}.txt allvar_${WIN}win_NA_${CHR}.txt | sed 's/[.]\t/NA\t/g' - > pi_peer_global_noout_log_${WIN}bwin_sites_${CHR}.txt
## divide the theta by number of SNPs in a window
awk '{if(/NA/)var="NA";else var=$4/$5;print var}' pi_peer_global_noout_log_${WIN}bwin_sites_${CHR}.txt | paste pi_peer_global_noout_log_${WIN}bwin_sites_${CHR}.txt - > pi_peer_global_noout_log_${WIN}bwin_sites_corrected_${CHR}.txt

done
```

An example showing the corrected theta (last column)

```sh
head pi_peer_global_noout_log_200bwin_sites_corrected_NC_035780.1.txt

NC_035780.1	0	200	NA	32	NA
NC_035780.1	200	400	NA	64	NA
NC_035780.1	400	600	NA	155	NA
NC_035780.1	600	800	NA	58	NA
NC_035780.1	800	1000	NA	90	NA
NC_035780.1	1000	1200	NA	0	NA
NC_035780.1	1200	1400	NA	0	NA
NC_035780.1	1400	1600	0.274983	34	0.00808774
NC_035780.1	1600	1800	NA	0	NA
NC_035780.1	1800	2000	NA	0	NA
```

### Part II. Pick the right window size using single chromosome as an example




