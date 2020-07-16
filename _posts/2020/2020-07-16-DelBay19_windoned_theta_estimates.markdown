---
comments: true
title: DelBay19 theta estimates at window scale
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

- Files needed for this analyses:

```sh
angsd_peer_global_noout.sh - global calling

angsd -b $HOME/wgr_peromyscus/sample_lists/peer_goodbam_all_noout.txt \
-anc ~/genomes/Peer/Peer1.7.2.fasta \
-out $HOME/wgr_peromyscus/angsd_results_peer/angsd_peer_global_noout \
-P 24 -minMapQ 20 -minQ 20 -setMinDepth 13 -minInd 13 -GL 1 -minMaf 0.01 -SNP_pval 1e-6 \
-doMaf 1 -doMajorMinor 1 -doPost 1 -doCounts 1 -doDepth 1 -doGlf 3 -doGeno 32 -dumpCounts 1 -doSaf 1 -fold 1

angsd_peer_global_noout_allvar.sh - variant and invariant calling

angsd -b $HOME/wgr_peromyscus/sample_lists/peer_goodbam_all_noout.txt \
-anc ~/genomes/Peer/Peer1.7.2.fasta \
-out $HOME/wgr_peromyscus/angsd_results_peer/angsd_peer_global_noout_allvar \
-P 24 -minMapQ 20 -minQ 20 -setMinDepth 13 -minInd 13 -GL 1 \
-doMaf 1 -doMajorMinor 1 -doPost 1 -doCounts 1 -doDepth 1

angsd_peer_pop.sh - recall snp within each popultion

angsd -b $HOME/wgr_peromyscus/sample_lists/peer_goodbam_${POP}.txt \
-anc ~/genomes/Peer/Peer1.7.2.fasta \
-out $HOME/wgr_peromyscus/angsd_results_peer/angsd_peer_${POP} \
-P 24 -minMapQ 20 -minQ 20 -setMinDepth $NUM -minInd $NUM -GL 1 \
-doMaf 1 -doMajorMinor 3 -doPost 1 -doCounts 1 -doDepth 1 -doGlf 3  -doGeno 32 -dumpCounts 1 -doSaf 1 -fold 1 -sites angsd_peer_global_noout_sites.txt
```

angsd_peer_global_noout_sites.txt - variants list from the global calling with stringent quality filter     
angsd_peer_global_noout_allvar.mafs - including all sites in the analysis (i.e. variant and invariant), need to to called again with angsd_peer_global_noout_allvar.sh
angsd_peer_global_noout_folded.sfs.thetas.site - file contailing pi for each loci (i.e. tsv file in theta estimation output)

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
tail -n +2 ref_doSAF.thetas.tsv > ref_doSAF.thetas.tsv_nohead
cut -f2 ref_doSAF.thetas.tsv_nohead | awk '{$1 = $1 + 1; print}' | paste ref_doSAF.thetas.tsv_nohead - | awk 'BEGIN {FS="\t"};{ print $1"\t"$2"\t"$8"\t"$4}' | sed 's/ //g' > ref_pi_global.bed

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

##################### for DelBay19 run ##################### 


```
