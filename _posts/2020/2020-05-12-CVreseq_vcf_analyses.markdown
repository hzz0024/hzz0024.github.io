---
comments: true
title: CVreseq vcf data analyese
date: '2020-05-13 11:00'
tags:
  - CVreseq
  - vcftools
  - plink
  - fst
  - dxy
  - pi
  - ADMIXTURE
  - Ho
  - He
  - Ne
categories:
  - CVreseq data analysis
---

This post will summarize some basic stats for CVreseq data

Some scripts are adopted from [here](https://github.com/grovesdixon/caveRAD/blob/master/cave_RAD_processing_walkthrough.txt)

Also I found a good paper to follow the ways of measuring genetic diversity, inbreeding, and effective size 
https://gsejournal.biomedcentral.com/articles/10.1186/s12711-019-0468-4

--- 

vcftools documentation: https://vcftools.github.io/examples.html

### Basic stats with vcftools

```sh
#MAKE POPULATION FILES: 

CS   
DEBY   
NEH   
OBOYS2   
SL   

# get pairwise fst between populations

vcftools --vcf LA.recode.vcf --maf 0.05 --recode --recode-INFO-all --out LA_maf
>After filtering, kept 284587 out of a possible 334011 Sites

vcftools --vcf DB_1.recode.vcf --maf 0.05 --recode --recode-INFO-all --out DB_1_maf
>After filtering, kept 282067 out of a possible 334011 Sites

vcftools --vcf DB_2.recode.vcf --maf 0.05 --recode --recode-INFO-all --out DB_2_maf
>After filtering, kept 293946 out of a possible 334011 Sites

vcftools --vcf LA_maf.recode.vcf --weir-fst-pop SL --weir-fst-pop OBOYS2 --out SL_OBOYS2
Weir and Cockerham mean Fst estimate: 0.024711
Weir and Cockerham weighted Fst estimate: 0.036727
After filtering, kept 284587 out of a possible 284587 Sites

vcftools --vcf DB_1_maf.recode.vcf --weir-fst-pop CS --weir-fst-pop NEH --out CS_NEH
Weir and Cockerham mean Fst estimate: 0.048347
Weir and Cockerham weighted Fst estimate: 0.063669
After filtering, kept 282067 out of a possible 282067 Sites

vcftools --vcf DB_2_maf.recode.vcf --weir-fst-pop CS --weir-fst-pop DEBY --out CS_DEBY
Weir and Cockerham mean Fst estimate: 0.02112
Weir and Cockerham weighted Fst estimate: 0.033042
After filtering, kept 293946 out of a possible 293946 Sites

#assemble results into single tsv
echo -e "CHROM\tPOS\tWEIR_AND_COCKERHAM_FST\tgroup" > all_fst_results.tsv
for file in *weir.fst
do echo "${file}..."
 tail -n +2 $file | awk -v f="${file/.weir.fst/}" '{print $0"\t"f}' >> all_fst_results.tsv
done
-----------------------

##get pairwise dxy between populations

./vcf_dxy.R --vcf LA_maf.recode.vcf --popX SL --popY OBOYS2 --out SL_OBOYS2.dxy
./vcf_dxy.R --vcf DB_1_maf.recode.vcf --popX CS --popY NEH --out CS_NEH.dxy
./vcf_dxy.R --vcf DB_2_maf.recode.vcf --popX CS --popY DEBY --out CS_DEBY.dxy

#assemble results into single tsv
echo -e "CHROM\tPOS\tdxy\tgroup" > all_dxy_results.tsv
for file in *.dxy
do echo "${file}..."
 tail -n +2 $file | awk -v f="${file/.dxy/}" '{print $0"\t"f}' >> all_dxy_results.tsv
done
```
### using plot_fst_dxy.R for fst and dxy plotting

# Groups:   group, spp [3]     
  group     spp   stat      mn           
1 CS-DEBY   n     fst   0.0211      
2 SL-OBOYS2 n     fst   0.0247      
3 CS-NEH    n     fst   0.0483      

# Groups:   group, spp [3]      
  group     spp   stat      mn            
1 CS-DEBY   n     dxy   0.383     
2 CS-NEH    n     dxy   0.384     
3 SL-OBOYS2 n     dxy   0.386      
 
```