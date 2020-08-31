---
comments: true
title: MAF evaluation
date: '2020-08-30 12:00'
tags:
  - DelBay19
  - MAF
  - WGS
categories:
  - WGS data analysis
---

In this post I collected some Angsd parameters from the related literatures. Some of the paper are inspiring for data analyses using Angsd outputs. The main purpose of doing it to figure out why allele frequency values in the mafs file do not go beyond 0.8. Is it due to some technical problems? Or just a refection to the data itself?

--- 

### Angsd parameters in the literatures

- Stahlke et al. 2020. Historical and contemporary signatures of selection in response to transmissible cancer in the Tasmanian Devil (Sarcophilus harrisii) [link](https://www.biorxiv.org/content/biorxiv/early/2020/08/07/2020.08.07.241885/DC1/embed/media-1.pdf?download=true)

<img src="https://hzz0024.github.io/images/MAF/TableS2.jpg" alt="img" width="800"/>

- Castruita et al. 2020. Analyses of key genes involved in Arctic adaptation in polar bears suggest selection on both standing variation and de novo mutations played an important role [link](https://link.springer.com/article/10.1186/s12864-020-06940-0)

> We used a genotype likelihood approach to construct the PCAs: input genotype likelihood files were constructed using ANGSD v0.929, with the SAMtools genotype likelihood algorithm (−GL 1), and specifying the following parameters: remove reads that have multiple mapping best hits (−unique_only), remove reads with a flag above 255/secondary hits (−remove_bads), include only read pairs with both mates mapping correctly (−only_proper_pairs), adjust mapQ for reads with excessive mismatches (−C 50), adjust quality scores around indels (−baq 1), a minimum mapping quality of 20 (−minMapQ 20), a minimum base quality of 20 (−minQ 20), discard sites where there is no data on at least 95% of the individuals (−minInd), skip tri-allelic sites (−skipTriallelic), and remove SNP sites with a p-value larger than 1e− 6 (−SNP_pval 1e-6). The ANGSD output beagle file was run through PCAngsd v0.95"

> Genotypes were called using ANGSD, specifying the same parameters as the PCA analyses with the additional parameters: write major and minor alleles and the genotype directly (−doGeno 5), estimate the posterior genotype probability based on the allele frequency as a prior (−doPost 1), use the reference allele as the major allele (−doMajorMinor 4), output as beagle likelihood file (−doGlf 2), and calculate allele frequencies assuming a fixed major allele and an unknown minor allele (−doMaf 2)"

- Fuller et al. 2020. Population genetics of the coral Acropora millepora: Toward genomic prediction of bleaching [link](https://science.sciencemag.org/content/sci/suppl/2020/07/15/369.6501.eaba4674.DC1/aba4674-Fuller-SM.pdf)

This science paper shows some interesting initial data processing/trimming parts with low-coverage sequencing data (~ 1.43 X). Maybe useful for my own analyses.

Not very informative in terms of Angsd -doMAF running though.
 
- Calfee et al 2020 (Graham Coop lab). Selection and hybridization shaped the Africanized honey bee invasion of the Americas [link](https://www.biorxiv.org/content/10.1101/2020.03.17.994632v2.full.pdf)

"Using the software ANGSD, we identified a set of SNPs with minor allele frequency ≥ 5% in the combined sample based on read counts (-doMajorMinor 2 -doCounts 1). We excluded unplaced scaffolds (<5Mb total) and applied standard quality filters for SNP calling (base quality ≥ 20, mapping quality ≥ 30, total read depth ≤ 5500 (∼2x mean), and coverage across individuals ≥ 50%). 

we used ANGSD to call genotypes (-doPost 1) using a minor allele frequency prior (-doMaf 1) and the SAMtools genotype likelihood (-GL 1), after quality filtering (map quality ≥ 30, reads matching major/minor allele ≥ 60%, and read depth ≥ 6x)"

- Delmore et al. 2020 The evolutionary history and genomics of European blackcap migration [link](https://elifesciences.org/articles/54462)

This elife paper used ΔPBS analysis to identify regions within each population that exhibit differences in allele frequencies. 

Not very informative in terms of Angsd -doMAF running though.

- Hooper et al. 2020 Runs of homozygosity in killer whale genomes provide a global record of demographic histories [link](https://www.biorxiv.org/content/10.1101/2020.04.08.031344v1.full.pdf)

Only some useful information about ROH,

"Runs of homozygous genotypes (ROH) were identified using the window-based approach implemented in PLINK v1.07 [16] from an input file of genotype likelihoods generated from autosomal scaffolds >10Mb by ANGSD [57] with the following filtering settings: removing reads of poor mapping quality (MAPQ < 30), removing sites with low base quality scores (q < 20), calling only SNPs inferred with a likelihood ratio test (LRT) of P < 0.000001, discarding reads that did not map uniquely, adjusting q-scores around indels, adjusting minimum quality score to 50 for excessive mismatches, and discarding bad reads (flag >=256). Inferred SNPs were lightly pruned based on linkage disequilibrium (LD) r2 >0.9 using PLINK, which has been found to improve the accuracy of detecting autozygous ROHs [58]. We estimated ROH from pruned and unpruned data and found minimal qualitative difference with our data. Sliding window size was set to 300 kb, with a minimum of 50 SNPs at a minimum density of 1 SNP per 50 kb required to call a ROH. To account for genotyping errors, we allowed up to 5 heterozygote sites per 300 kb window within called ROHs, as per ref. [59]. A length of 1,000 kb between two SNPs was required in order them to be considered in two different ROHs"

- Fang et al. 2020 On the causes of geographically heterogeneous parallel evolution in sticklebacks [link](https://www.nature.com/articles/s41559-020-1222-6)

All code are available in [DRYAD repository](https://datadryad.org/stash/dataset/doi:10.5061/dryad.b2rbnzsb1)--stored in cornell/test/

An example of using -sites option in this paper

```sh
angsd -bam $BAMLIST_3sp/pureEP_all_26IDs_bamlist \
-anc $REF_3sp -sites $SITES_3sp/EP_WP_ATL_overlap.txt -rf $SITES_3sp/chrs_3sp.txt \
-doMaf 1 -doMajorMinor 1 -GL 1 -P 15 \
-SNP_pval 1e-6 \
-out $SAF_3sp/pureEP_all_26IDs_bamlist
```

My angsd running code

```sh
angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doMajorMinor 3 -anc $ANC \
-remove_bads 1 -minMapQ 30 -minQ 20 -b $CH \
-sites ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 \
-out "~/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"
```

"ANGSD v.0.929......Bases with a q-score below 20 (-minQ 20) and reads with mapping quality below 25 (-minMapQ 25) were removed and variants were only retained if they had a P value smaller than 1 × 10–6 (-SNP_pval 1 × 10–6 flag in ANGSD). We retained sites with a minimum read depth of two (-minIndDepth 2) in at least 80% of the sampled individuals (-minInd 133). The sex chromosome (chromosome XIX; refs. 63,64) was excluded from downstream analyses due to sex-specific genomic heterogeneity65,66. The raw output of genotype likelihoods from all 166 individuals comprised 2,511,922 genome-wide loci."

- Wilder et al 2020 Footprints of local adaptation span hundreds of linked genes in the Atlantic silverside genome. [link](https://onlinelibrary.wiley.com/doi/full/10.1002/evl3.189)

"ANGSD version 0.912......Major and minor alleles, minor allele frequencies (MAF), and FST were estimated in ANGSD from genotype likelihoods for each SNP with data for at least 10 individuals per population (Korneliussen et al. 2014)." 
