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

In this post I collected some Angsd parameters from the related literatures. Some of the paper are inspiring. The main purpose of doing it to figure out why allele frequency values in the mafs file do not go beyond 0.8. Is it due to some technical problems? Or just a refection to the data itself?

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

> Using the software ANGSD, we identified a set of SNPs with minor allele frequency ≥ 5% in the combined sample based on read counts (-doMajorMinor 2 -doCounts 1). We excluded unplaced scaffolds (<5Mb total) and applied standard quality filters for SNP calling (base quality ≥ 20, mapping quality ≥ 30, total read depth ≤ 5500 (∼2x mean), and coverage across individuals ≥ 50%). 

> we used ANGSD to call genotypes (-doPost 1) using a minor allele frequency prior (-doMaf 1) and the SAMtools genotype likelihood (-GL 1), after quality filtering (map quality ≥ 30, reads matching major/minor allele ≥ 60%, and read depth ≥ 6x)"

- Delmore et al. 2020 The evolutionary history and genomics of European blackcap migration [link](https://elifesciences.org/articles/54462)

This elife paper used ΔPBS analysis to identify regions within each population that exhibit differences in allele frequencies. 

Not very informative in terms of Angsd -doMAF running though.

- Hooper et al. 2020 Runs of homozygosity in killer whale genomes provide a global record of demographic histories [link](https://www.biorxiv.org/content/10.1101/2020.04.08.031344v1.full.pdf)

Only some useful information about ROH,

> Runs of homozygous genotypes (ROH) were identified using the window-based approach implemented in PLINK v1.07 [16] from an input file of genotype likelihoods generated from autosomal scaffolds >10Mb by ANGSD [57] with the following filtering settings: removing reads of poor mapping quality (MAPQ < 30), removing sites with low base quality scores (q < 20), calling only SNPs inferred with a likelihood ratio test (LRT) of P < 0.000001, discarding reads that did not map uniquely, adjusting q-scores around indels, adjusting minimum quality score to 50 for excessive mismatches, and discarding bad reads (flag >=256). Inferred SNPs were lightly pruned based on linkage disequilibrium (LD) r2 >0.9 using PLINK, which has been found to improve the accuracy of detecting autozygous ROHs [58]. We estimated ROH from pruned and unpruned data and found minimal qualitative difference with our data. Sliding window size was set to 300 kb, with a minimum of 50 SNPs at a minimum density of 1 SNP per 50 kb required to call a ROH. To account for genotyping errors, we allowed up to 5 heterozygote sites per 300 kb window within called ROHs, as per ref. [59]. A length of 1,000 kb between two SNPs was required in order them to be considered in two different ROHs"

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

> ANGSD v.0.929......Bases with a q-score below 20 (-minQ 20) and reads with mapping quality below 25 (-minMapQ 25) were removed and variants were only retained if they had a P value smaller than 1 × 10–6 (-SNP_pval 1 × 10–6 flag in ANGSD). We retained sites with a minimum read depth of two (-minIndDepth 2) in at least 80% of the sampled individuals (-minInd 133). The sex chromosome (chromosome XIX; refs. 63,64) was excluded from downstream analyses due to sex-specific genomic heterogeneity65,66. The raw output of genotype likelihoods from all 166 individuals comprised 2,511,922 genome-wide loci."

- Wilder et al 2020 Footprints of local adaptation span hundreds of linked genes in the Atlantic silverside genome. [link](https://onlinelibrary.wiley.com/doi/full/10.1002/evl3.189)

> ANGSD version 0.912......Major and minor alleles, minor allele frequencies (MAF), and FST were estimated in ANGSD from genotype likelihoods for each SNP with data for at least 10 individuals per population (Korneliussen et al. 2014)." 

---

### Small test focusing on the potential outliers

A text file is created for -rf usage, which spans ~2M bp in chromosome NC_035784.1

```sh
cat rf_test.txt
> NC_035784.1:12000000-14000000
# create a snp list within this region
cat ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 | grep -A 18000 "NC_035784.1 12004160" > test_snp.list
cat test_snp.list | wc -l
> 18001
# start from NC_035784.1 12004160 to NC_035784.1 14075213
# bash command to estimate the mean, min and max allele frequency in the maf outputs
zcat test1.mafs.gz | awk '{print $6}' | tail -n +2| awk '{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1<min) {min=$1}; total+=$1; count+=1} END {print total/count, max, min}'
```

- test 1 MAF with −doMajorMinor 5, −doMaf 2 (similar to Castruita et al. 2020)

−doMajorMinor 5 ----- use the ancestral allele as the major allele      
−doMaf 2 ----- calculate allele frequencies assuming a fixed major allele and an unknown minor allele

```sh
angsd -P $NB_CPU -doMaf 2 -dosaf 1 -GL 1 -domajorminor 5 \
-anc $ANC -minQ 20 -b $CH \
-sites test_snp.list $REGIONS -out test/test1
```

Number of SNPs: 17419

Allele frequency range:

|   Mean   |    Max    |    Min    |
|----------|-----------|-----------|
|  0.263902| 0.999999  | 0.000000  |

Detailed maf file (frist five SNPs)

|   chromo   | position  |  major    |   minor   |   anc   |  unknownEM   |   nInd   |
|------------|-----------|-----------|-----------|---------|--------------|----------|
|NC_035784.1 | 12004160  | T         |      C    |    T    |   0.125158   |    39    |
|NC_035784.1 | 12021111  | A         |      T    |    A    |   0.211960   |    38    |
|NC_035784.1 | 12021201  | C         |      A    |    C    |   0.090226   |    42    |
|NC_035784.1 | 12021211  | T         |      A    |    T    |   0.443017   |    43    |
|NC_035784.1 | 12021217  | A         |      C    |    A    |   0.339384   |    43    |

- test 2 MAF with -doMajorMinor 3 and -doMaf 1. Same as [Claire's code](https://github.com/clairemerot/angsd_pipeline/blob/master/01_scripts/06_saf_maf_by_pop_maxdepth.sh) and [Nicolas's code](https://github.com/therkildsen-lab/genomic-data-analysis/blob/master/scripts/get_maf_per_pop.sh)

```sh
angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -domajorminor 3 \
-anc $ANC -minQ 20 -b $CH \
-sites test_snp.list $REGIONS -out test/test2
```

Number of SNPs: 17419

Allele frequency range:

|   Mean   |    Max    |    Min    |
|----------|-----------|-----------|
| 0.198844 | 0.707927  | 0.000000  |

Detailed maf file (frist five SNPs)

|   chromo   | position  |  major    |   minor   |   anc   |    knownEM   |   nInd   |
|------------|-----------|-----------|-----------|---------|--------------|----------|
|NC_035784.1 | 12004160  | T         |      C    |    T    |   0.125153   |    39    |
|NC_035784.1 | 12021111  | A         |      T    |    A    |   0.211954   |    38    |
|NC_035784.1 | 12021201  | C         |      A    |    C    |   0.090220   |    42    |
|NC_035784.1 | 12021211  | T         |      A    |    T    |   0.443018   |    43    |
|NC_035784.1 | 12021217  | A         |      C    |    A    |   0.339376   |    43    |


- Do SNP shows allele frequency larger then 0.8, or even fixed? 

Yes, SNPs with allele frequencies > 0.8 are observed in the test 1 resutls (with −doMajorMinor 5, −doMaf 2). I decide to take a look at those SNPs with high allele frequency

```sh
zcat test1.mafs.gz | grep "0.999998"
NC_035784.1	12926755	T	A	T	0.999998	39
NC_035784.1	13399467	T	C	T	0.999998	36
```

How about same SNPs from test 2?

```sh
# test 2
zcat test2.mafs.gz | grep "12926755"
NC_035784.1	12926755	A	T	T	0.000001	39
zcat test2.mafs.gz | grep "13399467"
NC_035784.1	13399467	C	A	T	0.182822	43
```

It seems for SNP NC_035784.1_13399467, Angsd will calculate allele frequency (p) for different alternative alleles. In test1 it estimate the p of C, while in test2 the A allele is used for p calculation. The reason why this happend is due the setting of −doMajorMinor 5 in test1 (i.e. force the major allelel according to your ancestral states) and −doMajorMinor 3 (major and minor alleles are predefined for the desired sites). Besides, result of test2 prefers to report the minor allele frequencies (see NC_035784.1_12926755)

- Compare between populations and different parameters

```sh
#column headers
chromo	position  major	minor anc  frequency  nInd
#SNP ID
snp NC_035784.1_12926755
#parameters and detailed information for specific SNP
-domajorminor 5 -domaf 2 
pop_CH  NC_035784.1	12926755	T	A	T	0.999998	39
pop_REF NC_035784.1	12926755	T	A	T	0.952048	39
-domajorminor 5 -domaf 1 
pop_CH  NC_035784.1	12926755	T	A	T	0.999998	39
pop_REF NC_035784.1	12926755	T	A	T	0.952048	39
-domajorminor 3 -domaf 1 
pop_CH  NC_035784.1	12926755	A	T	T	0.000001	39
pop_REF NC_035784.1	12926755	A	T	T	0.047948	39
-domajorminor 3 -domaf 2 
pop_CH  NC_035784.1	12926755	A	T	T	0.000001	39
pop_REF NC_035784.1	12926755	A	T	T	0.047804	39
-domajorminor 1 -domaf 1 
pop_CH  NC_035784.1	12926755	A	C	T	0.000001	39
pop_REF NC_035784.1	12926755	A	T	T	0.047948	39
-domajorminor 1 -domaf 2 
pop_CH  NC_035784.1	12926755	A	C	T	0.000001	39
pop_REF NC_035784.1	12926755	A	T	T	0.047804	39

#column headers
chromo	position  major	minor anc  frequency  nInd
#SNP ID
snp NC_035784.1_13399467
#parameters and detailed information for specific SNP
-domajorminor 5 -domaf 2 
pop_CH  NC_035784.1	13399467	T	C	T	0.999997	36
pop_REF NC_035784.1	13399467	T	C	T	1.000000	42
-domajorminor 5 -domaf 1 
pop_CH  NC_035784.1	13399467	T	C	T	0.999997	36
pop_REF NC_035784.1	13399467	T	C	T	1.000000	42
-domajorminor 3 -domaf 1 
pop_CH  NC_035784.1	13399467	C	A	T	0.182822	43
pop_REF NC_035784.1	13399467	C	A	T	0.108496	45
-domajorminor 3 -domaf 2 
pop_CH  NC_035784.1	13399467	C	A	T	0.182827	43
pop_REF NC_035784.1	13399467	C	A	T	0.108500	45
-domajorminor 1 -domaf 1 
pop_CH  NC_035784.1	13399467	C	A	T	0.182822	43
pop_REF NC_035784.1	13399467	C	A	T	0.108496	45
-domajorminor 1 -domaf 2 
pop_CH  NC_035784.1	13399467	C	A	T	0.182827	43
pop_REF NC_035784.1	13399467	C	A	T	0.108500	45
```

1) From the comparsion above it seems there is no big difference between the -doMaf 1 and 2.    
2) -domajorminor 3 and -domajorminor 5 may report for different alleles, and result from -domajorminor 3 always report the lower p value of fixed allele.     
3) -domajorminor 1 and -domajorminor 3 may infer different minor alleles with low p, which make sense.      
4) SNP like NC_035784.1_13399467 will show three allele in the maf result, which are major, minor and anc alleles. I tried to use -skipTriallelic 1 to remove this loci in the test data but failed, indicing this is not a triallelic locus.     
5) the major and minor allele status are consistant across populations (e.g. CH and REF here)    

Overall, it sounds the original methods for mafs data creation are correct. It would be interesting to test the results from -domajorminor 5 -domaf 2 settings.

I also asked questions [here](https://github.com/ANGSD/angsd/issues/342). I look farward to hearing some thoughts about the best way to approach allele frequency values.
