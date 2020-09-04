---
comments: true
title: DelBay19 Fishers' exact test update
date: '2020-09-01 12:00'
tags:
  - DelBay19
  - ouliter
  - WGS
  - Fisher's exact
categories:
  - WGS data analysis
---

In this post I made some modifications for Fisher's exact tests. 

1) The SNP lists in the last post [DelBay19 Fishers' exact test](https://github.com/hzz0024/hzz0024.github.io/blob/master/_posts/2020/2020-08-26-Outlier_detection_Fish.markdown) were created seperatly for the challenge and wild groups. I suspect this will reduce the number of common shared outliers in the test results. Here I create a single SNP list by combining all the populations in the Angsd global run, and generated the allele frequency dataset for each population.

2) I am still exploring the best way to obtain the allele frequency values, should be either doMajorMinor 3 plus doMaf 1 or doMajorMinor 5 plus doMaf 2. I made a comparsion in this post.

3) The most important change is the Fisher's test. Following the coral paper publisted by [matzlab](https://matzlab.weebly.com/), I performed Fisher’s exact tests on two "replicate" datasets. One is the challenge vs reference contrast, another one is the wild contrast from Cohansey (COH) vs Shell Rock (SR). P-values from the replicates were then combined using *Fisher’s Combined Probability method*. This approach requires the use of one-tailed tests so that information on the direction of effect is retained, with alternative hypothesis derived from the average of delta_p. More importantly, the p-value is by multipled by 2 to compensate for this post hoc use of the data and inform two-tailed tests. 

Whitlock made a good explaination for this method, which can be found in this paper [Combining probability from independent tests: The weighted Z-method is superior to Fisher’s approach](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1420-9101.2005.00917.x). This method, however, will treats large and small p-values asymmetrically. That is, the combined p-value is sensitive to small p-values compared to large p-values, and is likely to reject the null hypothesis in favour of contradictory one-tailed alternative hypothesis. Anyway, given that the code is already developed for this method, I modified the code to fit into my data and tried to identify the poential outliers from the "replicates" datasets. 

---

### Fisher's exact tests

The R script of Fisher's exact test for each repeated study is located in GitHub/DelBay_project/R_scripts/Fisher_exact/Fish_repeat

The R script of Fisher’s combined probability is located in GitHub/DelBay_project/R_scripts/Fisher_exact/Fish_combine_p

---

### SNPs and Allele Frequency Files (mafs.gz)

The regenerated dataset includes 1,934,038 SNPs after global Angsd calling. This SNP list is used for indivudal population SNP dataset creation.

Two strategies are used for allele frequency data creation,

1) doMajorMinor 3 + doMaf 1 (hereafter as do31)

doMajorMinor 3 ----- the major and minor allele can be predefined for the desired sites
doMaf 1  ----- calculate allele frequencies assuming known major and known minor, provided by SNP list

Same as [Claire's code](https://github.com/clairemerot/angsd_pipeline/blob/master/01_scripts/06_saf_maf_by_pop_maxdepth.sh) and [Nicolas's code](https://github.com/therkildsen-lab/genomic-data-analysis/blob/master/scripts/get_maf_per_pop.sh)

```sh
angsd -P $NB_CPU -doMaf 1 -dosaf 1 -GL 1 -doMajorMinor 3 -anc $ANC -minQ 20 -b $CH $REGIONS -sites ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 -out "/05_saf_maf_by_pop/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30"
```
2) doMajorMinor 5 + doMaf 2 (hereafter as do52)

doMajorMinor 5 ----- use the ancestral allele as the major allele      
doMaf 2 ----- calculate allele frequencies assuming a fixed major allele (ancestral allele) and an unknown minor allele

```sh
angsd -P $NB_CPU -doMaf 2 -dosaf 1 -GL 1 -doMajorMinor 5 -anc $ANC -minQ 20 -b $CH $REGIONS -sites ALL_sites_all_maf0.05_pctind0.7_maxdepth3dv_snplist_4col_cv30 -out "/05_saf_maf_by_pop/"$target"_maf"$MIN_MAF"_pctind"$PERCENT_IND"_cv30_anc"
```

- Results

For output comparsion between do31 and do52, the allele frequency results can be classied into three categories,

Using challenge population (n=50) as an example, 

1) do52$Major = do31$Major, do52$Minor = do31$Minor - No. SNP = 1646633  

2) do52$Major = do31$Minor, do52$Minor = do31$Major - No. SNP = 270369  

   | Condition                                        | Chr/Pos/do52$Major/Minor/Anc/p/do$31Major/Minor/Anc/p|
   |--------------------------------------------------|------------------------------------------------------|
   |do52$Major = do31$Minor, do52$Minor = do31$Major  | NC_035780.1 19095 T G T 0.89514 G T T 0.104855       | 

3) other condition such as

   | Condition                                        | Chr/Pos/do52$Major/Minor/Anc/p/do$31Major/Minor/Anc/p|
   |--------------------------------------------------|------------------------------------------------------|
   |do52$Major = do31$Major, do52$Minor != do31$Minor | NC_035780.1 27711 A T A 0.102502 A C A 0.036589      |
   |do52$Major != do31$Major, do52$Minor != do31$Minor| NC_035780.1 31928 A G A 0.926336 G T A 0.388542      |
   
One benefit of using do52 is that there is no constrain in allele frequency range (from 0-1). However, when I check the major/minor allele consistence across different populations, it appears that the same SNP will sometimes have different minor alleles due to the uncertainy in minor allele determination. 

On the other hand, for do31, because we fix the major and minor alleles using a SNP list, these alleles are always consistent across populations, which is good for delta_p comparsion. One drawback of using this result may be the constrain of frequency ranges (ususally less than 0.8). 

Overall, I deceide to continue using do31 method to produce allele frequency data. 

---

### Fisher’s combined probability test

- CH vs. REF & HC vs. SR manhattan plot

Manhattan plot of 5% FDR for CH vs. REF + HC vs. SR. Red points show markers at 5% FDR according to the Fisher’s combined probability test (n = 5234)

<img src="https://hzz0024.github.io/images/Fish/CH_REF_HC_SR_plot_0.05.jpg" alt="img" width="800"/>

Manhattan plot of 1% FDR for CH vs. REF + HC vs. SR. Red points show markers at 1% FDR according to the Fisher’s combined probability test (n = 87)

<img src="https://hzz0024.github.io/images/Fish/CH_REF_HC_SR_plot_0.01.jpg" alt="img" width="800"/>

- CH vs. REF & HC vs. SR delta_p

Let's look at the delta_p range. Because this is a combined p-value method, which means that the delta_p in one population may be small and but very large in the other population. In such case, I'd like to check the sum of absolute delta_p for each SNP. 

First check the overall delta_p patterns,

```R
delta_p_a = CH$knownEM - REF$knownEM
delta_p_b = SR$knownEM - HC$knownEM
sum( (abs(delta_p_a) + abs(delta_p_b)) > 0.3)
> 21006
length(delta_p_a)
> 1934038
sum( (abs(delta_p_a) + abs(delta_p_a)) > 0.3)/length(delta_p_a)
> 0.01086121
```

This means only 21006 out of 1934038 (1%) SNPs have a sum of delta_p larger than 0.3. What about this ratio in the outliers?

```R
sum( (abs(deltap) + abs(deltap2)) > 0.3)
> 85
length(deltap)
> 87
sum( (abs(deltap) + abs(deltap2)) > 0.3)/length(deltap)
> 0.9770115
```

This means only 85 out of 87 (98%) outliers SNPs have a sum of delta_p larger than 0.3. This ratio is pretty larger than the overall patters. I also randomly pickup up 87 SNPs and examined the propotion of SNPs with sum delta_p > 0.3 and the ratio is 0 (with only a few trial). Overall, the results suggest that Fisher’s combined probability method is useful to find out the outliers.  

HC vs. SR delta_p against p0. Here Black dots are SNPs with absolute delta_p < 0.1, blue dots with  0.1 <= abs(delta_p) < 0.2, gray dots with abs(delta_p) >= 0.2

<img src="https://hzz0024.github.io/images/Fish/HC_SR_fisher1.jpg" alt="img" width="800"/>

CH vs. REF delta_p against p0. This figure is used to trace the movement of SNPs in the HC vs. SR contrast (figure above).

<img src="https://hzz0024.github.io/images/Fish/HC_SR_fisher2.jpg" alt="img" width="800"/>


- CH vs. REF & HC vs. NB manhattan plot

Manhattan plot of 5% FDR for CH vs. REF + HC vs. NB. Red points show markers at 5% FDR according to the Fisher’s combined probability test (n = 5940)

<img src="https://hzz0024.github.io/images/Fish/CH_REF_HC_NB_plot_0.05.jpg" alt="img" width="800"/>

Manhattan plot of 1% FDR for CH vs. REF + HC vs. NB. Red points show markers at 1% FDR according to the Fisher’s combined probability test (n = 171)

<img src="https://hzz0024.github.io/images/Fish/CH_REF_HC_NB_plot_0.01.jpg" alt="img" width="800"/>

- CH vs. REF & HC vs. NB delta_p

HC vs. NB delta_p against p0. Here Black dots are SNPs with absolute delta_p < 0.1, blue dots with  0.1 <= abs(delta_p) < 0.2, gray dots with abs(delta_p) >= 0.2

<img src="https://hzz0024.github.io/images/Fish/HC_NB_fisher1.jpg" alt="img" width="800"/>

CH vs. REF delta_p against p0. This figure is used to trace the movement of SNPs in the HC vs. NB contrast (figure above).

<img src="https://hzz0024.github.io/images/Fish/HC_NB_fisher2.jpg" alt="img" width="800"/>
















