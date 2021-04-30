---
comments: true
title: Polygenic framework
date: '2021-04-23 12:00'
tags:
  - DelBay
  - Wild 
  - WGS
  - C-score
  - polygenic
categories:
  - WGS data analysis
--- 

### What is polygenic selection?

In general, I think the concept of polygenic selection of adaptation is synonyms for genetic redundancy, which infers "situation where more than one combination of genetic variants produces the same phenotype (or overall fitness) over an evolutionary timescale". See Láruson et al. [The Importance of Genetic Redundancy in Evolution](https://www.cell.com/trends/ecology-evolution/fulltext/S0169-5347(20)30116-6). 

### Polygenic framework

The reference paper that helps me develop this polygenic framework is [Ehrlich et al. 2020](https://academic.oup.com/gbe/article/13/2/evaa257/6031913). Using the teleost Fundulus heteroclitus as a model species, they have shown that "populations inhabiting distinct environmental niches exhibit subtle, yet significant changes at many positions in the genome within a single generation, and that these changes lead to genetic differentiation among niches. Subtle changes at multiple genes of small effect may allow organisms to adapt to specific niches over just a single generation."

In this paper they did not apply stringent significance during the early stage of outlier detection. Instead, they tried to combine the raw p-values from each of the outlier identification approaches 1) a permutation analysis (using Fst), 2) a simulation approach, and 3) Barnard’s exact test (similar to Fisher’s exact), and do FDR/Bonferroni correction afterward. In fact, only two SNPs within pond 1 remained significant after multiple test collection. After p-value combination, they only identify three outliers (< FDR 10%) in concordant allele frequency tests (i.e. parallel selection), but found a lot more significant SNPs due to subpopulation polygenic selection.

Below is my polygenic framework,

Populations used for this polygenic framework:

Del19 challenge: REF19-CHR19       
Del20 challenge: REF20-CHR20       
Wild contrasts: SR-HC and NB-HC       

Total SNPs used for this polygenic framework:

SNPs from coverage equalized dataset: 2032113 SNPs

---

1. perform single-generation selection (SGS) on each of the population contrast     
2. <s>perform Barnard's exact test</s>              
3. <s>combine p-value from the tests above for each of the SNP using Z-method, retain SNPs with p-value < 0.05, and SNPs with FDR < 0.05 </s>                         
4. <s>check the combined p-value and Spearman's rank correlation coefficient</s>               
5. calculate the C-score metric to quantify redundancy            
6. perform Cochran-Mantel_haenszel (CMH) test to assess the degree of concordance            
7. perform random forest on potential outliers in each contrast (perhaps SNPs with combined p-value < 0.05?) or the combined populations.         
8. perform gene annotation for potential candidates            

--- 

#### 1 perform single-generation selection on each of the population contrast

Table 1. Single-generation selection (SGS) outliers before (column 2-4) and after FDR correction (column 5-7).

| Populations | SGS (p<0.05) | Positive deltap | Negative deltap | SGS (FDR<0.05) | Positive deltap | Negative deltap |
|-------------|--------------|-----------------|-----------------|----------------|-----------------|-----------------|
| CHR19-REF19 | 257221       | 180780          | 76441           | 2979           | 2467            | 512             |
| CHR20-REF20 | 184137       | 139841          | 44296           | 182            | 164             | 18              |
| HC-SR       | 248994       | 184819          | 64175           | 2443           | 2129            | 314             |
| HC-NB       | 258788       | 192471          | 66317           | 2572           | 2183            | 389             |

Figure 1. P-value distribution before (left) and after FDR correction (right) for CHR19-REF19 SGS test.

<img src="https://hzz0024.github.io/images/polygenic/ps_Del19_challenge.txt.jpg" alt="img" width="800"/>

Table 2. Number of outlier overlaps among the contrasts. Left diagonal: outliers with p-value < 0.05. Right diagonal: outliers with FDR < 0.05 

|             | CHR19-REF19 | CHR20-REF20 | HC-SR | HC-NB |
|-------------|-------------|-------------|-------|-------|
| CHR19-REF19 | -           | 0           | 2     | 8     |
| CHR20-REF20 | 23650       | -           | 0     | 0     |
| HC-SR       | 32052       | 22803       | -     | 104   |
| HC-NB       | 32982       | 23831       | 66381 | -     |

Figure 2. Delta_p distribution (left) and trends (right) for CHR19-REF19 outliers with FDR < 0.05, n=2979).

<img src="https://hzz0024.github.io/images/polygenic/SGS_Del19_FDR_2979.jpg" alt="img" width="800"/>

I also checked the overlaps in outliers with p<0.05, below are the results:

```sh
> a <- intersect(Del19, Del20)
> length(a)
[1] 23650
> b <- intersect(Del19, HC_SR)
> length(b)
[1] 32052
> c <- intersect(Del19, HC_NB)
> length(c)
[1] 32982
> e <- intersect(Del20, HC_SR)
> length(e)
[1] 22803
> d <- intersect(Del20, HC_NB)
> length(d)
[1] 23831
> f <- intersect(HC_NB, HC_SR)
> length(f)
[1] 66381
> g <- intersect(HC_SR, (intersect(Del19, Del20)))
> length(g)
[1] 2966
> h <- intersect(HC_NB, (intersect(Del19, Del20)))
> length(h)
[1] 3077
> i <- intersect(HC_SR, (intersect(HC_NB, (intersect(Del19, Del20)))))
> length(i)
[1] 823
```

Delta_p distribution and trends for 823 common shared outliers is 

<img src="https://hzz0024.github.io/images/polygenic/SGS_shared_823.jpg" alt="img" width="800"/>

#### 2 perform Barnard's exact test on each of the population contrast

| Populations | BE (p<0.05) | Positive deltap | Negative deltap | BE (FDR<0.05) | Positive deltap | Negative deltap |
|-------------|-------------|-----------------|-----------------|---------------|-----------------|-----------------|
| REF19-CHR19 | 216940      | 101763          | 115177          | 266           | 140             | 126             |
| REF20-CHR20 | 134267      | 70468           | 63799           | 0             | 0               | 0               |
| SR-HC       | 199332      | 103834          | 95498           | 58            | 37              | 21              |
| NB-HC       | 207563      | 108923          | 98640           | 105           | 48              | 57              |

Figure 1. P-value distribution before (left) and after FDR correction (right) for CHR19-REF19 Barnard's exact test.

<img src="https://hzz0024.github.io/images/polygenic/barnard_REF19_CHR19.txt.jpg" alt="img" width="800"/>

Table 2. Number of outlier overlaps among the contrasts. Left diagonal: outliers with p-value < 0.05. Right diagonal: outliers with FDR < 0.05 

|             | CHR19-REF19 | CHR20-REF20 | HC-SR | HC-NB |
|-------------|-------------|-------------|-------|-------|
| CHR19-REF19 | -           | 0           | 0     | 0     |
| CHR20-REF20 | 14509       | -           | 0     | 0     |
| HC-SR       | 21706       | 13199       | -     | 3     |
| HC-NB       | 22365       | 13906       | 50494 | -     |

Again, I checked the overlaps in outliers with p<0.05, below are the results:

```sh
> a <- intersect(Del19, Del20)
> length(a)
[1] 14509
> b <- intersect(Del19, HC_SR)
> length(b)
[1] 21706
> c <- intersect(Del19, HC_NB)
> length(c)
[1] 22365
> e <- intersect(Del20, HC_SR)
> length(e)
[1] 13199
> d <- intersect(Del20, HC_NB)
> length(d)
[1] 13906
> f <- intersect(HC_NB, HC_SR)
> length(f)
[1] 50494
> g <- intersect(HC_SR, (intersect(Del19, Del20)))
> length(g)
[1] 1450
> h <- intersect(HC_NB, (intersect(Del19, Del20)))
> length(h)
[1] 1540
> i <- intersect(HC_SR, (intersect(HC_NB, (intersect(Del19, Del20)))))
> length(i)
[1] 396
```

#### 3&4. check the combined p-value and Spearman's rank correlation coefficient     

| Populations | Combined (p<0.05) | Positive deltap | Negative deltap | Combined (FDR<0.05) | Positive deltap | Negative deltap |
|-------------|-------------------|-----------------|-----------------|---------------------|-----------------|-----------------|
| CHR19-REF19 | 256917            | 122859          | 134058          | 5870                | 2980            | 2890            |
| CHR20-REF20 | 162062            | 84897           | 77165           | 441                 | 271             | 170             |
| HC-SR       | 238483            | 123696          | 114787          | 3994                | 2235            | 1759            |
| HC-NB       | 249527            | 130331          | 119196          | 5456                | 3003            | 2453            |

Figure 1. P-value distribution before (left) and after FDR correction (right) for CHR19-REF19 combined results.

<img src="https://hzz0024.github.io/images/polygenic/REF19_CHR19_out_all_z.txt.jpg" alt="img" width="800"/>

Table 2. Number of outlier overlaps among the contrasts. Left diagonal: outliers with p-value < 0.05. Right diagonal: outliers with FDR < 0.05 

|             | CHR19-REF19 | CHR20-REF20 | HC-SR | HC-NB |
|-------------|-------------|-------------|-------|-------|
| CHR19-REF19 | -           | 2           | 24    | 30    |
| CHR20-REF20 | 14509       | -           | 5     | 5     |
| HC-SR       | 21706       | 13199       | -     | 275   |
| HC-NB       | 22365       | 13906       | 50494 | -     |

Shared outliers with p<0.05:

```sh
> a <- intersect(Del19, Del20)
> length(a)
[1] 21136
> b <- intersect(Del19, HC_SR)
> length(b)
[1] 30620
> c <- intersect(Del19, HC_NB)
> length(c)
[1] 32115
> e <- intersect(Del20, HC_SR)
> length(e)
[1] 19553
> d <- intersect(Del20, HC_NB)
> length(d)
[1] 20546
> f <- intersect(HC_NB, HC_SR)
> length(f)
[1] 56854
> g <- intersect(HC_SR, (intersect(Del19, Del20)))
> length(g)
[1] 2585
> h <- intersect(HC_NB, (intersect(Del19, Del20)))
> length(h)
[1] 2766
> i <- intersect(HC_SR, (intersect(HC_NB, (intersect(Del19, Del20)))))
> length(i)
[1] 657
```

It looks that the p-values combination is not very accurate, as the Pearson correlation is low here,

<img src="https://hzz0024.github.io/images/polygenic/Del19_pearson.jpg" alt="img" width="800"/>

Some potentials,

1. Perhaps from the "0" p-value from SGS test. A slight replace might help this.

```sh
dat1 = read.delim("ps_Del19_challenge.txt", header = FALSE, sep='\t')
length(dat1[which(dat1$V6==0),]$V1)
dat1$V6[dat1$V6 == 0 ] <- 1e-5 
```

2. Something wrong with the combination process, need to extract those combination with low correlations and check what is going on.


#### 5. calculate the C-score metric to quantify redundancy

#### 6. perform Cochran-Mantel_haenszel (CMH) test           

### 7. perform random forest        

I followed the tutorial on [Brieuc et al 2015](https://onlinelibrary.wiley.com/doi/abs/10.1111/1755-0998.12773) for random forest.

One important thing is that analyses may be computational intensive - "the method may require extended computational time if purging includes a large number of loci and many trees per forest (e.g., a complete RF analysis that included the generalized purging approach described here for 400 individuals, 9,000 loci, and a continuous trait took approximately 3 weeks to run on a desktop computer with 24 GB of RAM and a processor speed of 3.47 GHz)."

Below is the summary from tutorial and Brieuc et al 2015

---

1. Two types of random forest are available from package randomForest: classification and regression random forest.        
2. The predictive ability of classification trees is measured by the out-of-bag (OOB) error rate.  An error rate is calculated for each tree within a forest. The OOB error rate from the last tree in the forest is usually reported. It takes all previous trees into account and thus represents the error rate after the model stabilizes/converges.        
3. Regression randon forest works on data with a continuous response variable. The proportion variation explained (PVE) will be used to measure the predictive ability of regression trees.   
4. The locus importance values measures how important the SNPs are. It is a key value to estimate the effect size.    
5. Some other key parameters are:    
    
ntree -- number of trees    
mtry -- number of features (i.e. SNPs)    
 
Figures from tutorial (classification forest):

<img src="https://hzz0024.github.io/images/polygenic/RF_1.jpeg" alt="img" width="800"/>

<img src="https://hzz0024.github.io/images/polygenic/RF_2.jpeg" alt="img" width="800"/>

---

#### Things keep in mind:

1. Previous probabilistic random forest (PRF) test using ~ 3000 FDR outliers in Del19 challenge samples has a accuracy of 100%. May need include more SNPs for parameter optimization.        
2. Probabilistic random forest can only support classifization forest but not Regression one, perhaps need to switch to classicial RF if necessary.      
3. If switched to RF, need genotypes without missing data (could use the resequencing data).      

#### Preliminary tests using probabilistic random forest (PRF)

- test with SGS outliers





