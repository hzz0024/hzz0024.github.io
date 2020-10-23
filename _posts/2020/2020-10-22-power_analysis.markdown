---
comments: true
title: Power analysis for Fisher and Z methods
date: '2020-10-22 12:00'
tags:
  - DelBay19
  - ouliter
  - WGS
  - Power
  - Fisher
  - Z
categories:
  - WGS data analysis
---

In this post I performed the power analysis for Fisher and weighted Z-test (also called ‘Stouffer’s method, see [Whitlock (2005)](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1420-9101.2005.00917.x)) approaches. The major difference between these two is that during the p-value combination step, the Z method favours symmetric rejection and is less sensitive to a single low p-value, requiring more consistently low p-values to yield a low combined p-value.

The idea comes from the simulation study in Rey et al. 2020 

<img src="https://hzz0024.github.io/images/power/Rey et al 2020.jpeg" alt="img" width="800"/>

Key features in the simulation,

- Delta_p Levene's model

The Levene's moldel was used to calculate the delta_p of selection, as a function of initial allele frequency before selection (P in the figures below) and the strength of selection (s). The fitness of AA is 1+s the fitness of Aa is 1 and the fitness of aa is 1-s. Therefore, allele A determines genotypes' fitness, and selection has an antagonistic and symmetrical effect on the fitness of homozygotes.

```R
### PARAMETERS ###
# Levene is a function to calculate the delta_p with specified selection coefficient. For each s, it will return the delta_p with p ranging from 0-1.
# (-P internal parameter) the vector of allele frequencies before selection
# -s the selection coefficient s

Levene <- function(s)
  {
  P <- seq(0,1,0.01)
  DpLev <- (s*(P-(P^2)))/(1-(2*s*P)+s)
  }
```
With this function we can create a figure to build the relationship among delta_p, s, and P. 

<img src="https://hzz0024.github.io/images/power/delta_p.jpeg" alt="img" width="800"/>

check the dataset of delta_p

```R
p <- rep(seq(0,1,0.01),101)
s <- sort(rep(seq(0,0.5,0.005),101))
delta <- unlist(lapply(seq(0,0.5,0.005),Levene))
data <- data.frame(p,s,delta)

> head(data)
     p s delta
1 0.00 0     0
2 0.01 0     0
3 0.02 0     0
4 0.03 0     0
5 0.04 0     0
6 0.05 0     0

> tail(data)
         p   s       delta
10196 0.95 0.5 0.043181818
10197 0.96 0.5 0.035555556
10198 0.97 0.5 0.027452830
10199 0.98 0.5 0.018846154
10200 0.99 0.5 0.009705882
10201 1.00 0.5 0.000000000
```

- Power detection

First, the initial allele frequency is used to build the pre-selection allele count of size N.

```R
RAC1 <- rbinom(100,(2*N),data[r,1])
```

After then, a random perturbation is applied to the allele frequencies post selection to model the effect of drift in a finite population of size 10000. Then these new post selection frequencies are used to build a post-selection allele count of size K.

```R
SAC1 <- rbinom(100,(2*K),(rbinom(1,20000,(data[r,1]+data[r,3]))/20000))
```

I did this twice to simulate two independent tests (similar to CH-REF and HC-SR). Then, the combined delta_p is calculated for later Fisher exact test. After fisher's exact test, the p-value list from two independent analyses were combined either using Fisher's approach or Z method.

```R
	  if (DP_mean[t]==0)
        {
        PVfish <- c(PVfish, 1)
        }
      if (DP_mean[t]>0)
        {
        CH1 = c(SAC1[t], 2*K-SAC1[t])
        REF1 = c(RAC1[t], (2*N)-RAC1[t])
        M1 = as.table(cbind(CH1, REF1))
        PV1 <- fisher.test(M1, alternative="greater")$p*2 #multiply the p value by 2 in order to account for fact that it should be a 2-tailed test
        CH2 = c(SAC2[t], (2*K)-SAC2[t])
        REF2 = c(RAC2[t], (2*N)-RAC2[t])
        M2 = as.table(cbind(CH2, REF2))
        PV2 <- fisher.test(M2, alternative="greater")$p*2 
        PV1[PV1 > 1] <- 1
        PV2[PV2 > 1] <- 1
        PVfish_cmp = combinePValues(PV1,PV2, method='fisher')
        PVfish <- c(PVfish, PVfish_cmp)
        PVz_cmp = combinePValues(PV1,PV2, method='z') 
        PVz <- c(PVz, PVz_cmp)
        }
      if (DP_mean[t]<0)
        {
        CH1 = c(SAC1[t], 2*K-SAC1[t])
        REF1 = c(RAC1[t], (2*N)-RAC1[t])
        M1 = as.table(cbind(CH1, REF1))
        PV1 <- fisher.test(M1, alternative="less")$p*2 
        CH2 = c(SAC2[t], (2*K)-SAC2[t])
        REF2 = c(RAC2[t], (2*N)-RAC2[t])
        M2 = as.table(cbind(CH2, REF2))
        PV2 <- fisher.test(M2, alternative="less")$p*2 
        PV1[PV1 > 1] <- 1
        PV2[PV2 > 1] <- 1
        PVfish_cmp = combinePValues(PV1,PV2, method='fisher') 
        PVfish <- c(PVfish, PVfish_cmp)
        PVz_cmp = combinePValues(PV1,PV2, method='z') 
        PVz <- c(PVz, PVz_cmp)
        }
      }
```

Lastly, the combined p-value vector was adjusted with FDR. I counted how many p-values with FDR < 0.05 in a 100-times test. Below are the results:

- Fisher's method, N = 100, K = 100

<img src="https://hzz0024.github.io/images/power/100_100_Fisher.jpeg" alt="img" width="800"/>

- Z method, N = 100, K = 100

<img src="https://hzz0024.github.io/images/power/100_100_Z.jpeg" alt="img" width="800"/>

- Z minus Fisher in power

<img src="https://hzz0024.github.io/images/power/100_100_diff.jpeg" alt="img" width="800"/> 

- Fisher's method, N = 50, K = 50

<img src="https://hzz0024.github.io/images/power/50_50_Fisher.jpeg" alt="img" width="800"/>

- Z method, N = 50, K = 50

<img src="https://hzz0024.github.io/images/power/50_50_Z.jpeg" alt="img" width="800"/>

- Z minus Fisher in power

<img src="https://hzz0024.github.io/images/power/50_50_diff.jpeg" alt="img" width="800"/> 

Overall, it looks that the Fisher's combined method has little difference to Z approach, wherears the sample size of N (start population size) and K (size after selection) have a large impact on the power.

We may have less power to detect the outliers under low start P and low selection coefficient situtions.

Some quesions remained in this analysis,

1) Why the delta_p in Levene's equation is alway positive? Where does equation come from? I borrowed this equation from Rey et al 2020 but no clues about its origin yet. Levene's paper can be found [here](https://www.journals.uchicago.edu/doi/pdfplus/10.1086/281792).

2) If negative delta_ps are shown due to drift effect, to what extent they appeared in the dataset?

3) Any other equation can be used for delta_p modeling? See some references [here](http://ww2.tnstate.edu/ganter/BIO416%2006%20Genotype.html).




