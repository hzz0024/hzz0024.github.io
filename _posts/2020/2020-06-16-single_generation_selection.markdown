---
comments: true
title: Single generation selection (SGS) test
date: '2020-06-17 12:00'
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

### Test for single-generation selection

In [Within-Generation Polygenic Selection Shapes Fitness-Related Traits across Environments in Juvenile Sea Bream](https://www.mdpi.com/2073-4425/11/4/398/htm#app1-genes-11-00398), the author estimated the expected distribution of allele frequency changes (deltap) due to finite sample sizes in the absence of selection to generate null expectations.

### Conceptional steps:

1) a panmictic common gene pool of finite size (i.e. the real population)

2) two samples of size N1 and N2 are drawn within the same generation

3) compare the observed allele frequency difference - ΔP=abs(P1−P2) with the distribution of ΔP expected from random sampling effects (i.e. due to finite sample size effects)

> Question: here the ΔP is estimated from absolute freqeuency changes. Shouldn't it to be the actual differences?

### Part I Build the Bayes’ theorem equation

<img src="https://hzz0024.github.io/images/SGS/SGS_1.jpg" alt="img" width="800"/>

Based on Bayes’ theorem equation above, the posterior probability distribution of Z in the common gene pool given k observed allele counts in N1 diploid individuals is could be estimated from three parts,

1. probability mass function of the binomial distribution, 

```R
# I used upper_left to represent this part of the Bayes’ equation
upper_left = factorial(n)/(factorial(n - k)*factorial(k)) * p^k * (1-p)^(n-k)
```

2. prior probability distribution of Z, given a genetic diversity parameter 𝜃=4𝑁e𝜇

> this is derived from equation (50) in Tajima (1989), Gn(i) = M (1/i + 1/(n-1)). Here i should be the seq(1:(n-1)).

> 𝜃 could be estimated by equation (5) in Tajima (1989), M = s/a1, where s is the number of segregating sites (snp), while a1 is the harmonic number.

```R
# I used tajima_points to represent this part of the Bayes’ equation
  tajima_points = c()
  for(i in seq(1:(n-1))){
    upper_right =  M * (1/i + 1/(n-i))
    tajima_points = c(tajima_points, upper_right)
  }
```

3. The probability Pr(X=k) of observing k allele counts was given a uniform value of 1(2N−1)

```R
1/(n-1)
```

4. because the in Figure S3 and S4 the y-axis is the probability values, each results with different p should be scaled probability. This could be done by 

```R
tajima_points = tajima_points/sum(tajima_points)
points=points/sum(points)
```

Overall, the code for Figure S3 and S4 are shown below,

```R
# these values are identical to data used for Figure S3 and S4, note that N = 50 in Figure S4 equals n=100 DNA sequences
n = 100
k = 25
M = 0.0037
  
# create the function
draw_distribution <- function(n,k,M){
  tajima_points = c()
  for(i in seq(1:(n-1))){
    upper_right =  M * (1/i + 1/(n-i))
    tajima_points = c(tajima_points, upper_right)
  }

  tajima_points = tajima_points/sum(tajima_points) 
  points = c()
  for(i in seq(1:(n-1))){
    p = i/n
    upper_left = factorial(n)/(factorial(n - k)*factorial(k)) * p^k * (1-p)^(n-k)
    upper_right =  M * (1/i + 1/(n-i))/sum(tajima_points)
    res = upper_left * upper_right / (1/(n-1))
    points = c(points, res)
  }
  points=points/sum(points)

  return(list("tajima" = tajima_points, "points" = points))
}


res1 = draw_distribution(n, k, 0.0037)
```

5. Draw the figures

- Figure S3

```R
plot( seq(1:(n-1))/n, res1$points, xlab="p", ylab="probability")
```

<img src="https://hzz0024.github.io/images/SGS/S3.jepg" alt="img" width="800"/>


- Figure S4

```R
plot( seq(1:(n-1))/n, res1$tajima, xlab="p", ylab="probability")
```

<img src="https://hzz0024.github.io/images/SGS/S4.jepg" alt="img" width="800"/>

In this part the remained question is 

- where should I obtain the theta (M) values, from the angsd result or from the equation (5) in Tajima (1989)

### Part II Assessment of the deltap and power of the Bayesian test

I am still confused about this part, the key questions are:

1) How can I compare the ΔP from N1 and N2 in the null model, this sentence is really confusing, 

> The posterior probability distribution of 𝑍 in the common gene pool given X=𝑘 observed allele counts in N1 diploid individuals was finally used to obtain the null distribution of ΔP between the two samples of size N1 and N2. For that, we assumed that the two samples are drawn from the same common gene pool, and used the posterior probability distribution of Z conditioned the first sample of size N1 to predict the null distribution of allele counts in the second sample of size N2. Finally, the distribution of allele frequency differences between the two samples was computed and compared to the observed value of ΔP to estimate a P-value.

2) How should I estimte the p-value, a simple t-test?

3) What is the difference between SGS and Fisher's exact test (in terms of coding)


### Part III Incorprate the uncertainty of frequency (p) in the model

No investigation yet. Perhaps some hints from paper below

[Experimental evidence for ecological selection on genome variation in the wild](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12238)



   

