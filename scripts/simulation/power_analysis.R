---
  title: "AlleFreqPower"
author: "Katie Lotterhos"
date: "4/28/2019"
output:
  pdf_document: default
  html_document: default
---
  
#setwd("/Users/lotterhos/Documents/GitHub/LifeCyclePower/src")

#{r setup, include=FALSE}
#knitr::opts_chunk$set(echo = TRUE)
# if(!require(devtools)) install.packages("devtools")
# devtools::install_github("kassambara/ggpubr")
library(ggplot2)
library(ggpubr)

## FST function

#This function returns FST for a haploid sample with 2 alleles. Its not designed for pool-seq data but just used here as an approximation.


WC_FST_FiniteSample_Haploids_2AllelesB_MCW<-function(AllCounts){
  #Input a matrix of the counts of each allele (columns) in each population (rows)
  #returns vector instead of list of Fst values, according to Weir
  AllCounts = 1:12
  dim(AllCounts) <- c(6,2)
  AllCounts
  # counts of the rows
  n_pops<-dim(AllCounts)[1]
  r<-n_pops
  # first column
  counts1 <- AllCounts[,1]
  # sum values for each row, i.e. sample size (remember this is for haploid simulation)
  sample_sizes <- rowSums(AllCounts)
  # mean number of sample sizes
  n_ave <- mean(as.numeric(sample_sizes))
  # ?
  n_c = (n_pops*n_ave - sum(sample_sizes^2)/(n_pops*n_ave))/(n_pops-1)  
  # p freq for each sample
  p_freqs = counts1/sample_sizes
  # average p freq across samples
  p_ave = sum(sample_sizes*p_freqs)/(n_ave*n_pops)
  #note: this differs slightly from mean(p_freqs) in R
  # Expected heterozygosity
  He <- 2*p_ave*(1-p_ave)
  
  s2 = sum(sample_sizes*(p_freqs - p_ave)^2)/((n_pops-1)*n_ave)
  #note: this differs slightly from var(p_freqs) in R
  
  T1 <- s2 - 1/(n_ave-1)*(p_ave*(1-p_ave) -(s2*(r-1)/r))
  T2 <- (n_c - 1)*(p_ave*(1-p_ave))/(n_ave-1) + (1 + (r-1)*(n_ave-n_c)/(n_ave-1))*s2/r
  
  FST <- T1/T2 
  
  return(c(He,p_ave, FST, T1, T2))
}

## Allele frequency functions

allelefreqchange.pooled <- function(p0, p1, numind, numreads){
  # p0 true allele frequency at time 0
  # p1 true allele frequency at time 1
  # numind is number of individual sampled at time 0 and again at time 1
  # numreads is the number of reads
  
  p0.samp1 <- rbinom(1, numind, p0)/numind
  p0.samp2 <- rbinom(1, numreads, p0.samp1)/numreads
  
  p1.samp1 <- rbinom(1, numind, p1)/ numind
  p1.samp2 <- rbinom(1, numreads, p1.samp1)/numreads
  
  dp.true<- p1-p0
  dp.1<- p1.samp1-p0.samp1
  dp.2<- p1.samp2-p0.samp2 # this is observed allele frequency
  
  p.bar<- mean(p0,p1)
  p.bar.samp2<- mean(p0.samp2,p1.samp2)
  
  FST.0<-WC_FST_FiniteSample_Haploids_2AllelesB_MCW(
    matrix(c(p0, p1, (1-p0)*numreads, (1-p1)*numreads), ncol=2))
  
  FST.2<-WC_FST_FiniteSample_Haploids_2AllelesB_MCW(
    matrix(c(p0.samp2*numreads, p1.samp2*numreads, (1-p0.samp2)*numreads, (1-p1.samp2)*numreads), ncol=2))
  
  if(p0==p1){return(NA)}else{
    
    return(list(p0=p0,p1=p1, p0.samp2=p0.samp2, p1.samp2=p1.samp2, dp.true=dp.true, dp.1=dp.1, dp.2=dp.2, FST.0=FST.0[3], FST.2=FST.2[3], He.2=FST.2[1], T1=FST.2[4], T2= FST.2[5]))
  }
}


numreads_levels = c(1,10,100)
numind_levels <- c(50,100,1000)

dp <- seq(0.01, 0.2, 0.01)
power<- c()
dat <- data.frame()
length(p0)
k <- 0
for (a in 1: length(numreads_levels)){
  numreads=numreads_levels[a]
  for (b in 1: length(numind_levels)){
    numind = numind_levels[b]
    for (i in 1:length(dp)){
      print(i)
      for (p_0 in seq(0.01, 0.2, 0.01)){
        k <- k+1
        p0 <- rep(p_0, 1000)
        delta_p <- dp[i]
        p1 <- p0 + delta_p
        
        x <- data.frame(t(mapply(allelefreqchange.pooled, p0,p1,numind=numind ,numreads=numreads)))
        x <- sapply( x, as.numeric )
        x <- data.frame(x)
        #head(x)
        null <- data.frame(t(mapply(allelefreqchange.pooled, p0+0.00001,p0,numind=numind ,numreads=numreads)))
        null <- sapply( null, as.numeric )
        null <- data.frame(null)
        #head(null)
        
        cutoff <- quantile(null$dp.2,0.95, type=1)
        power <- sum(x$dp.2 > cutoff)/length(x$dp.2)
        dat0 <- data.frame(numind, numreads, p_0, delta_p, power)
        dat <- rbind(dat, dat0)
      }
    }
  }
}
head(dat)
saveRDS(dat,file = "PowerTable.RDS")

dat_full <- readRDS("PowerTable.RDS")
str(dat_full)
levels(factor(dat_full$numind))
levels(factor(dat_full$numreads))
makeplot <- function(n_ind, n_reads){
  dat1 <- dat_full[dat_full$numind==n_ind & dat_full$numreads==n_reads,]
  e <- ggplot(dat1, aes(delta_p, power))
  #e + geom_point()#, x, y, alpha, color, fill, shape, size, stroke
  #e + geom_point(aes(color=p_0))
  #e + geom_point(aes(color=p_0)) + theme_bw()
  #e + geom_point(aes(color=p_0)) + theme_classic() 
  #e + geom_point(aes(color=p_0)) + theme_classic() + labs(x="Change in allele frequency", y="Power")
  #e + geom_point(aes(color=p_0)) + theme_classic() + labs(x="Change in allele frequency", y="Power") + geom_hline(yintercept=0.8, color="grey")
  #e + geom_point(aes(color=p_0)) + theme_classic() + labs(x="Change in allele frequency", y="Power", color="Intial\nallele\frequency") + geom_hline(yintercept=0.8, color="grey")
  #e + geom_point(aes(color=p_0)) + theme_classic() + labs(x="Change in allele frequency", y="Power", color="Intial\nallele\nfrequency") + geom_hline(yintercept=0.8, color="grey")
  e + geom_point(aes(color=p_0)) + theme_classic() + labs(x="Change in allele frequency", y="Power", color="Intial MAF") + geom_hline(yintercept=0.8, color="grey") + ylim(0,1) + ggtitle(paste0("# Ind.=",n_ind, ", # Reads=", n_reads))
  #+ geom_vline(xintercept=0.1, color="grey")
}

numind_levels
numreads_levels
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
i50_r1 <-   makeplot(50, 1)
i100_r1 <-   makeplot(100, 1)

i50_r10 <-  makeplot(50, 10)
i100_r10 <-  makeplot(100, 10)

i50_r100 <-  makeplot(50, 100)
i100_r100 <-   makeplot(100, 100)
ggarrange(i50_r1, i100_r1,
          i50_r10, i100_r10,
          i50_r100, i100_r100,
          labels = LETTERS[1:6],
          ncol = 2, nrow = 3)

