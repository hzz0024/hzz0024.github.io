#BiocManager::install("LEA")
library(LEA)
genotype = lfmm2geno("genotype1.lfmm")
obj.snmf = snmf(genotype, K = 1:6, entropy = T, ploidy = 2, project = "new")
plot(obj.snmf)
barplot(t(Q(obj.snmf, K = 2)), col = 2:3, ylab = "Ancestry coefficients", xlab = "Sampled individuals")
#----------------------------------------------------------------
# define the function fst function as follows.
fst = function(project,run = 1, K, ploidy = 2){
  library(LEA)
  l = dim(G(project, K = K, run = run))[1]
  q = apply(Q(project, K = K, run = run), MARGIN = 2, mean)
  if (ploidy == 2) {
    G1.t = G(project, K = K, run = run)[seq(2,l,by = 3),]
    G2.t = G(project, K = K, run = run)[seq(3,l,by = 3),]
    freq = G1.t/2 + G2.t
    }
  else {
    freq = G(project, K = K, run = run)[seq(2,l,by = 2),]}
  H.s = apply(freq*(1-freq), MARGIN = 1, FUN = function(x) sum(q*x))
  P.t = apply(freq, MARGIN = 1, FUN = function(x) sum(q*x))
  H.t = P.t*(1-P.t)
  return(1-H.s/H.t)
}
# Then we compute the FST statistics as follows, here fst.values is just fst values for 500 snps
fst.values = fst(obj.snmf, K=2)
#convert fst values into z-scores (absolute values), here n is the number of individuals (200) 
n = dim(Q(obj.snmf, K = 2))[1]
# here change the negative fst value into 1e-6
fst.values[fst.values<0] = 0.000001
K = 2
z.scores = sqrt(fst.values*(n-K)/(1-fst.values))
#Compute the GIF
K=2
lambda = median(z.scores^2)/qchisq(1/2, df = K-1)
lambda
# compute adjusted p-values from the combined z-scores
adj.p.values = pchisq(z.scores^2/lambda, df = K-1, lower = FALSE)
#histogram of p-values
hist(adj.p.values, col = "red")

# compute adjusted p-values from the combined z-scores
adj.p.values = pchisq(z.scores^2/1.6, df = 1, lower = FALSE)
#histogram of p-values
hist(adj.p.values, col = "green")
## FDR control: Benjamini-Hochberg at level q
L = 500
q = 0.1
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates = order(adj.p.values)[w]
#Now, our list of candidate loci is recorded in the R object “candidates”.
plot(-log10(adj.p.values), main="Manhattan plot", xlab = "Locus", cex = .7, col = "grey")
points(candidates, -log10(adj.p.values)[candidates], pch = 19, cex = .7, col = "red")

#-----------------------------------------------------------------------
# start to measure the p-value in challenge dataset
library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)

# first make qq-plot for challenge data
DT <- fread("ch_ref_1kb_fold.fst")
print(DT)
fst.values = DT$angsd_Fst
# here change the negative fst value into 0
fst.values[fst.values<0] = 0
# make qqplot for fst
qqplot(rexp(length(fst.values)),
       fst.values, xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)

#---------------test challenge data in window---------------#
# load the data set
DT <- fread("ch_ref_1kb_fold.fst")
print(DT)
# here fst.values is just fst values for snps
fst.values = DT$angsd_Fst
n = 98
# here change the negative fst value into 1e-6
fst.values[fst.values<0] = 1e-10
#fst.values = fst.values[fst.values!=0]
fst.values
K = 2
z.scores = sqrt(fst.values*(n-K)/(1-fst.values))
z.scores
#Compute the GIF
K=2
lambda = median(z.scores^2)/qchisq(1/2, df = K-1)
lambda
# compute adjusted p-values from the combined z-scores
adj.p.values = pchisq(z.scores^2/lambda, df = K-1, lower = FALSE)
#histogram of p-values
hist(adj.p.values, col = "red")
# make qqplot for pvalue
qqplot(rexp(length(adj.p.values), rate = log(10)),
       -log10(adj.p.values), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)

plot(-log10(adj.p.values), main="Manhattan plot", xlab = "Locus", cex = .7, col = "grey")
L = length(adj.p.values)
L 
q = 1e-7
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates = order(adj.p.values)[w]
candidates
points(candidates, -log10(adj.p.values)[candidates], pch = 19, cex = .7, col = "red")

# --------------------------------------------------
# load the data set
DT <- fread("ch_ref_1kb_fold.fst")
print(DT)
# here fst.values is just fst values for snps
fst.values = DT$angsd_Fst
n = 98
# here change the negative fst value into 1e-6
#fst.values[fst.values<0] = 1e-10
fst.values[fst.values<0] = 0
#fst.values = fst.values[fst.values!=0]
fst.values
K = 2
z.scores = sqrt(fst.values*(n-K)/(1-fst.values))
z.scores
#Compute the GIF
K=2
lambda = median(z.scores^2)/qchisq(1/2, df = K-1)
lambda
# compute adjusted p-values from the combined z-scores
adj.p.values = pchisq(z.scores^2/lambda, df = K-1, lower = FALSE)
#histogram of p-values
hist(adj.p.values, col = "red")
# make qqplot for pvalue
qqplot(rexp(length(adj.p.values), rate = log(10)),
       -log10(adj.p.values), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)
## FDR control: Benjamini-Hochberg at level q
plot(-log10(adj.p.values), main="Manhattan plot", xlab = "Locus", cex = .7, col = "grey")
L = length(adj.p.values)
q = 1e-10
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates = order(adj.p.values)[w]
candidates
points(candidates, -log10(adj.p.values)[candidates], pch = 19, cex = .7, col = "red")

#-----------------------------------------------
# load the data set
DT <- fread("ARN_COH_15kb_fold.fst")
print(DT)
# here fst.values is just fst values for snps
fst.values = DT$angsd_Fst
# number of individuals in the comparison
n = 91
# here change the negative fst value into 0
fst.values[fst.values<0] = 0
#uncomment the line below will retain the SNPs with Fst >0
#fst.values = fst.values[fst.values!=0] 
fst.values
K = 2
z.scores = sqrt(fst.values*(n-K)/(1-fst.values))
z.scores
#Compute the GIF
K=2
lambda = median(z.scores^2)/qchisq(1/2, df = K-1)
lambda
# compute adjusted p-values from the combined z-scores
adj.p.values = pchisq(z.scores^2/lambda, df = K-1, lower = FALSE)
#histogram of p-values
hist(adj.p.values, col = "red")
# make qqplot for pvalue
qqplot(rexp(length(adj.p.values), rate = log(10)),
       -log10(adj.p.values), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)
## FDR control: Benjamini-Hochberg at level q
plot(-log10(adj.p.values), main="Manhattan plot", xlab = "Locus", cex = .7, col = "grey")
L = length(adj.p.values)
q = 0.5
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates = order(adj.p.values)[w]
candidates
points(candidates, -log10(adj.p.values)[candidates], pch = 19, cex = .7, col = "red")

