library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)
library(LEA)

# load the data set
DT <- fread("ARN_COH_1kb_fold.fst")
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
q = 1e-5
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates = order(adj.p.values)[w]
candidates
points(candidates, -log10(adj.p.values)[candidates], pch = 19, cex = .7, col = "red")

