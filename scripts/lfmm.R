library(lfmm)

## Simulated phenotypes for Arabidopsis thaliana SNP data
data("example.data")
## Simulated (and real) methylation levels for sun-exposed tissue sampled
data("skin.exposure")
## pca plotting
Y <- example.data$genotype
pc <- prcomp(Y)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(6,pc$sdev[6]^2, type = "h", lwd = 3, col = "blue")
## For the A. thaliana samples, the screeplot indicates that there are around K=6 main components in the data. 
## We will use K=6 latent factors in subsequent analyses of the A. thaliana example
Y <- skin.exposure$beta.value
pc <- prcomp(Y)
plot(pc$sdev[1:20]^2, xlab = 'PC', ylab = "Variance explained")
points(2,pc$sdev[2]^2, type = "h", lwd = 3, col = "blue")
## Ridge estimates and GWAS tests for the A.thaliana example
Y <- example.data$genotype
X <- example.data$phenotype #scaled phenotype
## Fit an LFMM, i.e, compute B, U, V estimates
mod.lfmm <- lfmm_ridge(Y = Y, 
                       X = X, 
                       K = 6)
## performs association testing using the fitted model:
pv <- lfmm_test(Y = Y, 
                X = X, 
                lfmm = mod.lfmm, 
                calibrate = "gif")
## The histogram of test significance values is expected to be flat, with a peak near zero. A QQ-plot is displayed as follows
pvalues <- pv$calibrated.pvalue 
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)
# A Manhattan plot can be shown as follows.
plot(-log10(pvalues), 
     pch = 19, 
     cex = .2, 
     xlab = "SNP", ylab = "-Log P",
     col = "grey")
points(example.data$causal.set, 
       -log10(pvalues)[example.data$causal.set], 
       type = "h", 
       col = "blue")