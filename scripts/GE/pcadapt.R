#install.packages("pcadapt")
library(pcadapt)
library(qvalue)
# 1. reading genotype data (“pcadapt”, “lfmm”, “vcf”, “bed”, “ped”, “pool”)
#path_to_file <- "path_to_directory/foo.lfmm"
#filename <- read.pcadapt(path_to_file, type = "lfmm")
path_to_file <- system.file("extdata", "geno3pops.bed", package = "pcadapt")
filename <- read.pcadapt(path_to_file, type = "bed")

# 2. Choosing the number K of Principal Components
x <- pcadapt(input = filename, K = 20)
plot(x, option = "screeplot")
plot(x, option = "screeplot", K = 10)

# Another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
# With integers
poplist.int <- c(rep(1, 50), rep(2, 50), rep(3, 50))
print(poplist.int)
# With names
poplist.names <- c(rep("POP1", 50),rep("POP2", 50),rep("POP3", 50))
print(poplist.names)
plot(x, option = "scores", pop = poplist.int)
plot(x, option = "scores", pop = poplist.names)
# Looking at population structure beyond K = 2 confirms the results of the scree plot. The third and the fourth principal components do not ascertain population structure anymore.
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)

# 3 Computing the test statistic based on PCA
x <- pcadapt(filename, K = 2)
summary(x)

# 4 Graphical tools
# Manhattan Plot
plot(x , option = "manhattan")
# Q-Q Plot
plot(x, option = "qqplot")
# Histograms of the test statistic and of the p-values
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
# The presence of outliers is also visible when plotting a histogram of the test statistic 𝐷𝑗.
plot(x, option = "stat.distribution")

# 5. Choosing a cutoff for outlier detection
# q-values
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.1
outliers <- which(qval < alpha)
length(outliers)
# Benjamini-Hochberg Procedure
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)
# Bonferroni correction
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)

# 6 Linkage Disequilibrium (LD) thinning
path_to_file <- system.file("extdata", "SSMPG2017.rds", package = "pcadapt")
genotypes <- readRDS(path_to_file)
matrix <- read.pcadapt(genotypes, type = "pcadapt")
res <- pcadapt(matrix, K = 20)
plot(res, option = "screeplot")

res <- pcadapt(matrix, K = 4)
plot(res)

par(mfrow = c(2, 2))
for (i in 1:4)
  plot(res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))

res <- pcadapt(matrix, K = 20, LD.clumping = list(size = 200, thr = 0.1))
plot(res, option = "screeplot")

res <- pcadapt(matrix, K = 2, LD.clumping = list(size = 200, thr = 0.1))
par(mfrow = c(1, 2))
for (i in 1:2)
  plot(res$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
plot(res)


# 7 Detecting local adaptation with pooled sequencing data
# load your own dataset
#pool.data <- read.table("path_to_directory/foo")
#filename <- read.pcadapt(pool.data, type = "pool")

pool.data <- system.file("extdata", "pool3pops", package = "pcadapt")
filename <- read.pcadapt(pool.data, type = "pool")
res <- pcadapt(filename)
#The same as res <- pcadapt(filename,K=2)
summary(res)
plot(res,option="screeplot")
plot(res, option = "manhattan")
padj <- p.adjust(res$pvalues, method = "BH")
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers)
get.pc(res,1:150)->aux
print(aux[,2])

# 8. Miscellaneous
# Run pcadapt on matrices loaded in memory
path_to_file <- system.file("extdata", "SSMPG2017.rds", package = "pcadapt")
genotypes <- readRDS(path_to_file)
print(dim(genotypes))
matrix <- read.pcadapt(genotypes, type = "pcadapt")
res <- pcadapt(matrix, K = 20)
plot(res, option = "screeplot")
#  Association between PCs and outliers
snp_pc <- get.pc(x, outliers)
# Component-wise genome scans
path_to_file <- system.file("extdata", "geno3pops.bed", package = "pcadapt")
filename <- read.pcadapt(path_to_file, type = "bed")
x_cw <- pcadapt(filename, K = 2, method = "componentwise")
summary(x_cw$pvalues)








