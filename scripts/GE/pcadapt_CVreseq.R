install.packages("pcadapt")
devtools::install_github(repo="knausb/vcfR")
install.packages("vcfR")
library(pcadapt)
library(qvalue)
library(vcfR)
# 1. reading genotype data (“pcadapt”, “lfmm”, “vcf”, “bed”, “ped”, “pool”)
path_to_file <- "./LA_prun.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
#path_to_file <- system.file("extdata", "geno3pops.bed", package = "pcadapt")
#filename <- read.pcadapt(path_to_file, type = "bed")

# 2. Choosing the number K of Principal Components
x <- pcadapt(input = filename, K = 10)
plot(x, option = "screeplot")
plot(x, option = "screeplot", K = 5)

# Another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
# With integers
#poplist.int <- c(rep(1, 50), rep(2, 50), rep(3, 50))
#print(poplist.int)
# With names
poplist.names <- c(rep("POP1", 6),rep("POP2", 6), rep("POP3", 6),rep("POP4", 6), rep("POP5", 6))
print(poplist.names)
plot(x, option = "scores", pop = poplist.names)
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)
# Looking at population structure beyond K = 2 confirms the results of the scree plot. The third and the fourth principal components do not ascertain population structure anymore.
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)

# 3 Computing the test statistic based on PCA
x <- pcadapt(filename, K = 5)
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
alpha <- 0.01
outliers <- which(qval < alpha)
length(outliers)
# Benjamini-Hochberg Procedure
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.01
outliers <- which(padj < alpha)
length(outliers)
# Bonferroni correction
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.01
outliers <- which(padj < alpha)
length(outliers)

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


############################# DB_1

library(pcadapt)
library(qvalue)
library(vcfR)
# 1. reading genotype data (“pcadapt”, “lfmm”, “vcf”, “bed”, “ped”, “pool”)
path_to_file <- "./DB_1_prun.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
x <- pcadapt(input = filename, K = 10)
plot(x, option = "screeplot")
plot(x, option = "screeplot", K = 2)
poplist.names <- c(rep("CS", 6),rep("NEH", 6))
print(poplist.names)
plot(x, option = "scores", pop = poplist.names)
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)
x <- pcadapt(filename, K = 2)
summary(x)
plot(x , option = "manhattan")
plot(x, option = "qqplot")
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
plot(x, option = "stat.distribution")
### Bonferroni correction
path_to_file <- "./DB_1_prun.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
x <- pcadapt(filename, K = 2)
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.01DB_1_outliers <- which(padj < alpha)
length(DB_1_outliers)
DB_1_file = 'DB_1_prun.bim'
DB_1_SNP = read.delim(DB_1_file, header = FALSE, sep = "\t", dec = ".")
DB_1_list <- DB_1_SNP$V2[DB_1_outliers]

path_to_file <- "./DB_2_prun.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
x <- pcadapt(filename, K = 2)
padj <- p.adjust(x$pvalues,method="bonferroni")
alpha <- 0.01DB_2_outliers <- which(padj < alpha)
length(DB_2_outliers)
DB_2_file = 'DB_2_prun.bim'
DB_2_SNP = read.delim(DB_2_file, header = FALSE, sep = "\t", dec = ".")
DB_2_list <- DB_2_SNP$V2[DB_2_outliers]

######################################## formal run
library(export)
path_to_file <- "./LA_prun.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
par(mfrow = c(1, 2))
x <- pcadapt(input = filename, K = 10)
plot(x, option = "screeplot")
graph2ppt(file = "LA_scteeplot",width=3.2,height=2.0)
x <- pcadapt(filename, K = 2)
plot(x , option = "manhattan")
graph2ppt(file = "LA_manhattan",width=3.2,height=2.0)
# Q-Q Plot
plot(x, option = "qqplot")
graph2ppt(file = "LA_qqplot",width=3.2,height=2.0)
# pca
poplist.names <- c(rep("OBOYS2", 6),rep("SL", 6))
print(poplist.names)
plot(x, option = "scores", pop = poplist.names)
graph2ppt(file = "LA_pca",width=3.2,height=2.0)
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.1
LA_outliers <- which(padj < alpha)
length(LA_outliers)
# pc loading
#par(mfrow = c(1, 2))
#for (i in 1:2)
#  plot(x$loadings[, i], pch = 19, cex = .1, ylab = paste0("Loadings PC", i), main=paste("No. oultiers",bf,"(bonferroni),",bh,"(FDR)") )
LA_file = 'LA_prun.bim'
LA_SNP = read.delim(LA_file, header = FALSE, sep = "\t", dec = ".")
LA_list <- LA_SNP$V2[LA_outliers]

########### DB_1
library(export)
path_to_file <- "./DB_1_prun.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
par(mfrow = c(1, 2))
x <- pcadapt(input = filename, K = 10)
plot(x, option = "screeplot")
graph2ppt(file = "DB_1_scteeplot",width=3.2,height=2.0)
x <- pcadapt(filename, K = 2)
plot(x , option = "manhattan")
graph2ppt(file = "DB_1_manhattan",width=3.2,height=2.0)
# Q-Q Plot
plot(x, option = "qqplot")
graph2ppt(file = "DB_1_qqplot",width=3.2,height=2.0)
# pca
poplist.names <- c(rep("CS", 6),rep("NEH", 6))
print(poplist.names)
plot(x, option = "scores", pop = poplist.names)
graph2ppt(file = "DB_1_pca",width=3.2,height=2.0)
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.1
DB_1_outliers <- which(padj < alpha)
length(DB_1_outliers)
DB_1_file = 'DB_1_prun.bim'
DB_1_SNP = read.delim(DB_1_file, header = FALSE, sep = "\t", dec = ".")
DB_1_list <- DB_1_SNP$V2[DB_1_outliers]
########### DB_2
library(export)
path_to_file <- "./DB_2_prun.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
par(mfrow = c(1, 2))
x <- pcadapt(input = filename, K = 10)
plot(x, option = "screeplot")
graph2ppt(file = "DB_2_scteeplot",width=3.2,height=2.0)
x <- pcadapt(filename, K = 2)
plot(x , option = "manhattan")
graph2ppt(file = "DB_2_manhattan",width=3.2,height=2.0)
# Q-Q Plot
plot(x, option = "qqplot")
graph2ppt(file = "DB_2_qqplot",width=3.2,height=2.0)
# pca
poplist.names <- c(rep("CS", 6),rep("DEBY", 6))
print(poplist.names)
plot(x, option = "scores", pop = poplist.names)
graph2ppt(file = "DB_2_pca",width=3.2,height=2.0)
padj <- p.adjust(x$pvalues,method="BH")
alpha <- 0.1
DB_2_outliers <- which(padj < alpha)
length(DB_2_outliers)
DB_2_file = 'DB_2_prun.bim'
DB_2_SNP = read.delim(DB_2_file, header = FALSE, sep = "\t", dec = ".")
DB_2_list <- DB_2_SNP$V2[DB_2_outliers]

# read file from bayescan results
Bas_LA = read.table("Bas_LA.txt", sep="\t")
Bas_DB_1 = read.table("Bas_DB_1.txt", sep="\t")
Bas_DB_2 = read.table("Bas_DB_2.txt", sep="\t")
# check the common shared snps
intersect(Bas_DB_1$V1,Bas_DB_2$V1)
intersect(Bas_DB_1$V1,Bas_LA$V1)
intersect(Bas_LA$V1,Bas_DB_2$V1)
intersect(intersect(Bas_LA$V1,Bas_DB_1$V1), Bas_DB_2$V1)

intersect(DB_1_list,Bas_DB_1$V1)
length(intersect(DB_1_list,Bas_DB_1$V1))

intersect(DB_2_list,Bas_DB_2$V1)
length(intersect(DB_2_list,Bas_DB_2$V1))

intersect(LA_list,Bas_LA$V1)
length(intersect(LA_list,Bas_LA$V1))

intersect(DB_1_list,DB_2_list)
length(intersect(DB_1_list,DB_2_list))
intersect(DB_1_list,LA_list)
length(intersect(DB_1_list,LA_list))
intersect(DB_2_list,LA_list)
length(intersect(DB_2_list,LA_list))

intersect(intersect(DB_2_list,LA_list), intersect(DB_1_list,Bas_DB_1$V1))

intersect(intersect(DB_2_list,LA_list), intersect(DB_2_list,Bas_DB_2$V1))

intersect(intersect(DB_2_list,LA_list), intersect(LA_list,Bas_LA$V1))

