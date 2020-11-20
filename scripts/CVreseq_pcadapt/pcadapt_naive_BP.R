library(OutFLANK)
library(vcfR)
library(bigsnpr)
library(ggplot2)
library(plyr)
library(related)
source("manhattan.R")

############################ start pcadapt BP #########################

toMatrix <- function(G){
  Gm = matrix(nrow=length(G[,1]), ncol=length(G[1,]))
  for(i in seq(length(G[,1]))){
    Gm[i,]=G[i,]
  }
  return(Gm)
}

remove.packages("bigsnpr")
install.packages("/Users/ryan/Downloads/bigsnpr", repos = NULL, type="source")
library(bigsnpr)
setwd("~/Documents/Ryan_workplace/CVreseq_clumping")
bedfile = "SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos

### Calculate Relatedness ####
setwd("~/Documents/Ryan_workplace/CVreseq_pcadapt")
vcf <- read.vcfR("SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.all.thin.recode.vcf.gz")
geno <- vcf@gt[,-1] 
genosub <- geno
G_relate <- matrix(NA, nrow=ncol(genosub), ncol=nrow(genosub)*2)
abob <- t(substr(genosub,start = 1, stop=1))
class(abob) <- "numeric"
abob = abob + 1
bbob <- t(substr(genosub,start = 3, stop=3))
class(bbob) <- "numeric"
bbob = bbob + 1
odd <- seq(1, ncol(G_relate), by=2)
G_relate[,odd] <- abob
G_relate[,odd+1] <- bbob
rownames(G_relate) <- rownames(abob)
position = getPOS(vcf)
colnames(G_relate) <- rep(position, each=2)
class(G_relate) <- "numeric"
G_relate0 = G_relate

ind_group = rep(c(1,2,3,4,5),each=6)
G_relate <- data.frame(Individual = rownames(G_relate), population=ind_group, G_relate)
head(G_relate[,1:11])
G_relate2 <- G_relate[,-2]

head(G_relate2[,1:11])
G_relate2$Individual <- as.character(G_relate$Individual)
#G_relate2[,2:ncol(G_relate2)] <- G_relate2[,2:ncol(G_relate2)]+1 # 0s read as missing data
head(G_relate2[,1:11])
library('related')
poprelate2 <- coancestry(G_relate2,  lynchrd =1 )
(rem_ind <- unique(poprelate2$relatedness$ind2.id[which(poprelate2$relatedness$lynchrd>0.5)]))
ind_keep <- which(!(G_relate2$Individual %in% rem_ind))
print(c("final sample size:", length(ind_keep)))
# relatedness plotting
relateness = poprelate2$relatedness
write.table(relateness, file = "relatedness.csv", sep = ",", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
filename = '/out.relatedness2'
dat = read.csv(filename, sep='\t')
relateness  = matrix(relateness, nrow=30)
rownames(relateness) =  dat[,2][1:30]
colnames(relateness) =  dat[,2][1:30]
#heatmap(relateness, Colv = NA, Rowv = NA)
library('glots')
heatmap.2(relateness, trace="none", density.info='none')

# pruning
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, k=10)
#which_pruned = attr(newpc, 'subset')
#Gm = toMatrix(G)
#Gm = Gm[,which_pruned]
#G_coded <- add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)

######### test all vs pruned #########
newpc1 <- big_randomSVD(G)
# result of snp_autoSVD
newpc$u
# result of big_SVD
which_pruned = attr(newpc, 'subset')
Gm = toMatrix(G)
Gm = Gm[,which_pruned]
G_coded <- add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
newpc2 <- big_randomSVD(G_coded, fun.scaling = snp_scaleBinom())
newpc2$u
#######################################
na_idx = c(102120,432510,452617,470761,492280,492281,506875,696154,696155,947788,1045640,1281332,1285977,1292155,1387378,1387379,1489399,1502236,1502237,1566254,1571089,1704001,1704049,1705835,1705836,1706681,1751737,1933633,2015470,2040086,2040108,2086156,2503937,2503981,2556820,2556822,2569353,2650550,2730314,2730315,2730316,2771093,2785790,2842281,2862515,2873343,2978492,3023826,3023827,3076380,3076381,3076382,3126955,3134065,3261668,3428280,3636806,3977518,3977519,4240259,4272699,4374683,4374745,4635250,4635251,4735050,4836671,4836672,4991820,4991990,4991991,4991992,5000999,5095498,5104855,5128374,5189110,5189111,5410894)
Gm = toMatrix(G)
Gm = Gm[,-na_idx]
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:5]))
result_p <- predict(test,log10 = F)
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)

print(CHR[-na_idx][outliers])
#tmp = as.integer(CHR[-na_idx])*104144668+POS[-na_idx]
#(maxp <- max(result, na.rm=TRUE))
#plot(tmp, result,
#     ylim= c(0, maxp*1.19), pch=19, col=rgb(0,0,0,0.5), xaxs="i", yaxs="i",
#     cex=0.5, ylab= "-log10(p-values) PCAdapt", xaxt="n", type="h", bty="l", main="PCAdapt", adj=0)

daf = data.frame(CHR=CHR[-na_idx], POS=POS[-na_idx], SNP=paste0(CHR[-na_idx],"_",POS[-na_idx]), Ps=-log10(result_p))
outlier_SNP = paste0(CHR[-na_idx][outliers],'_',POS[-na_idx][outliers])

jpeg("Mahattan_BP_K3.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 20),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "PCAdapt Best Practice")
dev.off()

############################ start pcadapt naive #########################

library("pcadapt")
path_to_file <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./SNP.MASKED.TRSdp5g75.nDNA.g1.maf05.max2alleles.FIL.format.dom_wild.maf05.nomissing.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
x <- pcadapt(input = filename, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("naive_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("Louisiana_SL", 6),rep("Louisiana_Sel", 6), rep("DelBay_Sel_NEH", 6),rep("Ches_Sel_DEBY", 6), rep("DelBay_CS", 6))
print(poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(filename, K = 3)
jpeg("naive_pca_k3.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
dev.off()
plot(x, option = "scores", i =2, j = 3, pop = poplist.names)
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
jpeg("naive_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x , option = "manhattan")
dev.off()

result_p <- x$pvalues
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("Mahattan_naive_K3.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 160),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "PCAdapt Naive",)
dev.off()



######################### start pairwise ##############################

############################ start pcadapt BP #########################

toMatrix <- function(G){
  Gm = matrix(nrow=length(G[,1]), ncol=length(G[1,]))
  for(i in seq(length(G[,1]))){
    Gm[i,]=G[i,]
  }
  return(Gm)
}

remove.packages("bigsnpr")
install.packages("/Users/ryan/Downloads/bigsnpr", repos = NULL, type="source")
library(bigsnpr)
setwd("~/Documents/Ryan_workplace/CVreseq_pcadapt/")
bedfile = "CS_DEBY.bed"
snp_readBed(bedfile)
# this will create a .rds file
obj.bigSNP <- snp_attach("CS_DEBY.rds")
G <- obj.bigSNP$genotypes
SNPs <- obj.bigSNP$map$marker.ID
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos

# pruning
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS, k=5)
#which_pruned = attr(newpc, 'subset')
#Gm = toMatrix(G)
#Gm = Gm[,which_pruned]
#G_coded <- add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)

na_idx = c(102120,432510,452617,470761,492280,492281,506875,696154,696155,947788,1045640,1281332,1285977,1292155,1387378,1387379,1489399,1502236,1502237,1566254,1571089,1704001,1704049,1705835,1705836,1706681,1751737,1933633,2015470,2040086,2040108,2086156,2503937,2503981,2556820,2556822,2569353,2650550,2730314,2730315,2730316,2771093,2785790,2842281,2862515,2873343,2978492,3023826,3023827,3076380,3076381,3076382,3126955,3134065,3261668,3428280,3636806,3977518,3977519,4240259,4272699,4374683,4374745,4635250,4635251,4735050,4836671,4836672,4991820,4991990,4991991,4991992,5000999,5095498,5104855,5128374,5189110,5189111,5410894)
Gm = toMatrix(G)
Gm = Gm[,-na_idx]
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1:5]))
result_p <- predict(test,log10 = F)
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)

print(CHR[-na_idx][outliers])
#tmp = as.integer(CHR[-na_idx])*104144668+POS[-na_idx]
#(maxp <- max(result, na.rm=TRUE))
#plot(tmp, result,
#     ylim= c(0, maxp*1.19), pch=19, col=rgb(0,0,0,0.5), xaxs="i", yaxs="i",
#     cex=0.5, ylab= "-log10(p-values) PCAdapt", xaxt="n", type="h", bty="l", main="PCAdapt", adj=0)

daf = data.frame(CHR=CHR[-na_idx], POS=POS[-na_idx], SNP=paste0(CHR[-na_idx],"_",POS[-na_idx]), Ps=-log10(result_p))
outlier_SNP = paste0(CHR[-na_idx][outliers],'_',POS[-na_idx][outliers])

jpeg("Mahattan_BP_K3.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 20),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "PCAdapt Best Practice")
dev.off()

############################ start pcadapt naive #########################

library("pcadapt")
path_to_file <- "./CS_DEBY.bed"
filename <- read.pcadapt(path_to_file, type = "bed")
path_to_bim <- "./CS_DEBY.bim"
bim_file = read.delim(path_to_bim, header = FALSE, sep='\t')
# Scree plot
x <- pcadapt(input = filename, K = 10)
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure lie on a steep curve. It is recommended to keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule).
jpeg("naive_screeplot.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "screeplot")
dev.off()
# Score plot: another option to choose the number of PCs is based on the ‘score plot’ that displays population structure.
poplist.names <- c(rep("Louisiana_SL", 6),rep("Louisiana_Sel", 6), rep("DelBay_Sel_NEH", 6),rep("Ches_Sel_DEBY", 6), rep("DelBay_CS", 6))
print(poplist.names)
# Computing the test statistic based on PCA
x <- pcadapt(filename, K = 3)
jpeg("naive_pca_k3.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
dev.off()
plot(x, option = "scores", i =2, j = 3, pop = poplist.names)
# count how many "NA" values in the p-value column
sum(is.na(x$pvalues))
# Manhattan Plot
jpeg("naive_manhattan.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x , option = "manhattan")
dev.off()

result_p <- x$pvalues
padj <- p.adjust(result_p,method="bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
CHR=bim_file$V1
POS=bim_file$V4
outlier_SNP = paste0(CHR[outliers],'_',POS[outliers])
daf = data.frame(CHR=bim_file$V1, POS=bim_file$V4, SNP=bim_file$V2, Ps=-log10(result_p))
daf = daf[complete.cases(daf), ]
jpeg("Mahattan_naive_K3.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 160),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "PCAdapt Naive",)
dev.off()