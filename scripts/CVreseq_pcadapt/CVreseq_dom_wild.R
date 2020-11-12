library(OutFLANK)
library(vcfR)
library(bigsnpr)
library(ggplot2)
library(plyr)
source("manhattan.R")

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

# pruning
newpc <- snp_autoSVD(G, infos.chr = CHR, infos.pos = POS)
#which_pruned = attr(newpc, 'subset')
#Gm = toMatrix(G)
#Gm = Gm[,which_pruned]
#G_coded <- add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)

na_idx = c(102120,432510,452617,470761,492280,492281,506875,696154,696155,947788,1045640,1281332,1285977,1292155,1387378,1387379,1489399,1502236,1502237,1566254,1571089,1704001,1704049,1705835,1705836,1706681,1751737,1933633,2015470,2040086,2040108,2086156,2503937,2503981,2556820,2556822,2569353,2650550,2730314,2730315,2730316,2771093,2785790,2842281,2862515,2873343,2978492,3023826,3023827,3076380,3076381,3076382,3126955,3134065,3261668,3428280,3636806,3977518,3977519,4240259,4272699,4374683,4374745,4635250,4635251,4735050,4836671,4836672,4991820,4991990,4991991,4991992,5000999,5095498,5104855,5128374,5189110,5189111,5410894)
Gm = toMatrix(G)
Gm = Gm[,-na_idx]
G_coded = add_code256(big_copy(Gm, type="raw"), code=bigsnpr:::CODE_012)
test <- snp_gc(snp_pcadapt(G_coded, U.row = newpc$u[,1]))
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

jpeg("Mahattan_BP_PC1.jpg", width = 16, height = 9, units = 'in', res = 300)
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
jpeg("naive_PC1_PC2.jpg", width = 16, height = 9, units = 'in', res = 300)
plot(x, option = "scores", pop = poplist.names)
dev.off()
plot(x, option = "scores", i = 3, j = 4, pop = poplist.names)
x <- pcadapt(filename, K = 5)
plot(x, option = "scores", pop = poplist.names)
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
jpeg("Mahattan_naive_PC1-5.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(1,1))
manhattan(daf, chr="CHR",bp="POS", p="Ps", highlight1=outlier_SNP, logp=FALSE, cex.axis = 0.8, ylim = c(0, 160),
          col=c("grey","black"),genomewideline=F, suggestiveline=F,
          ylab="-log10(p-values)", cex.lab=1.4, main = "PCAdapt Naive",)
dev.off()
