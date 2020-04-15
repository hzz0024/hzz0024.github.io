# manhattan plot of Fst
#install.packages("qqman")
#install.packages("caret")
#install.packages("animation")

library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)

require(data.table)

watterson_estimator <- function(Sn, n) {
  a_n = 0
  for (i in 1:(n-1)) { 
    a_n = a_n + 1/i
  }
  theta = Sn/a_n
  return(theta)
}

theta = watterson_estimator(Sn=1,n=12)
print(theta)



##### script for single-SNP plot (due to difficulty in opening the pdf, I export the jepg plot here)

DT <- fread("plot/NEH_no6inv.thetas.ch10.tsv")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1))
DT$phi <- exp(DT$Pairwise)
jpeg("Mahattan_NEH_no6inv_fold_chr10.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="#Chromo",bp="Pos",p="phi", logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 1),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Tajima Theta", cex.lab=1.4, main = "NEH chr10")
dev.off()

DT <- fread("plot/NEH_no6inv.thetas.ch10.tsv")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1))
DT$phi <- exp(DT$Watterson)
jpeg("Mahattan_NEH_no6inv_fold_chr10_Wt.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="#Chromo",bp="Pos",p="phi", logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 1),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Watterson Theta", cex.lab=1.4, main = "NEH chr10")
dev.off()


##### script for windowed Tajima D plot (due to difficulty in opening the pdf, I export the jepg plot here)

DT <- fread("NEH.thetas.window.idx.pestPG")
print(DT)
max(DT$Tajima, na.rm = TRUE)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1))
jpeg("Mahattan_NEH_no6inv_fold_TajimaD.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="Chr",bp="WinCenter",p="Tajima", logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Tajima D", cex.lab=1.4, main = "NEH TajimaD")
dev.off()

DT <- fread("CS.thetas.window.idx.pestPG")
print(DT)
max(DT$Tajima, na.rm = TRUE)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1))
jpeg("Mahattan_CS_no6inv_fold_TajimaD.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="Chr",bp="WinCenter",p="Tajima", logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Tajima D", cex.lab=1.4, main = "CS TajimaD")
dev.off()

DT <- fread("DEBY.thetas.window.idx.pestPG")
print(DT)
max(DT$Tajima, na.rm = TRUE)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1))
jpeg("Mahattan_DEBY_no6inv_fold_TajimaD.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="Chr",bp="WinCenter",p="Tajima", logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Tajima D", cex.lab=1.4, main = "DEBY TajimaD")
dev.off()

DT <- fread("OBOYS2.thetas.window.idx.pestPG")
print(DT)
max(DT$Tajima, na.rm = TRUE)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1))
jpeg("Mahattan_OBOYS2_no6inv_fold_TajimaD.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="Chr",bp="WinCenter",p="Tajima", logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Tajima D", cex.lab=1.4, main = "OBOYS2 TajimaD")
dev.off()

DT <- fread("SL.thetas.window.idx.pestPG")
print(DT)
max(DT$Tajima, na.rm = TRUE)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1))
jpeg("Mahattan_SL_no6inv_fold_TajimaD.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="Chr",bp="WinCenter",p="Tajima", logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Tajima D", cex.lab=1.4, main = "SL TajimaD")
#dev.off()
