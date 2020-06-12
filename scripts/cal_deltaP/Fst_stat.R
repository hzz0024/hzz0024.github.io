# manhattan plot of Fst
# install.packages("qqman")
# install.packages("caret")
install.packages("animation")

library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)

##### script for single-SNP plot (due to difficulty in opening the pdf, I export the jepg plot here)
DT <- fread("obs.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#pdf("Mahattan_ch_ref_1kb_fold.pdf",width=15,height=10)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_obs.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="V1",bp="V2",p="V3",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.2),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="CH vs REF Fst ", cex.lab=1.4, main = "Observed Fst values with MAF > 0.2",)
dev.off()

DT <- fread("obs.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#pdf("Mahattan_ch_ref_1kb_fold.pdf",width=15,height=10)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_obs.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="V1",bp="V2",p="V3",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.2),
          col=c("blue4","orange3"),genomewideline=0.02996, suggestiveline=F,
          ylab="CH vs REF Fst ", cex.lab=1.4,  main = "Observed Fst values with MAF > 0.2",)

dev.off()

DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_neutral.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="V1",bp="V2",p="V3",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.2),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="CH vs REF Fst ", cex.lab=1.4, main = "Neutral Fst values with MAF > 0.2",)
dev.off()

# plot the Fst distribution

datA = read.table(file = "obs.fst", header=FALSE)$V5
datB = read.table(file = "out.txt", header=FALSE)$V3
    
b <- min(c(datA,datB))
e <- max(c(datA,datB))

ax<-seq(b-0.01, e+0.005,0.005)
    
hgA <- hist(datA, breaks=ax, plot=FALSE)
hgB <- hist(datB, breaks=ax, plot=FALSE)
    
c1 <- rgb(173,216,230,max=255,alpha=80,names="lt.blue")
c2 <- rgb(255,192,203,max=255,alpha=80,names="lt.pink")
plot(hgA, col=c1, xlab="Fst", main=paste0("Observed vs. neutral Fst distribution"))
plot(hgB, col=c2, add=TRUE, xlab="Fst")


##### deltaP for observed data
DT <- fread("ch_ref_obs.txt")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_obs_deltap.jpg", width = 16, height = 9, units = 'in', res = 300)
p1 <- manhattan(DT,chr="chromo",bp="position",p="deltaP",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="CH vs REF deltap ", cex.lab=1.4, main = "Observed dxy",)
p1
dev.off()

##### deltaP for neutral data
DT <- fread("delta_p_out.txt")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_neu_deltap.jpg", width = 16, height = 9, units = 'in', res = 300)
p1 <- manhattan(DT,chr="chromo",bp="position",p="deltaP",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.5),
                col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
                ylab="CH vs REF deltap ", cex.lab=1.4, main = "Observed dxy",)
p1
dev.off()

# plot the allele shift
setwd("~/Documents/Ryan_workplace/DelBay19/11_permutation/deltap/alleleshift")
datA = read.table(file = "obs.txt", header=FALSE)$V3
datB = read.table(file = "neu0.txt", header=FALSE)$V3

b <- min(c(datA,datB))
e <- max(c(datA,datB))

ax<-seq(b-0.01, e+0.01,0.01)

hgA <- hist(datA, breaks=ax, plot=FALSE)
hgB <- hist(datB, breaks=ax, plot=FALSE)

c1 <- rgb(173,216,230,max=255,alpha=80,names="lt.blue")
c2 <- rgb(255,192,203,max=255,alpha=80,names="lt.pink")
plot(hgA, col=c1, xlab="Frequency shift", main=paste0("Frequency shift of observed (blue) and neutral(red) datasets"))
plot(hgB, col=c2, add=TRUE, xlab="Frequency shift")
