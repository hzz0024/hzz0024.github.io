library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)

DT <- fread("CH_REF_maf0.05_pctind0.7_cv30_nochr56invers_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_cv30_single.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="pos",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.4, main = "Genome not masked, single SNP")
dev.off()

DT <- fread("CH_REF_maf0.05_pctind0.7_mask_nochr56invers_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_mask_single.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="pos",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.4, main = "Genome masked, single SNP")
dev.off()
##### script for chromosome-wide plots at 1kb, 5kb, and 15kb (in one window)
#setwd("~/Documents/Ryan_workplace/oyster/WGS/WGS_fst/Fst_challenge/fold/plot")

DT1 <- fread("CH_REF_cv30_100_fold.fst")
DT2 <- fread("CH_REF_cv30_500_fold.fst")
DT3 <- fread("CH_REF_cv30_1kb_fold.fst")
DT4 <- fread("CH_REF_mask_100_fold.fst")
DT5 <- fread("CH_REF_mask_500_fold.fst")
DT6 <- fread("CH_REF_mask_1kb_fold.fst")
DT1$chr <- as.numeric(DT1$chr)
DT2$chr <- as.numeric(DT2$chr)
DT3$chr <- as.numeric(DT3$chr)
DT4$chr <- as.numeric(DT4$chr)
DT5$chr <- as.numeric(DT5$chr)
DT6$chr <- as.numeric(DT6$chr)
#pdf("Mahattan_ch_ref_fold.pdf",width=15,height=10)
jpeg("Mahattan_ch_ref_cv30.jpg", width = 16, height = 9, units = 'in', res = 300)
#png("Mahattan_ch_ref_fold.png", width = 6, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1)) 
manhattan(DT1,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.4, main = "Genome not masked, window size from top to bottom: 100bp, 500bp, 1kb")

manhattan(DT2,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.4) #main = "Chromosome",

manhattan(DT3,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.4) #main = "Chromosome",
dev.off()
jpeg("Mahattan_ch_ref_mask.jpg", width = 16, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1)) 
manhattan(DT4,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.4, main = "Genome masked, window size from top to bottom: 100bp, 500bp, 1kb")

manhattan(DT5,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.4) #main = "Chromosome",

manhattan(DT6,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="Fst", cex.lab=1.4) #main = "Chromosome",
dev.off()