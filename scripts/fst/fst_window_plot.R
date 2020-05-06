library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)


jpeg("Mahattan_ch_ref_window.jpg", width = 16, height = 9, units = 'in', res = 300)
DT <- fread("ch_ref_100_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(4,1)) 
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="CH vs. REF 100bp-window Fst ", cex.lab=1.4) 
DT <- fread("ch_ref_500_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="CH vs. REF 500bp-window Fst ", cex.lab=1.4) 
DT <- fread("ch_ref_1kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="CH vs REF 1kb-window Fst ", cex.lab=1.4)
DT <- fread("ch_ref_15kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="CH vs REF 15kb-window Fst ", cex.lab=1.4)
dev.off()

jpeg("Mahattan_ARN_COH_window.jpg", width = 16, height = 9, units = 'in', res = 300)
DT <- fread("ARN_COH_100_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(4,1)) 
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_COH 100bp-window Fst ", cex.lab=1.4) 
DT <- fread("ARN_COH_500_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_COH 500bp-window Fst ", cex.lab=1.4) 
DT <- fread("ARN_COH_1kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_COH 1kb-window Fst ", cex.lab=1.4)
DT <- fread("ARN_COH_15kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_COH 15kb-window Fst ", cex.lab=1.4)
dev.off()


jpeg("Mahattan_ARN_HC_window.jpg", width = 16, height = 9, units = 'in', res = 300)
DT <- fread("ARN_HC_100_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(4,1)) 
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_HC 100bp-window Fst ", cex.lab=1.4) 
DT <- fread("ARN_HC_500_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_HC 500bp-window Fst ", cex.lab=1.4) 
DT <- fread("ARN_HC_1kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_HC 1kb-window Fst ", cex.lab=1.4)
DT <- fread("ARN_HC_15kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_HC 15kb-window Fst ", cex.lab=1.4)
dev.off()


jpeg("Mahattan_ARN_NB_window.jpg", width = 16, height = 9, units = 'in', res = 300)
DT <- fread("ARN_NB_100_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(4,1)) 
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_NB 100bp-window Fst ", cex.lab=1.4) 
DT <- fread("ARN_NB_500_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_NB 500bp-window Fst ", cex.lab=1.4) 
DT <- fread("ARN_NB_1kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_NB 1kb-window Fst ", cex.lab=1.4)
DT <- fread("ARN_NB_15kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_NB 15kb-window Fst ", cex.lab=1.4)
dev.off()


jpeg("Mahattan_ARN_SR_window.jpg", width = 16, height = 9, units = 'in', res = 300)
DT <- fread("ARN_SR_100_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(4,1)) 
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_SR 100bp-window Fst ", cex.lab=1.4) 
DT <- fread("ARN_SR_500_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_SR 500bp-window Fst ", cex.lab=1.4) 
DT <- fread("ARN_SR_1kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_SR 1kb-window Fst ", cex.lab=1.4)
DT <- fread("ARN_SR_15kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_SR 15kb-window Fst ", cex.lab=1.4)
dev.off()


jpeg("Mahattan_COH_HC_window.jpg", width = 16, height = 9, units = 'in', res = 300)
DT <- fread("COH_HC_100_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(4,1)) 
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="COH_HC 100bp-window Fst ", cex.lab=1.4) 
DT <- fread("COH_HC_500_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="COH_HC 500bp-window Fst ", cex.lab=1.4) 
DT <- fread("COH_HC_1kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="COH_HC 1kb-window Fst ", cex.lab=1.4)
DT <- fread("COH_HC_15kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="COH_HC 15kb-window Fst ", cex.lab=1.4)
dev.off()


jpeg("Mahattan_COH_NB_window.jpg", width = 16, height = 9, units = 'in', res = 300)
DT <- fread("COH_NB_100_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(4,1)) 
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="COH_NB 100bp-window Fst ", cex.lab=1.4) 
DT <- fread("COH_NB_500_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="COH_NB 500bp-window Fst ", cex.lab=1.4) 
DT <- fread("COH_NB_1kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="COH_NB 1kb-window Fst ", cex.lab=1.4)
DT <- fread("COH_NB_15kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="COH_NB 15kb-window Fst ", cex.lab=1.4)
dev.off()


jpeg("Mahattan_COH_SR_window.jpg", width = 16, height = 9, units = 'in', res = 300)
DT <- fread("COH_SR_100_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(4,1)) 
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="COH_SR 100bp-window Fst ", cex.lab=1.4) 
DT <- fread("COH_SR_500_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="COH_SR 500bp-window Fst ", cex.lab=1.4) 
DT <- fread("COH_SR_1kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="COH_SR 1kb-window Fst ", cex.lab=1.4)
DT <- fread("COH_SR_15kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="COH_SR 15kb-window Fst ", cex.lab=1.4)
dev.off()


jpeg("Mahattan_HC_SR_window.jpg", width = 16, height = 9, units = 'in', res = 300)
DT <- fread("HC_SR_100_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(4,1)) 
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="HC_SR 100bp-window Fst ", cex.lab=1.4) 
DT <- fread("HC_SR_500_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="HC_SR 500bp-window Fst ", cex.lab=1.4) 
DT <- fread("HC_SR_1kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="HC_SR 1kb-window Fst ", cex.lab=1.4)
DT <- fread("HC_SR_15kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="HC_SR 15kb-window Fst ", cex.lab=1.4)
dev.off()


jpeg("Mahattan_NB_SR_window.jpg", width = 16, height = 9, units = 'in', res = 300)
DT <- fread("NB_SR_100_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
par(mfrow=c(4,1)) 
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="NB_SR 100bp-window Fst ", cex.lab=1.4) 
DT <- fread("NB_SR_500_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="NB_SR 500bp-window Fst ", cex.lab=1.4) 
DT <- fread("NB_SR_1kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="NB_SR 1kb-window Fst ", cex.lab=1.4)
DT <- fread("NB_SR_15kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.25),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="NB_SR 15kb-window Fst ", cex.lab=1.4)
dev.off()

