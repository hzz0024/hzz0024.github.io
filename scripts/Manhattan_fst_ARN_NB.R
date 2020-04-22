library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)

setwd("/Volumes/cornell/DelBay19/angsd/WGS_fst/plot_w_mt/ARN_NB")
#combine 3 plot into one 
jpeg("Mahattan_ARN_NB_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
DT <- fread("ARN_NB_1kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#pdf("Mahattan_ARN_NB_1kb_fold.pdf",width=15,height=10)
par(mfrow=c(3,1)) 
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.1),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_NB 1kb Window Fst ", cex.lab=1.4) 
DT <- fread("ARN_NB_5kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#jpeg("Mahattan_ARN_NB_5kb_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.1),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_NB 5kb Window Fst ", cex.lab=1.4) 

DT <- fread("ARN_NB_15kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#jpeg("Mahattan_ARN_NB_15kb_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.1),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ARN_NB 15kb Window Fst ", cex.lab=1.4)
dev.off()

