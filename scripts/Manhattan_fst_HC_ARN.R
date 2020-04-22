library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)

setwd("/Volumes/cornell/DelBay19/angsd/WGS_fst/plot_w_mt/HC_ARN")
#combine 3 plot into one 
jpeg("Mahattan_HC_ARN_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
DT <- fread("HC_ARN_1kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#pdf("Mahattan_HC_ARN_1kb_fold.pdf",width=15,height=10)
par(mfrow=c(3,1)) 
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.1),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="HC_ARN 1kb Window Fst ", cex.lab=1.4) 
DT <- fread("HC_ARN_5kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#jpeg("Mahattan_HC_ARN_5kb_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.1),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="HC_ARN 5kb Window Fst ", cex.lab=1.4) 

DT <- fread("HC_ARN_15kb_fold.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#jpeg("Mahattan_HC_ARN_15kb_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.1),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="HC_ARN 15kb Window Fst ", cex.lab=1.4)
dev.off()

