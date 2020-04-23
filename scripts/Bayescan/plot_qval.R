library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)

DT <- fread("DB_2_bayescan_all.txt")
print(DT)
DT$chr <- as.numeric(DT$chr)
DT$qval <- -log10(DT$qval)
par(mfrow=c(1,1)) 
jpeg("Bayescan_DB_2.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="pos",p="qval",logp=FALSE, cex = 1.5, cex.axis = 0.8, ylim = c(0, 3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=T, 
          ylab="DB_2 -log10(qval) ", cex.lab=1.4) 
dev.off()

DT <- fread("DB_1_bayescan_all.txt")
print(DT)
DT$chr <- as.numeric(DT$chr)
DT$qval <- -log10(DT$qval)
par(mfrow=c(1,1)) 
jpeg("Bayescan_DB_1.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="pos",p="qval",logp=FALSE, cex = 1.5, cex.axis = 0.8, ylim = c(0, 3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=T,
          ylab="DB_1 -log10(qval) ", cex.lab=1.4)
dev.off()

DT <- fread("LA_bayescan_all.txt")
print(DT)
DT$chr <- as.numeric(DT$chr)
DT$qval <- -log10(DT$qval)
par(mfrow=c(1,1)) 
jpeg("Bayescan_LA.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="pos",p="qval",logp=FALSE, cex = 1.5, cex.axis = 0.8, ylim = c(0, 3),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=T,
          ylab="LA -log10(qval) ", cex.lab=1.4)
dev.off()