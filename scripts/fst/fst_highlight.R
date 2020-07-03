install.packages("readr") 

library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)

DT <- fread("ch_ref_doSAF_fold.fst")
DT$chr <- as.numeric(DT$chr)
DT$SNP <- paste(DT$chr,DT$pos,sep="_")
OT <- read.delim("outlier.txt", header = FALSE)
OT <- lapply(OT, as.character)
# convert the list into vectors
OT <- unlist(OT, use.names=FALSE)

DT <- DT[DT$chr==5]
#pdf("Mahattan_ch_ref_1kb_fold.pdf",width=15,height=10)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_fold_chr5.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="pos",p="angsd_Fst", logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.2),
          col=c("blue4","orange3"),genomewideline=F, highlight=OT, 
          ylab="CH vs REF Single SNP Fst ", cex.lab=1.4) #main = "Chromosome",
dev.off()