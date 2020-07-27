library(ggplot2)
require(data.table)
library(plyr)

SL_pi_100 = 'SL_100.windowed.pi'
OBOYS2_pi_100 = 'OBOYS2_100.windowed.pi'
CS_pi_100 = 'CS_100.windowed.pi'
NEH_pi_100 = 'NEH_100.windowed.pi'
DEBY_pi_100 = 'DEBY_100.windowed.pi'

SL_OBOYS2_fst_100 = 'SL_OBOYS2_100.windowed.weir.fst'
CS_NEH_fst_100 = 'CS_NEH_100.windowed.weir.fst'
CS_DEBY_fst_100 = 'CS_DEBY_100.windowed.weir.fst'

SL_pi_100 <- read.delim(SL_pi_100, header = TRUE, sep='\t')
OBOYS2_pi_100 <- read.delim(OBOYS2_pi_100, header = TRUE, sep='\t')
CS_pi_100 <- read.delim(CS_pi_100, header = TRUE, sep='\t')
NEH_pi_100 <- read.delim(NEH_pi_100, header = TRUE, sep='\t')
DEBY_pi_100 <- read.delim(DEBY_pi_100, header = TRUE, sep='\t')

SL_OBOYS2_fst_100 <- read.delim(SL_OBOYS2_fst_100, header = TRUE, sep='\t') 
SL_OBOYS2_fst_100$WEIGHTED_FST[SL_OBOYS2_fst_100$WEIGHTED_FST<0] <- 0
CS_NEH_fst_100 <- read.delim(CS_NEH_fst_100, header = TRUE, sep='\t') 
CS_DEBY_fst_100 <- read.delim(CS_DEBY_fst_100, header = TRUE, sep='\t') 



############################# 5000
SL_pi_5000 = 'SL_5000.windowed.pi'
OBOYS2_pi_5000 = 'OBOYS2_5000.windowed.pi'
CS_pi_5000 = 'CS_5000.windowed.pi'
NEH_pi_5000 = 'NEH_5000.windowed.pi'
DEBY_pi_5000 = 'DEBY_5000.windowed.pi'

SL_pi_5000 <- read.delim(SL_pi_5000, header = TRUE, sep='\t')
OBOYS2_pi_5000 <- read.delim(OBOYS2_pi_5000, header = TRUE, sep='\t')
CS_pi_5000 <- read.delim(CS_pi_5000, header = TRUE, sep='\t')
NEH_pi_5000 <- read.delim(NEH_pi_5000, header = TRUE, sep='\t')
DEBY_pi_5000 <- read.delim(DEBY_pi_5000, header = TRUE, sep='\t')

SL_OBOYS2_fst_5000 = 'SL_OBOYS2_5000.windowed.weir.fst'
CS_NEH_fst_5000 = 'CS_NEH_5000.windowed.weir.fst'
CS_DEBY_fst_5000 = 'CS_DEBY_5000.windowed.weir.fst'

SL_OBOYS2_fst_5000 <- read.delim(SL_OBOYS2_fst_5000, header = TRUE, sep='\t') 
SL_OBOYS2_fst_5000$WEIGHTED_FST[SL_OBOYS2_fst_5000$WEIGHTED_FST<0] <- 0
CS_NEH_fst_5000 <- read.delim(CS_NEH_fst_5000, header = TRUE, sep='\t') 
CS_DEBY_fst_5000 <- read.delim(CS_DEBY_fst_5000, header = TRUE, sep='\t') 

# finding common rows in between two populations see https://stackoverflow.com/questions/30312371/finding-common-rows-in-r
SL_OBOYS2_100 <- join_all(list(SL_pi_100,OBOYS2_pi_100), by = c('CHROM', 'BIN_START'), type = 'inner')
# change the name of headers
colnames(SL_OBOYS2_100) <- c("CHROM", "START_wild", "END_wild", "N1", "P1", "END_dom", "N2", "P2" )
# obtain the ratio
SL_OBOYS2_100$ratio <- SL_OBOYS2_100$P1/SL_OBOYS2_100$P2

SL_OBOYS2_5000 <- join_all(list(SL_pi_5000,OBOYS2_pi_5000), by = c('CHROM', 'BIN_START'), type = 'inner')
# change the name of headers
colnames(SL_OBOYS2_5000) <- c("CHROM", "START_wild", "END_wild", "N1", "P1", "END_dom", "N2", "P2" )
# obtain the ratio
SL_OBOYS2_5000$ratio <- SL_OBOYS2_5000$P1/SL_OBOYS2_5000$P2



######################### Mahattan Plot #########################
# manhattan plot of theta
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

#watterson_estimator <- function(Sn, n) {
#  a_n = 0
#  for (i in 1:(n-1)) { 
#    a_n = a_n + 1/i
#  }
#  theta = Sn/a_n
#  return(theta)
#}

#theta = watterson_estimator(Sn=1,n=12)
#print(theta)

##### script for windowed theta plot (due to difficulty in opening the pdf, I export the jepg plot here)

SL_OBOYS2_100$ratio <- as.numeric(SL_OBOYS2_100$ratio)
par(mfrow=c(1,1))
jpeg("SL_OBOYS2_pi_100.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(SL_OBOYS2_100,chr="CHROM",bp="START_wild",p="ratio", logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab= expression(paste(theta, pi, " wild/", theta, pi, " domestic")), cex.lab=1.4, main = "Louisiana wild (SL) vs domestic (OBOYS2), sliding window = 100 bp")
dev.off()

SL_OBOYS2_5000$ratio <- as.numeric(SL_OBOYS2_5000$ratio)
par(mfrow=c(1,1))
jpeg("SL_OBOYS2_pi_5000.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(SL_OBOYS2_5000,chr="CHROM",bp="START_wild",p="ratio", logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab= expression(paste(theta, pi, " wild/", theta, pi, " domestic")), cex.lab=1.4, main = "Louisiana wild (SL) vs domestic (OBOYS2), sliding window = 5000 bp")
dev.off()


jpeg("SL_OBOYS2_fst_100.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(SL_OBOYS2_fst_100,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST", logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 1),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab= expression(paste("Fst")), cex.lab=1.4, main = "Louisiana wild (SL) vs domestic (OBOYS2), sliding window = 100 bp")
dev.off()

jpeg("SL_OBOYS2_fst_100.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(SL_OBOYS2_fst_100,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST", logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 1),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab= expression(paste("Fst")), cex.lab=1.4, main = "Louisiana wild (SL) vs domestic (OBOYS2), sliding window = 100 bp")
dev.off()

jpeg("SL_OBOYS2_fst_5000.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(SL_OBOYS2_fst_5000,chr="CHROM",bp="BIN_START",p="WEIGHTED_FST", logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 1),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab= expression(paste("Fst")), cex.lab=1.4, main = "Louisiana wild (SL) vs domestic (OBOYS2), sliding window = 5000 bp")
dev.off()
