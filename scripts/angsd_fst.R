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
##### script for chromosome-wide plots
#DT1 <- fread("ch_ref_1kb_fold.fst")
#print(DT1)
#DT1$chr <- as.numeric(DT1$chr)
#pdf("Mahattan_ch_ref_1kb_fold.pdf",width=15,height=10)
#jepg("Mahattan_ch_ref_1kb_fold_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
#par(mfrow=c(2,1)) 
#manhattan(DT1,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.2),
#          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
#          ylab="ch_ref angsd Fst") #main = "Chromosome",
#dev.off()

##### script for single-SNP plot (due to difficulty in opening the pdf, I export the jepg plot here)
DT <- fread("ch_ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#pdf("Mahattan_ch_ref_single_SNP_fold.pdf",width=15,height=10)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_no16inv_minI8D8maxD16_MQ20_nochr6_invers.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="pos",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 1),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref angsd Fst", cex.lab=1.4) #main = "Chromosome",
dev.off()

setwd("~/Documents/Ryan_workplace/oyster/WGS/WGS_fst/Fst_challenge/unfold/plot/single_snp")
DT <- fread("ch_ref_MQ20_minMAF05_SNPe6.fst")
print(DT)
DT$chr <- as.numeric(DT$chr)
#pdf("Mahattan_ch_ref_single_SNP_unfold.pdf",width=15,height=10)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_single_SNP_unfold.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="pos",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.5),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref angsd Fst", cex.lab=1.4) #main = "Chromosome",
dev.off()

##### script for chromosome-wide plots at 1kb, 5kb, and 15kb (in one window)
setwd("~/Documents/Ryan_workplace/oyster/WGS/WGS_fst/Fst_challenge/fold/plot")

DT1 <- fread("ch_ref_1kb_fold.fst")
DT2 <- fread("ch_ref_5kb_fold.fst")
DT3 <- fread("ch_ref_15kb_fold.fst")
DT1$chr <- as.numeric(DT1$chr)
DT2$chr <- as.numeric(DT2$chr)
DT3$chr <- as.numeric(DT3$chr)
#pdf("Mahattan_ch_ref_fold.pdf",width=15,height=10)
jpeg("Mahattan_ch_ref_fold.jpg", width = 16, height = 9, units = 'in', res = 300)
#png("Mahattan_ch_ref_fold.png", width = 6, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1)) 
manhattan(DT1,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref 1kb Fst", cex.lab=1.4) #main = "Chromosome",

manhattan(DT2,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref 5kb Fst", cex.lab=1.4) #main = "Chromosome",

manhattan(DT3,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref 15kb Fst", cex.lab=1.4) #main = "Chromosome",
dev.off()

setwd("~/Documents/Ryan_workplace/oyster/WGS/WGS_fst/Fst_challenge/unfold/plot")
DT1 <- fread("ch_ref_1kb_unfold.fst")
DT2 <- fread("ch_ref_5kb_unfold.fst")
DT3 <- fread("ch_ref_15kb_unfold.fst")
DT1$chr <- as.numeric(DT1$chr)
DT2$chr <- as.numeric(DT2$chr)
DT3$chr <- as.numeric(DT3$chr)
#pdf("Mahattan_ch_ref_unfold.pdf",width=15,height=10)
jpeg("Mahattan_ch_ref_unfold.jpg", width = 16, height = 9, units = 'in', res = 300)
#png("Mahattan_ch_ref_unfold.png", width = 6, height = 9, units = 'in', res = 300)
par(mfrow=c(3,1)) 
manhattan(DT1,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref 1kb Fst", cex.lab=1.4) #main = "Chromosome",

manhattan(DT2,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref 5kb Fst", cex.lab=1.4) #main = "Chromosome",

manhattan(DT3,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.15),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="ch_ref 15kb Fst", cex.lab=1.4) #main = "Chromosome",
dev.off()
















##### Leo script for wide plots, loop to make 10 individual LARGE pdf files
##### I couldn't figure out how to combine them onto one page, but they are too large to combine anyway


#Necessary lines for all plot scripts below
LG.chr = scan("LG.list",what=" ") # empty quotes indicates character string
DT <- fread("ARN_COH_1kb_fold.fst.bak")
DT_complete <- DT[complete.cases(DT),]
testDT<-DT_complete

##### Leo script for wide plots, no loop
LG <- "NC_035780.1"
currDT <- testDT[testDT$chr == 1,]#using the bracket "]" notation which designates the indices of the data set. The first index is for the rows and the second for the columns. leaving the index for the columns blank indicates that we want currDT to contain all the variables (columns) of the original data frame.
print(currDT)
currDT$chr <- as.character("1")
currDT$chr <- as.numeric(currDT$chr)
pdf("Mahattan ARN_COH_1kb_fold_chr6.pdf",width=15,height=5)
par(mfrow=c(1,1)) 
manhattan(currDT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.2),
          col=c("#3D3D3D","#B0B0B0"),genomewideline=F, suggestiveline=F,
          main = LG, ylab="ARN_COH_1kb_fold")
dev.off()

##### Leo script for wide plots, loop to make 10 individual LARGE pdf files
##### I couldn't figure out how to combine them onto one page, but they are too large to combine anyway

pdf('ARN_COH 1 KB Fst Manhattan by chr', width=15, height=5)
par(mfrow=c(10,1)) 
i=0
for (LG in LG.chr){
  i <- i+1
  currDT <- testDT[testDT$chr == LG,]  #using the bracket "]" notation which designates the indices of the data set. The first index is for the rows and the second for the columns. leaving the index for the columns blank indicates that we want currDT to contain all the variables (columns) of the original data frame.
  currDT$chr <- as.character(i)
  currDT$chr <- as.numeric(currDT$chr)
  
  pdf(file=paste("ARN_COH ", LG, ".pdf", sep=""),width=15,height=5)
  
  manhattan(currDT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 0.20),
            col=c("#3D3D3D","#B0B0B0"),genomewideline=F, suggestiveline=F,
            main = LG, ylab="ARN_COH angsd Fst")
  dev.off()
}
##### make 10 individual square Fst manhattan plots, png format
i=0
for (LG in LG.chr){
  i <- i+1
  currDT <- testDT[testDT$chr == LG,]  #using the bracket "]" notation which designates the indices of the data set. The first index is for the rows and the second for the columns. leaving the index for the columns blank indicates that we want currDT to contain all the variables (columns) of the original data frame.
  currDT$chr <- as.character(i)  # replace with counter number
  currDT$chr <- as.numeric(currDT$chr)  #convert back to numeric
  
  options(device=png)   # plot to screen (used to work)
  par(mar=c(4, 4, 1, 1) + 0.001)
  manhattan(currDT,chr="chr",bp="window_start",p="angsd_Fst",logp=FALSE, ylab="NB_SR angsd Fst",
            main = LG, cex = 0.5, cex.axis = 0.8)
  
  # save .png
  plotname <- sprintf("Mahattan NB_SR angsd Fst chr%d.png", i) # sprintf = string print function with %s for text and %d for digit
  dev.copy(png, filename=plotname)  #copy the contents of the graph window to a file without having to re-enter the commands.
  dev.off()
}


