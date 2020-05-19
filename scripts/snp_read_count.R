library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)

DT <- fread("ch_ref_98_pca_minI69D49maxD347_minQ25_minMAF05_SNPe6_nochr56invers_70.pos")
print(length(DT$chr[DT$chr == "10"]))
DT$totDepth <- as.numeric(DT$totDepth)
par(mfrow=c(1,1)) 
jpeg("Mahattan_ch_ref_snpdepth.jpg", width = 16, height = 9, units = 'in', res = 300)
manhattan(DT,chr="chr",bp="pos",p="totDepth",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(0, 500),
          col=c("blue4","orange3"),genomewideline=F, suggestiveline=F,
          ylab="read counts ", cex.lab=1.4) #main = "Chromosome",
dev.off()


#Necessary lines for all plot scripts below
LG.chr = scan("LG.list",what=" ") # empty quotes indicates character string
DT <- fread("ch_ref_98_pca_minI69D49maxD347_minQ25_minMAF05_SNPe6_nochr56invers_70.pos")
DT_complete <- DT[complete.cases(DT),]
testDT<-DT_complete
##### Leo script for wide plots, loop to make 10 individual LARGE pdf files
##### make 10 individual square Fst manhattan plots, png format
i=0
LG.chr = scan("LG.list",what=" ") # empty quotes indicates character string
for (LG in LG.chr){
  i <- i+1
  currDT <- testDT[testDT$chr == i,]  #using the bracket "]" notation which designates the indices of the data set. The first index is for the rows and the second for the columns. leaving the index for the columns blank indicates that we want currDT to contain all the variables (columns) of the original data frame.
  currDT$chr <- as.character(i)  # replace with counter number
  currDT$chr <- as.numeric(currDT$chr)  #convert back to numeric
  
  options(device=png)   # plot to screen (used to work)
  jpeg(paste("ch_ref_98indi_minI69D49maxD347_minQ25_minMAF05_SNPe6_nochr56invers_70_chr", i, ".jpg", sep = ""), width = 16, height = 9, units = 'in', res = 300)
  manhattan(currDT,chr="chr",bp="pos",p="totDepth",logp=FALSE, cex = 0.5, cex.axis = 0.8, ylim = c(50, 500),
            col=c("#3D3D3D","#B0B0B0"),genomewideline=F, suggestiveline=F,
            main = LG, ylab="snp read counts")  
  #plotname <- sprintf("Mahattan_ch_ref_no16inv_minI8D8maxD32_unfold_%d.png", i) # sprintf = string print function with %s for text and %d for digit
  #dev.copy(png, filename=plotname, width = 16, height = 9, units = 'in', res = 300)  #copy the contents of the graph window to a file without having to re-enter the commands.
  dev.off()
}
