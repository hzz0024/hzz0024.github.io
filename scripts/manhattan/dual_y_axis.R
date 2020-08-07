#library
library(latticeExtra)
library(dplyr)
# create data
set.seed(1)

# 8_36715103
pos = 36715103
low = pos - 250000
high = pos + 250000
d1<-read.table("CS.sites.pi", header=T, as.is=T, sep="\t")
d2<-read.table("NEH.sites.pi", header=T, as.is=T, sep="\t")
CS <- subset(d1, CHROM==8 & POS > low & POS < high)$PI
NEH <- subset(d2, CHROM==8 & POS > low & POS < high)$PI
pi <- subset(d2, CHROM==8 & POS > low & POS < high)$POS
data <- data.frame(pi,CS,NEH)
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(pi,CS, pch = 19, col = 2)              # Create first plot
par(new = TRUE)                             # Add new plot
plot(pi,NEH, pch = 17, col = 3,              # Create second plot without axes
     axes = FALSE, xlab = "Position (Chr 8)", ylab = expression(pi))
legend(36900000, 0.1,legend=c("CS (wild)", "NEH (dom)"),
       col=c("red", "green"), pch=c(19,17), cex=1)
abline(v = 36715103, col="red", lwd=1, lty=2)
library(export)
graph2ppt(file = "8_36715103",width=10,height=7)

# 8_36715103
pos = 36715103
low = pos - 250000
high = pos + 250000
d1<-read.table("CS.sites.pi", header=T, as.is=T, sep="\t")
d2<-read.table("DEBY.sites.pi", header=T, as.is=T, sep="\t")
CS <- subset(d1, CHROM==8 & POS > low & POS < high)$PI
DEBY <- subset(d2, CHROM==8 & POS > low & POS < high)$PI
pi <- subset(d2, CHROM==8 & POS > low & POS < high)$POS
data <- data.frame(pi,CS,DEBY)
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(pi,CS, pch = 19, col = 2)              # Create first plot
par(new = TRUE)                             # Add new plot
plot(pi,DEBY, pch = 17, col = 3,              # Create second plot without axes
     axes = FALSE, xlab = "Position (Chr 8)", ylab = expression(pi))
legend(36900000, 0.1,legend=c("CS (wild)", "DEBY (dom)"),
       col=c("red", "green"), pch=c(19,17), cex=1)
abline(v = 36715103, col="red", lwd=1, lty=2)
library(export)
graph2ppt(file = "8_36715103",width=10,height=7)

###################################
# 8_63266137_DB_1
pos = 63266137
low = pos - 100000
high = pos + 100000
d1<-read.table("CS.sites.pi", header=T, as.is=T, sep="\t")
d2<-read.table("NEH.sites.pi", header=T, as.is=T, sep="\t")
CS <- subset(d1, CHROM==8 & POS > low & POS < high)$PI
NEH <- subset(d2, CHROM==8 & POS > low & POS < high)$PI
pi <- subset(d2, CHROM==8 & POS > low & POS < high)$POS
data <- data.frame(pi,CS,NEH)
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(pi,CS, pch = 19, col = 2)              # Create first plot
par(new = TRUE)                             # Add new plot
plot(pi,NEH, pch = 17, col = 3,              # Create second plot without axes
     axes = FALSE, xlab = "Position (Chr 8)", ylab = expression(pi))
legend(pos, 0.1,legend=c("CS (wild)", "NEH (dom)"),
       col=c("red", "green"), pch=c(19,17), cex=1)
abline(v = pos, col="red", lwd=1, lty=2)
library(export)
graph2ppt(file = "8_63266137_DB_1",width=10,height=7)

# 8_63266137_LA
pos = 63266137
low = pos - 100000
high = pos + 100000
d1<-read.table("SL.sites.pi", header=T, as.is=T, sep="\t")
d2<-read.table("OBOYS2.sites.pi", header=T, as.is=T, sep="\t")
SL <- subset(d1, CHROM==8 & POS > low & POS < high)$PI
OBOYS2 <- subset(d2, CHROM==8 & POS > low & POS < high)$PI
pi <- subset(d2, CHROM==8 & POS > low & POS < high)$POS
data <- data.frame(pi,SL,OBOYS2)
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(pi,SL, pch = 19, col = 2)              # Create first plot
par(new = TRUE)                             # Add new plot
plot(pi,OBOYS2, pch = 17, col = 3,              # Create second plot without axes
     axes = FALSE, xlab = "Position (Chr 8)", ylab = expression(pi))
legend(pos, 0.1,legend=c("SL (wild)", "OBOYS2 (dom)"),
       col=c("red", "green"), pch=c(19,17), cex=1)
abline(v = pos, col="red", lwd=1, lty=2)
library(export)
graph2ppt(file = "8_63266137_LA",width=10,height=7)

###################################################
### 5_36361325_DB_1
name = "5_36361325_DB_1"
pos = 36361325
low = pos - 2000000
high = pos + 2000000
d1<-read.table("CS.sites.pi", header=T, as.is=T, sep="\t")
d2<-read.table("NEH.sites.pi", header=T, as.is=T, sep="\t")
CS <- subset(d1, CHROM==5 & POS > low & POS < high)$PI
NEH <- subset(d2, CHROM==5 & POS > low & POS < high)$PI
pi <- subset(d2, CHROM==5 & POS > low & POS < high)$POS
data <- data.frame(pi,CS,NEH)
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(pi,CS, pch = 19, col = 2)              # Create first plot
par(new = TRUE)                             # Add new plot
plot(pi,NEH, pch = 17, col = 3,              # Create second plot without axes
     axes = FALSE, xlab = "Position (Chr 5)", ylab = expression(pi))
legend(pos, 0.1,legend=c("CS (wild)", "NEH (dom)"),
       col=c("red", "green"), pch=c(19,17), cex=1)
abline(v = pos, col="red", lwd=1, lty=2)
library(export)
graph2ppt(file=name,width=10,height=7)


### 5_36361325_LA
name = "5_36361325_LA"
pos = 36361325
low = pos - 2000000
high = pos + 2000000
d1<-read.table("SL.sites.pi", header=T, as.is=T, sep="\t")
d2<-read.table("OBOYS2.sites.pi", header=T, as.is=T, sep="\t")
SL <- subset(d1, CHROM==5 & POS > low & POS < high)$PI
OBOYS2 <- subset(d2, CHROM==5 & POS > low & POS < high)$PI
pi <- subset(d2, CHROM==5 & POS > low & POS < high)$POS
data <- data.frame(pi,SL,OBOYS2)
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(pi,SL, pch = 19, col = 2)              # Create first plot
par(new = TRUE)                             # Add new plot
plot(pi,OBOYS2, pch = 17, col = 3,              # Create second plot without axes
     axes = FALSE, xlab = "Position (Chr 5)", ylab = expression(pi))
legend(pos, 0.1,legend=c("SL (wild)", "OBOYS2 (dom)"),
       col=c("red", "green"), pch=c(19,17), cex=1)
abline(v = pos, col="red", lwd=1, lty=2)
library(export)
graph2ppt(file=name,width=10,height=7)

### 5_36361325_DB_2
name = "5_36361325_DB_2"
pos = 36361325
low = pos - 2000000
high = pos + 2000000
d1<-read.table("CS.sites.pi", header=T, as.is=T, sep="\t")
d2<-read.table("DEBY.sites.pi", header=T, as.is=T, sep="\t")
CS <- subset(d1, CHROM==5 & POS > low & POS < high)$PI
DEBY <- subset(d2, CHROM==5 & POS > low & POS < high)$PI
pi <- subset(d2, CHROM==5 & POS > low & POS < high)$POS
data <- data.frame(pi,CS,DEBY)
par(mar = c(5, 4, 4, 4) + 0.3)              # Additional space for second y-axis
plot(pi,CS, pch = 19, col = 2)              # Create first plot
par(new = TRUE)                             # Add new plot
plot(pi,DEBY, pch = 17, col = 3,              # Create second plot without axes
     axes = FALSE, xlab = "Position (Chr 5)", ylab = expression(pi))
legend(pos, 0.1,legend=c("CS (wild)", "DEBY (dom)"),
       col=c("red", "green"), pch=c(19,17), cex=1)
abline(v = pos, col="red", lwd=1, lty=2)
library(export)
graph2ppt(file=name,width=10,height=7)




# --> construct separate plots for each series
obj1 <- xyplot(CS ~ pi, data, type = "l" , lwd=2)
obj2 <- xyplot(NEH ~ pi, data, type = "l", lwd=2)

plot(pi,CS, pch = 16, col = 2)
plot(pi,NEH, pch = 16, col = 3)
# --> Make the plot with second y axis AND legend:
doubleYScale(obj1, obj2, text = c("obj1", "obj2") , add.ylab2 = TRUE)