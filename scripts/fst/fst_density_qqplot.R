library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)
library(ggpubr)
library(export)
DT5 <- fread("ch_ref_100_50.fst")
DT7 <- fread("ch_ref_100_70.fst")

# here convert the negative fst value into 0
DT5$angsd_Fst[DT5$angsd_Fst<0] = 0
DT7$angsd_Fst[DT7$angsd_Fst<0] = 0
DT5_fst = DT5$angsd_Fst
DT7_fst = DT7$angsd_Fst
# plot fst distribution
# quantile(fst.values, 0.99) will return the 99% quantile position point
hist(DT7_fst, freq = TRUE, breaks = 1000, xlim = c(0, quantile(DT7_fst, 0.99)), xlab='Fst in 100bp windows', main = "ch_ref_minInd70_minMAPQ25")
hist(DT5_fst, freq = TRUE, breaks = 1000, xlim = c(0, quantile(DT5_fst, 0.99)), xlab='Fst in 100bp windows', main = "ch_ref_minInd50_minMAPQ20")

# start qq-plot
# qexp will return the quantile position for exponential distribution
x = qexp(seq(0,1 ,l=50000))
y1 = DT5_fst
y2 = DT7_fst



plotq <- function(x,y, outfile){
  qqplot(x, y, xlab = "Expected quantile",
         pch = 19, cex = 1)
  # plot the fitting line
  qy=quantile(y,seq(0,1,l=length(x)))
  qx=x
  qx = qx[1:length(qx)-1]
  qy = qy[1:length(qy)-1]
  l.1 <- lm(sort(qy)~sort(qx))
  abline(coef(l.1)[1], coef(l.1)[2])
  print(coef(l.1))
  graph2ppt(file = outfile,width=5.6,height=5.6)
}

plotq(x,y1, "test1")
plotq(x,y2, "test2")

#----------end-------------#


require(vcd)
require(MASS)

# data generation
ex <- rexp(length(fst.values),100) # generate some exponential distribution
control <- fst.values # generate some other distribution

# estimate the parameters
fit1 <- fitdistr(ex, "exponential") 
fit2 <- fitdistr(control, "exponential")

# goodness of fit test
ks.test(ex, "pexp", fit1$estimate) # p-value > 0.05 -> distribution not refused
ks.test(control, "pexp", fit2$estimate) #  significant p-value -> distribution refused

# plot a graph
hist(control, freq = TRUE, breaks = 1000, xlim = c(0, quantile(control, 0.99)))
curve(dexp(x, rate = fit2$estimate), from = 0, col = "red", add = TRUE)




# create qqplot for fst, dataset = ch_ref_1kb_fold.fst
x = qexp(seq(0,1 ,l=10))
y = fst.values
qqplot(x, y, xlab = "Expected quantile",
    pch = 19, cex = 1)

qx = x
qy=quantile(sort(y),seq(0,1,l=length(x)))


for(idx in c(1,8,9)){
  abline(h=qy[idx])
  abline(v=qx[idx])
  axis(1, at=qx[idx], labels=floor(idx/length(qx)*100), cex=.1)
  axis(2, at=qy[idx], labels=floor(idx/length(qx)*100), cex=.1)
}



l.1 <- lm(sort(y)~sort(x))
abline(coef(l.1)[1], coef(l.1)[2])
qqline(qexp(fst.values), distribution=qexp, col="blue", lty=2)


x <- qnorm(seq(0,1,l=1002))  # Theoretical normal quantiles
x <- x[-c(1, length(x))]  # Drop ends because they are -Inf and Inf
y <- rnorm(1000)  # Actual data. 1000 points drawn from a normal distribution
l.1 <- lm(sort(y)~sort(x))
qqplot(x, y, xlab="Theoretical Quantiles", ylab="Actual Quantiles")
abline(coef(l.1)[1], coef(l.1)[2])
