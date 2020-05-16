library(qqman)
library(dplyr)
library(caret)
library(tidyverse)
library(animation)
library(stringr)
require(data.table)

DT <- fread("ch_ref_1kb_fold.fst")
DT
# here convert the negative fst value into 0
DT$angsd_Fst[DT$angsd_Fst<0] = 0

# Change line color and fill color
p <- ggplot(DT, aes(x=angsd_Fst))+
  geom_density(color="darkblue", fill="lightblue", alpha=0.4)+
  labs(x="1kb windowed Fst", y = "Density")
p
# make qqplot for fst
qqplot(rexp(length(fst.values)),
       fst.values, xlab = "Expected quantile",
       pch = 19, cex = .4)

abline(0,1)