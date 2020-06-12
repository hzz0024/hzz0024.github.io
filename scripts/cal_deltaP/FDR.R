require(data.table)
library(ggplot2)
# load the p-values from python script
DT <- fread("out.txt")
print(DT)
DT$V3 <- as.numeric(DT$V3)
DT$V3[DT$V3<0.01]
ggplot(data=DT, aes(DT$V3)) + 
  geom_histogram(bins = 200)+
  xlab("p-value")
p_fdr <- round(p.adjust((DT$V3), "fdr"), 3)

#create a data containing the loci with p-value < 0.01
new = data.frame(DT$V1[DT$V3<0.01], DT$V2[DT$V3<0.01])
names(new) <- c('chr','pos')
write.table(new, file = "my_data.txt", sep = "\t", row.names = FALSE, col.names = FALSE)
#################################






x <- rnorm(50, mean = c(rep(0, 25), rep(3, 25)))
p <- 2*pnorm(sort(-abs(x)))
DT$p = p 
ggplot(data=DT, aes(p)) + 
  geom_histogram()

round(p.adjust(sort(p), "fdr"), 6)

ggplot(data=DT, aes(DT$V3)) + 
  geom_histogram()
p_fdr <- round(p.adjust(sort(DT$V3), "fdr"), 3)

t = seq(0.01,1,length.out=1000)
tt = round(p.adjust(t, "fdr"), 3)
tt

p_fdr[p_fdr<0.1]
p_fdr
