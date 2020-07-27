require(data.table)
library(ggplot2)
# load the p-value results
DT <- fread("p_value_list_all.txt")
print(DT)
DT$V3 <- as.numeric(DT$V3)
# measure the number of snp with p-value < 0.01
length(DT$V3[DT$V3 < 0.001])
# plot the distribution of p-values across the whole genome
ggplot(data=DT, aes(DT$V3)) + 
  geom_histogram(bins = 200)+
  xlab("p-value")
# obtain the p-value after fdr correction
DT$V3 <- round(p.adjust((DT$V3), "fdr"), 5)
length(DT$V3[DT$V3 < 0.001])

library(ggplot2)
require(data.table)
DT <- fread("pi_peer_global_noout_log_200bwin_sites_corrected_NC_035780.1.txt")
print(DT)

theta = DT[DT$V6 != "NA"]
mean(theta$V6)


