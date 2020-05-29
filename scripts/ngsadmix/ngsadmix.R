logs <- as.data.frame(read.table("logfile"))
logs$K <- c(rep("1", 7), rep("2", 7), rep("3", 7), rep("4", 7))
write.table(logs[, c(2, 1)], "logfile_formatted", row.names = F, 
            col.names = F, quote = F)

# start plot
admix <- t(as.matrix(read.table("chr_k4_run1.qopt")))
head(admix[, 1:98])
ref <- read.table('challenge_98.list', head=TRUE)
populations <- as.vector(ref[,3])
# sort by population
orders <- order(populations, decreasing=FALSE)
admix <- admix[,orders]
populations <- populations[orders]
# draw abline
pops = c("CH1", "CH2", "CH3", "CH4", "REF")
idx = c(rep(0,5))
for (i in c(1:length(populations))){
  for (j in c(1:length(pops))) {
    if (populations[i] == pops[j]){
      idx[j] = i
    }
    
  }
}

orders1 <- order(idx, decreasing=FALSE)
pops <- pops[orders1]
idx <- idx[orders1]
text_loc = c()
runner = 0
for (i in idx){
    text_loc <- c(text_loc, (runner+i)/2)
    runner = i
}

barplot(admix, col = c("blue", "green",'black','red'), space = 0, border = NA, 
        ylab = "Admixture", main = "ch-ref group (K=4)")
text(text_loc, -0.05, unique(populations), 
     xpd = T)
abline(v = idx, lty = 5, lwd = 2, col = "white")


