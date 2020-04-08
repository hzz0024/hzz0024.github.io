# usage: Rscript mean_sd.R test 2> /dev/null

fin <- commandArgs(T)
cat("", file=paste(fin,"_mean_sd.log", sep="", collapse=""))

dep <- as.numeric(scan(paste(fin,".depthGlobal", sep = "", collapse=""), what = "char", quiet = T))
cnt = seq(0,length(dep)-1)
Mean = sum((dep * cnt))/(sum(dep))
Deviation = sum((cnt-Mean)^2*dep)/(sum(dep)-1)
Sd = sqrt(Deviation)

cat("\nMean\tDeviation\tSD\n", file=paste(fin,"_mean_sd.log", sep="", collapse=""),append=T)
write.table(cbind(Mean, Deviation, Sd), row.names = F, col.names = F, quote = F, sep = "\t", file=paste (fin,"_mean_sd.log", sep="", collapse=""), append = T)

