
for(MU in c('1e-5','1e-6')){
  for(R in c('1e-5','1e-6','1e-7')){
    filename= paste0('output_m1/',MU,'_',R,'.txt')
    datA = read.table(filename, header=FALSE)$V1
    filename= paste0('output_m2/',MU,'_',R,'.txt')
    datB = read.table(filename, header=FALSE)$V1
    
    b <- min(c(datA,datB))
    e <- max(c(datA,datB))
    ax <- pretty(b:e, n = 12)
    ax<-seq(b-0.05, e+0.05 ,0.01)
  
    hgA <- hist(datA, breaks=ax, plot=FALSE)
    hgB <- hist(datB, breaks=ax, plot=FALSE)
    
    
    c1 <- rgb(173,216,230,max=255,alpha=80,names="lt.blue")
    c2 <- rgb(255,192,203,max=255,alpha=80,names="lt.pink")
    plot(hgA, col=c1, main=paste0('MU:',MU,' R:',R), xlab = "frequency change")
    plot(hgB, col=c2, add=TRUE)
  }
}




