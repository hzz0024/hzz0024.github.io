library(ggplot2)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#ff0000", "#ff7878", "#D55E00", "#CC79A7")

name <- "ch_ref_98_pca_minI69D49maxD347_minQ30_minMAF05_SNPe6_nochr56invers_70.ibsMat"
m <- as.matrix(read.table(name))
mds <- cmdscale(as.dist(m))
anno <- read.table('challenge_98.list', sep="\t", header=T)
mds <- as.data.frame(mds)
mds$population <- anno$CLUSTER

#plot(mds,lwd=2,ylab="Dimension 2",xlab="Dimension 1",main="Multidimensional Scaling Plot",col=rep(1:5,each=1))
ggplot() + geom_point(data=mds, aes_string(x='V1', y='V2', color="population"))+ 
ggtitle('Multidimensional Scaling')+
theme(plot.title = element_text(hjust = 0.5))+
xlab("Dimension 1") + ylab("Dimension 2") +
scale_colour_manual(values=cbPalette, breaks=c("CH1","CH2","CH3","CH4","REF"))

ggsave("ch_ref_MDS.jpg")