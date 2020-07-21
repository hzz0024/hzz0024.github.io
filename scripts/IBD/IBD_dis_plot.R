library(ggplot2)
IBD <- read.csv("IBD_dis_plot.csv", header = TRUE)
IBD$Fst
ggplot(data = IBD, aes(x = Distance, y = Fst)) +
  geom_point(shape=20, fill="red", color="black", size=5)+
  theme(axis.title.y=element_text(size=20))+
  theme(axis.title.x=element_text(size=20))+
  scale_y_continuous(name="Fst/(1-Fst)", limits=c(0, 0.1)) +
  scale_x_continuous(name="Distance(km)", limits=c(0, 30))+
  geom_abline(intercept = 0.068903807201, slope =-0.00035320186051, color="red", 
              linetype="dashed", size=2)
