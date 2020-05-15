# Read values from csv
dt <- read.table("plot.csv", header=T, sep=",")
#dt <- read.table("plot_no_invers.csv", header=T, sep=",")
dt
max_y <- max(dt)
plot_colors <- c("#e41a1c",
                 "#377eb8",
                 "#4daf4a",
                 "#984ea3",
                 "#ff7f00",
                 "#e6beff",
                 "#a65628",
                 "#f781bf",
                 "#808080",
                 "#000080")
jpeg("com_snps.jpg", width = 8, height = 6, units = 'in', res = 300)
#jpeg("com_snps_no_invers.jpg", width = 8, height = 6, units = 'in', res = 300)
plot(NULL, type="o", col=plot_colors[1], 
     ylim=c(0,max_y), xlim=c(1,10), axes=FALSE, ann=FALSE)
# Make x axis using 1-10
axis(1, at=1:10)
# Make y axis with horizontal labels that display ticks at 
# every 4 marks. 4*0:max_y is equivalent to c(0,4,8,12).
axis(2, las=1, at=4*0:max_y)
# Create box around plot
box()
pchs = c(21,22,23,21,22,23,21,22,23,21)
ltys = c(1,2,3,1,2,3,1,2,3,1)

for (i in c(1: length(names(dt)))){
   lines(dt[,i], type="o", pch=pchs[i], lty=ltys[i], col=plot_colors[i],cex=1)
}


# Create a title with a red, bold/italic font
title(main="Common shared SNPs in challenge vs. wild comparsions (including inversions)", col.main="red", font.main=4)

# Label the x and y axes with dark green text
title(xlab= "Chromosome", col.lab=rgb(0,0,0))
title(ylab= "Counts", col.lab=rgb(0,0,0))

# Create a legend at (1, max_y) that is slightly smaller 
# (cex) and uses the same line colors and points used by 
# the actual plots
legend("topright", max_y, names(dt), cex=0.8, col=plot_colors,
       pch=pchs, lty=ltys);

# Turn off device driver (to flush output to png)
dev.off()