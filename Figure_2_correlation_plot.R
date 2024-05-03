### packages 
library(gclus)


### correlation plot 
setwd("Path_to_data")
data <- read.table("Infection_data.txt", header = T, fill=T,as.is=T)
accessions <- as.character(levels(as.factor(data[,1])))
Means_heatmap <- c()
for(i in c(1:length(accessions))){
  Means_heatmap <- rbind(Means_heatmap,c(data[which(data[,1]==accessions[i])[1],1],colMeans(data[which(data[,1]==accessions[i]),c(3:12)])))
}
row.names(Means_heatmap) <- Means_heatmap[,1]
Means_heatmap <- Means_heatmap[,-1]
Means_heatmap <- Means_heatmap[-which(is.na(Means_heatmap[,10])),]

class(Means_heatmap) <- "numeric"

shoot.r <- abs(abs(cor(Means_heatmap))-1)
shoot.col <- dmat.color(shoot.r, colors=heat.colors(10, alpha=0.5))
png("/media/benji/data_benji/005_manuscript/AGENT/Figures/Fig_2/correlation_plot.png", width = 1000, height = 1000)
cpairs(Means_heatmap, panel.colors=shoot.col, gap=.5,pch=16)
dev.off()
### generating scale for legend

plot(0, xlim=c(1,10), ylim=c(0,11))
points(cbind(1, c(1:10)), pch=15, cex=4, col=heat.colors(10))
shoot.r[order(shoot.r)]
