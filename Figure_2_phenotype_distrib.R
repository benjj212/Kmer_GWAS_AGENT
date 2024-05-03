### packages
library(vroom)
library(ggplot2)
library(gplots)
library(grid)
library(tidyverse)
### loading the data
setwd("path_to_data")
data <- read.table("Infection_data.txt", header = T, fill=T,as.is=T)
accessions <- as.character(levels(as.factor(data[,1])))
Means_heatmap <- c()
for(i in c(1:length(accessions))){
  Means_heatmap <- rbind(Means_heatmap,c(data[which(data[,1]==accessions[i])[1],1],colMeans(data[which(data[,1]==accessions[i]),c(3:12)])))
}
row.names(Means_heatmap) <- Means_heatmap[,1]
Means_heatmap <- Means_heatmap[,-1]
### removing NA 
Means_heatmap <- as.data.frame(Means_heatmap[-which(is.na(Means_heatmap[,10])),])

Means_heatmap_num <- data.frame(sapply(Means_heatmap, function(x) as.numeric(as.character(x))))
row.names(Means_heatmap_num) <- row.names(Means_heatmap)

for(j in c(1:10)){
  order_plot_acc <- rownames(Means_heatmap_num)[order(Means_heatmap_num[,j])]
  data_order <- c()
  for(i in c(1:length(order_plot_acc))){
    data_order <- rbind(data_order, cbind(i,data[which(data[,1]==order_plot_acc[i]),]))
  }
  
  data_order <- data_order %>% mutate(Line = fct_reorder(Line, data_order[,5]))
  plot <- ggplot(data_order,aes(x=i, y=data_order[,j+3])) +
    stat_summary(fun.data = "mean_cl_boot", geom = "smooth", se = TRUE) +
    xlab("Genotypes") + ylab("Infection rate") +
    scale_x_continuous(breaks = seq(1, max(data_order$i), by = 2)) + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                                                                        axis.text.x=element_blank(),axis.ticks.x=element_blank())
  ggsave(paste("/media/benji/data_benji/005_manuscript/AGENT/Figures/Fig_2/",colnames(Means_heatmap_num)[j], ".png",sep=""),plot=plot,units="px", height = 1000, width=1700)
}

