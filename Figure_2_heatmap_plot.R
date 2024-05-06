### packages
library(pheatmap) ## for heatmap generation

### loading the data
load("heatmap_fig2.RData")

accessions <- as.character(levels(as.factor(data[,1])))
Means_heatmap <- c()
for(i in c(1:length(accessions))){
  Means_heatmap <- rbind(Means_heatmap,c(data[which(data[,1]==accessions[i])[1],1],colMeans(data[which(data[,1]==accessions[i]),c(3:12)])))
}

#### fully resistant accessions

all_res <- Means_heatmap[which(as.numeric(Means_heatmap[,2])<20 & as.numeric(Means_heatmap[,3])<20 & as.numeric(Means_heatmap[,4])<20 & as.numeric(Means_heatmap[,5])<20 & 
                      as.numeric(Means_heatmap[,6])<20 & as.numeric(Means_heatmap[,7])<20 & as.numeric(Means_heatmap[,8])<20 & as.numeric(Means_heatmap[,9])<20 & 
                      as.numeric(Means_heatmap[,10])<20 & as.numeric(Means_heatmap[,11])<20),]
#### 28 accessions below 20

AG_num_pm4 <- pm4_allele$Line[which(pm4_allele[,8]=="Pm4b")]
Means_heatmap <- cbind(Means_heatmap,25)
for(j in c(1:length(AG_num_pm4))){
  Means_heatmap[which(Means_heatmap[,1]==AG_num_pm4[j]),12] <- 0
}


AG_num <- c()
list_AG_pm2 <- c()
Means_heatmap <- cbind(Means_heatmap,25)
for(i in c(1:length(pm2_allele[,1]))){
  list_AG_pm2 <- c(list_AG_pm2,DAN_num_AG[which(DAN_num_AG[,2]==pm2_allele[i,1]),1])
  AG_temp <- DAN_num_AG[which(DAN_num_AG[,2]==pm2_allele[i,1]),1]
  print(AG_temp)
  AG_num <- c(AG_num,which(winterspring[,1]==AG_temp & winterspring[,2]=="winter"))
  Means_heatmap[which(Means_heatmap[,1]==AG_temp),13] <- 0
}
### nunmber of the fully resistance genotyped as Pm4b of Pm2
### 11 accessions with pm4
for(i in c(1:length(all_res[,1]))){
  print(Means_heatmap[grep(all_res[i,1], Means_heatmap[,1]),12])
}
### 6 accessions with pm2
for(i in c(1:length(all_res[,1]))){
  print(Means_heatmap[grep(all_res[i,1], Means_heatmap[,1]),13])
}

#### adding the field data

field_spring <- read.table("/media/benji/data_benji/005_manuscript/AGENT/Figures/Fig_2/data/field_Spring_wheat.txt", header = T)
field_winter <- read.table("/media/benji/data_benji/005_manuscript/AGENT/Figures/Fig_2/data/field_Winter_wheat.txt", header = T)

Means_heatmap <- cbind(Means_heatmap, 25)
Means_heatmap <- cbind(Means_heatmap, 25)
### scaling 
for(i in c(1:length(Means_heatmap[,1]))){
  if(length(grep(Means_heatmap[i,1], field_spring[,1]))==1){
  Means_heatmap[i,14] <- (field_spring[grep(Means_heatmap[i,1], field_spring[,1]),2]/max(field_spring[,2]))*100
  Means_heatmap[i,15] <- 100
  }
  if(length(grep(Means_heatmap[i,1], field_winter[,1]))==1){
  Means_heatmap[i,14] <- (field_winter[grep(Means_heatmap[i,1], field_winter[,1]),2]/max(field_winter[,2]))*100
  Means_heatmap[i,15] <- 0
  }
}

row.names(Means_heatmap) <- Means_heatmap[,1]
Means_heatmap <- Means_heatmap[,-1]
Means_heatmap <- as.data.frame(Means_heatmap[-which(is.na(Means_heatmap[,10])),])
Means_heatmap_num <- data.frame(sapply(Means_heatmap, function(x) as.numeric(as.character(x))))
colnames(Means_heatmap_num)[11] <- "Pm4b"
colnames(Means_heatmap_num)[12] <- "Pm2"
length(which(Means_heatmap_num[,11]==0))

Means_heatmap <- cbind(Means_heatmap, 0)
for(h in c(1:length(Means_heatmap[,1]))){
  Means_heatmap[h,16] <- info_landrace[grep(rownames(Means_heatmap)[h], info_landrace[,1])[1],2]
}
Means_heatmap

# ordering the heatmap only base on the 10 isolate phenotype
heatmap_order <- pheatmap(Means_heatmap_num[,c(1:10)])
temp_heatmap <- Means_heatmap[heatmap_order[[1]][[3]],]
row.names(Means_heatmap_num) <- ""

### plotting the whole heatmap
png("/media/benji/data_benji/005_manuscript/AGENT/Figures/Fig_2/heatmap_pm4b_pm2.png", height = 1200, width=800)
pheatmap(Means_heatmap_num[heatmap_order[[1]][[3]],], show_rownames=F, fontsize=20,
         cluster_col = FALSE,cluster_rows = FALSE, col=colorRampPalette(c("darkgreen", "orange", "orange", "darkred"))(30))
dev.off()

png("/media/benji/data_benji/005_manuscript/AGENT/Figures/Fig_2/heatmap_dendo.png", height = 1200, width=800)
pheatmap(Means_heatmap_num[,c(1:10)], treeheight_row=100, cluster_col = FALSE, fontsize=20)
dev.off()



