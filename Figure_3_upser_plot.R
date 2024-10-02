
library(VennDiagram)
library(UpSetR)
library(tidyverse)

### comparing the kmer shared between isolates 
isolates <- c("CHE_96224", "CHE_97251", "GRB_JIW2", "ARG_4.2", "IRN_GOR5", "JPN.Chikara", "KAZ_1b", "THUN12", "TUR_1C")
path_to_data <- "/media/benji/data_benji/005_manuscript/AGENT/scripts_github/data/kmers/"
y <- list(
  "CHE_96224"= as.character(read.table(paste(path_to_data,isolates[1],"/kmers/", isolates[1],".kmer_list.txt", sep=""))[,1]),
  "CHE_97251"= as.character(read.table(paste(path_to_data,isolates[2],"/kmers/", isolates[2],".kmer_list.txt", sep=""))[,1]),
  "GRB_JIW2"= as.character(read.table(paste(path_to_data,isolates[3],"/kmers/", isolates[3],".kmer_list.txt", sep=""))[,1]),
  "ARG_4.2"= as.character(read.table(paste(path_to_data,isolates[4],"/kmers/", isolates[4],".kmer_list.txt", sep=""))[,1]),
  "IRN_GOR5"= as.character(read.table(paste(path_to_data,isolates[5],"/kmers/", isolates[5],".kmer_list.txt", sep=""))[,1]),
  "JPN.Chikara"= as.character(read.table(paste(path_to_data,isolates[6],"/kmers/", isolates[6],".kmer_list.txt", sep=""))[,1]),
  "KAZ_1b"= as.character(read.table(paste(path_to_data,isolates[7],"/kmers/", isolates[7],".kmer_list.txt", sep=""))[,1]),
  "THUN12"= as.character(read.table(paste(path_to_data,isolates[8],"/kmers/", isolates[8],".kmer_list.txt", sep=""))[,1]),
  "TUR_1C"= as.character(read.table(paste(path_to_data,isolates[9],"/kmers/", isolates[9],".kmer_list.txt", sep=""))[,1])
)

# Plot
pdf("Upset_plot.pdf", width = 12, height = 8)
upset(fromList(y), 
      nintersects = 30, 
      nsets = 10, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 2, 
      point.size = 2.8,
      line.size = 1, main.bar.color=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6"))
dev.off()
    
df2 <- data.frame(gene=unique(unlist(y)))
df1 <- lapply(y,function(x){
  data.frame(gene = x)
}) %>% 
  bind_rows(.id = "path")

df_int <- lapply(df2$gene,function(x){
  # pull the name of the intersections
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  
  # build the dataframe
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

df_int %>% 
  group_by(int) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

first_selected <- df_int[which(test[,"TUR_1C"]==0 & test[,"CHE_96224"]==0 & test[,"IRN_GOR5"]==0 & test[,"ARG_4.2"]==1 & test[,"JPN.Chikara"]==1 & test[,"GRB_JIW2"]==1 & 
                                 test[,"CHE_97251"]==1 & test[,"THUN12"]==1 & test[,"KAZ_1b"]==1),1]
second_selected <- df_int[which(test[,"TUR_1C"]==0 & test[,"CHE_96224"]==0 & test[,"IRN_GOR5"]==0 & test[,"ARG_4.2"]==0 & test[,"JPN.Chikara"]==0 & test[,"GRB_JIW2"]==1 & 
                                  test[,"CHE_97251"]==1 & test[,"THUN12"]==1 & test[,"KAZ_1b"]==1),1]
third_selected <- df_int[which(test[,"TUR_1C"]==0 & test[,"CHE_96224"]==0 & test[,"IRN_GOR5"]==0 & test[,"ARG_4.2"]==0 & test[,"JPN.Chikara"]==1 & test[,"GRB_JIW2"]==1 & 
                                 test[,"CHE_97251"]==1 & test[,"THUN12"]==1 & test[,"KAZ_1b"]==1),1]
fourth_selected <- df_int[which(test[,"TUR_1C"]==0 & test[,"CHE_96224"]==0 & test[,"IRN_GOR5"]==0 & test[,"ARG_4.2"]==0 & test[,"JPN.Chikara"]==0 & test[,"GRB_JIW2"]==1 & 
                                  test[,"CHE_97251"]==1 & test[,"THUN12"]==0 & test[,"KAZ_1b"]==1),1]
five_selected <- df_int[which(test[,"TUR_1C"]==0 & test[,"CHE_96224"]==0 & test[,"IRN_GOR5"]==0 & test[,"ARG_4.2"]==0 & test[,"JPN.Chikara"]==0 & test[,"GRB_JIW2"]==0 & 
                                test[,"CHE_97251"]==0 & test[,"THUN12"]==0 & test[,"KAZ_1b"]==1),1]
six_selected <- df_int[which(test[,"TUR_1C"]==1 & test[,"CHE_96224"]==1 & test[,"IRN_GOR5"]==1 & test[,"ARG_4.2"]==1 & test[,"JPN.Chikara"]==1 & test[,"GRB_JIW2"]==0 & 
                               test[,"CHE_97251"]==0 & test[,"THUN12"]==1 & test[,"KAZ_1b"]==0),1]
seven_selected <- df_int[which(test[,"TUR_1C"]==0 & test[,"CHE_96224"]==0 & test[,"IRN_GOR5"]==1 & test[,"ARG_4.2"]==0 & test[,"JPN.Chikara"]==0 & test[,"GRB_JIW2"]==0 & 
                                 test[,"CHE_97251"]==0 & test[,"THUN12"]==0 & test[,"KAZ_1b"]==0),1]
eight_selected <- df_int[which(test[,"TUR_1C"]==1 & test[,"CHE_96224"]==1 & test[,"IRN_GOR5"]==1 & test[,"ARG_4.2"]==0 & test[,"JPN.Chikara"]==0 & test[,"GRB_JIW2"]==0 & 
                                 test[,"CHE_97251"]==0 & test[,"THUN12"]==1 & test[,"KAZ_1b"]==0),1]
nine_selected <- df_int[which(test[,"TUR_1C"]==0 & test[,"CHE_96224"]==0 & test[,"IRN_GOR5"]==0 & test[,"ARG_4.2"]==0 & test[,"JPN.Chikara"]==0 & test[,"GRB_JIW2"]==1 & 
                                test[,"CHE_97251"]==0 & test[,"THUN12"]==0 & test[,"KAZ_1b"]==1),1]
ten_selected <- df_int[which(test[,"TUR_1C"]==0 & test[,"CHE_96224"]==0 & test[,"IRN_GOR5"]==0 & test[,"ARG_4.2"]==0 & test[,"JPN.Chikara"]==0 & test[,"GRB_JIW2"]==0 & 
                               test[,"CHE_97251"]==1 & test[,"THUN12"]==0 & test[,"KAZ_1b"]==1),1]
eleven_selected <- df_int[which(test[,"TUR_1C"]==0 & test[,"CHE_96224"]==1 & test[,"IRN_GOR5"]==1 & test[,"ARG_4.2"]==0 & test[,"JPN.Chikara"]==0 & test[,"GRB_JIW2"]==0 & 
                                  test[,"CHE_97251"]==0 & test[,"THUN12"]==0 & test[,"KAZ_1b"]==0),1]
twelve_selected <- df_int[which(test[,"TUR_1C"]==0 & test[,"CHE_96224"]==0 & test[,"IRN_GOR5"]==0 & test[,"ARG_4.2"]==0 & test[,"JPN.Chikara"]==0 & test[,"GRB_JIW2"]==0 & 
                                  test[,"CHE_97251"]==0 & test[,"THUN12"]==1 & test[,"KAZ_1b"]==0),1]
thirteen_selected <- df_int[which(test[,"TUR_1C"]==0 & test[,"CHE_96224"]==0 & test[,"IRN_GOR5"]==0 & test[,"ARG_4.2"]==1 & test[,"JPN.Chikara"]==1 & test[,"GRB_JIW2"]==0 & 
                                    test[,"CHE_97251"]==0 & test[,"THUN12"]==0 & test[,"KAZ_1b"]==0),1]
fourteen_selected <- df_int[which(test[,"TUR_1C"]==0 & test[,"CHE_96224"]==1 & test[,"IRN_GOR5"]==0 & test[,"ARG_4.2"]==0 & test[,"JPN.Chikara"]==0 & test[,"THUN12"]==0 & 
                                    test[,"GRB_JIW2"]==0 & test[,"KAZ_1b"]==0 & test[,"CHE_97251"]==0),1]


write.table(first_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/first_bar.txt", quote=F, col.names = F, row.names = F)
write.table(second_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/second_bar.txt", quote=F, col.names = F, row.names = F)
write.table(third_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/third_bar.txt", quote=F, col.names = F, row.names = F)
write.table(fourth_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/fourth_bar.txt", quote=F, col.names = F, row.names = F)
write.table(five_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/five_bar.txt", quote=F, col.names = F, row.names = F)
write.table(six_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/six_bar.txt", quote=F, col.names = F, row.names = F)
write.table(seven_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/seven_bar.txt", quote=F, col.names = F, row.names = F)
write.table(eight_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/eight_bar.txt", quote=F, col.names = F, row.names = F)
write.table(nine_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/nine_bar.txt", quote=F, col.names = F, row.names = F)
write.table(ten_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/ten_bar.txt", quote=F, col.names = F, row.names = F)
write.table(eleven_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/eleven_bar.txt", quote=F, col.names = F, row.names = F)
write.table(twelve_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/twelve_bar.txt", quote=F, col.names = F, row.names = F)
write.table(thirteen_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/thirteen_bar.txt", quote=F, col.names = F, row.names = F)
write.table(fourteen_selected, "/media/benji/data_benji/005_manuscript/AGENT/GWAS_data/Upset_plot/fourteen_bar.txt", quote=F, col.names = F, row.names = F)
