
## loading the packages
library(tidyverse)
library(plotly)
### setting up the working directory and path to the data
setwd("/media/benji/data_benji/002_Analysis/022_circular_plot/script/")
path_to_fai <- "/media/benji/data_benji/005_manuscript/AGENT/scripts_github/data/GWAS_comparison/Genome_Triticum_aestivum.IWGSC.dna.toplevel.fa.fai"
path_to_data <- "/media/benji/data_benji/005_manuscript/AGENT/scripts_github/data/GWAS_comparison/"
coordinate_chr <- read.table(path_to_fai)
### loading the data
CS_kmer <- read.table(paste(path_to_data, "CHE_96224.kmers.full.bed",sep=""), header=T)
SYMattis_kmer <- read.table(paste(path_to_data, "SYMattis.CHE_96224.full.bed",sep=""), header=T)[-1,]
dartseq_snp <- read.table(paste(path_to_data, "CHE_96224.DartSeq.SNPs.txt",sep=""), header=T)
agroscope_snp <- read.table(paste(path_to_data, "CHE_96224.SNPs.chip.txt",sep=""), header=T)
### selecting on the chromosome 2A for ploting
CS_2a <- CS_kmer[grep("2A", CS_kmer[,1]),]
SYM_chr2a <- SYMattis_kmer[grep("2A", SYMattis_kmer[,1]),]
dartseq_snp_2A <- dartseq_snp[grep("2A", dartseq_snp[,1]),]
agroscope_snp_2A <- agroscope_snp[which(agroscope_snp[,2]==4),]

### merge the data frame into one dataframe needed for the plot
A <- CS_2a[,c(2,6)]
B <- cbind(dartseq_snp_2A[,3], -log10(as.numeric(dartseq_snp_2A[,4])))
C <- cbind(agroscope_snp_2A[,3], -log10(as.numeric(agroscope_snp_2A[,4])))
D <- SYM_chr2a[,c(2,6)]

### scaling the different genomes as different genomes are used. The scaling was done based on the end of the chromosome
#to add to kmer CS 6983525 A
#to remove from SYMattis 6368278 D
A[,1] <- A[,1]+6983525
D[,1] <- D[,1]-6368278
colnames(A) <- c("pos", "pval")
colnames(B) <- c("pos", "pval")
colnames(C) <- c("pos", "pval")
colnames(D) <- c("pos", "pval")
### merging the data
E <- merge(A,B, by="pos", all=T)
F_ <- merge(E,C, by="pos", all=T)
G <- merge(F_,D, by="pos", all=T)
# setting up the coordinated for the z axis 
G <- cbind(G,0)
G[which(G[,2]>0),6] <- 2
G[which(G[,3]>0),6] <- 4
G[which(G[,4]>0),6] <- 3
G[which(G[,5]>0),6] <- 1

G <- cbind(G,0)
G[which(G[,2]>0),7] <- G[which(G[,2]>0),2]
G[which(G[,3]>0),7] <- G[which(G[,3]>0),3]
G[which(G[,4]>0),7] <- G[which(G[,4]>0),4]
G[which(G[,5]>0),7] <- G[which(G[,5]>0),5]

colnames(G) <- c("pos", "kmer", "SNP_dart", "SNP_ag","SYM", "coord", "pvalue")
#selecting only the pm4 region
G_pm4 <- G[which(as.numeric(G[,1])>750000000),]
### plotting
fig <- plot_ly(G_pm4, x = ~coord, y = ~pos, z = ~pvalue, color = ~coord, colors = c("black","#542788", "#b35806", "#abd9e9"))
fig <- fig %>% add_markers()
fig <- fig %>% add_trace(x = 1, y = 788738188-6368278, mode = "lines")
fig <- fig %>% add_trace(z = 0, y = 788738188-6368278, mode = "lines")
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Genotype classes'),
                                   yaxis = list(title = 'Coordinate in Mbp'),
                                   zaxis = list(title = '-log10(pvalue)')))
### the 3D plot can then be modify and rotated as wish and the plot saved
fig


