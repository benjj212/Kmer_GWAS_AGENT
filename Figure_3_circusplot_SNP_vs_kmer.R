library("circlize")
setwd("/media/benji/data_benji/002_Analysis/022_circular_plot/script/")
path_to_fai <- "/media/benji/data_benji/005_manuscript/AGENT/scripts_github/data/GWAS_comparison/Genome_Triticum_aestivum.IWGSC.dna.toplevel.fa.fai"
path_to_data <- "/media/benji/data_benji/005_manuscript/AGENT/scripts_github/data/GWAS_comparison/"
coordinate_chr <- read.table(path_to_fai)
### loading the data
CS_kmer <- read.table(paste(path_to_data, "CHE_96224.kmers.txt",sep=""), header=T)[-1,]
dartseq_snp <- read.table(paste(path_to_data, "CHE_96224.DartSeq.SNPs.txt",sep=""), header=T)
agroscope_snp <- read.table(paste(path_to_data, "CHE_96224.SNPs.chip.txt",sep=""), header=T)
#saving the plot
pdf(paste(path_to_data,"circus_manhattan_plot.pdf",sep=""), height = 7, width = 7)
## creating the circus plot 
chrom_name <- as.character(coordinate_chr[,1])
cytoband.wheat <- as.data.frame(cbind(paste("Chr",chrom_name,sep=""),-50000000, coordinate_chr[,2]+50000000, "a", "gneg"))
circos.par("clock.wise"=TRUE,start.degree=90)
circos.initializeWithIdeogram(cytoband.wheat)
## changing chr
dartseq_snp[,1] <- substr(dartseq_snp[,1], 4,5)
## changing_chr
chr_list <- c("1A", "1B", "1D", "2A", "2B", "2D","3A", "3B", "3D","4A", "4B", "4D","5A", "5B", "5D","6A", "6B", "6D","7A", "7B", "7D")
for(i in c(1:21)){
  agroscope_snp[which(agroscope_snp[,2]==i),2] <- chr_list[i]
}
agro_circus <- agroscope_snp[,c(2,3,4)]
agro_circus[,3] <- -log10(agro_circus[,3])
dartseq_circus <- dartseq_snp[,c(1,3,4)]
dartseq_circus[,3] <- -log10(dartseq_circus[,3])
CS_circus <- CS_kmer[,c(1,2,4)]
CS_circus[,3] <- -log10(CS_circus[,3])
circus_all <- list(agro_circus, dartseq_circus, CS_circus)
#setting up each color 
colors_list <- c("#b35806", "#abd9e9", "#542788")
## plotting the SNPs and kmer in each sectors
for(i in c(3,1,2)){
X <- circus_all[[i]]
  circos.track(factors=NULL, ylim=c(0,20),track.height=0.2)
  Chr_in_file <- as.character(levels(as.factor(X[,1])))
  for(j in c(1:length(Chr_in_file))){
    selected_X <- X[which(X[,1]==Chr_in_file[j]),]
    circos.points(x=selected_X[,2], y=selected_X[,3], sector.index = paste("Chr",Chr_in_file[j],sep=""), pch=16, cex=0.6, col = colors_list[i])
  }
}
dev.off()

