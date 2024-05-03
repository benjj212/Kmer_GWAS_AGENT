args <- commandArgs(trailingOnly = TRUE)
path_gwas <- as.character(args[1])
name_iso <- c("CHE_96224","CHE_97251","IRN_GOR5","GRB_JIW2","KAZ_1b","THUN12","ARG_4.2","JPN.Chikara")
for(g in c(1:8)){
  isolate <- name_iso[g]
  list_kmers <- read.table(paste(path_GWAS,"/", isolate, "/kmers/",isolate,".kmer_list.txt", sep=""), header=F, as.is=T)
  data_kmers <- cbind(paste(isolate, c(1:length(list_kmers[,1])), sep="_"), list_kmers[,1])
  write.table(data_kmers,paste(path_GWAS,"/", isolate, "/kmers/",isolate,".kmer_list.fasta", sep=""), quote = F, col.names = F, row.names = F)
}