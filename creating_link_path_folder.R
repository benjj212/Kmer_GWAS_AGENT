args <- commandArgs(trailingOnly = TRUE)
path_to_fastq <- as.character(args[1])
path_new_folder <- as.character(args[2])
setwd(path_to_fastq)
files <- dir()

for(i in c(1:length(files))){
  name_acc <- strsplit(files[i], ".", fixed=T)[[1]][1]
  setwd(path_new_folder)
  dir.create(name_acc)
  setwd(paste(path_new_folder,"/", name_acc, sep=""))
  link <-  paste(path_to_fastq,"/",files[i], sep="")
  write.table(link, paste(path_new_folder,"/", name_acc,"/",name_acc, ".txt", sep=""), row.names = F, col.names = F, quote = F)
}

