library(sp)
library(raster)
library(ncdf4)
library(RColorBrewer)
library(sf)
library(tmap)
library(maps)
library(plotrix)
ACCESSION <- read.table("Path_to_coordinates/coordinate_accessions", row.names = NULL, sep="\t", header=T)
accessions_name <- as.character(levels(as.factor(ACCESSION[,2])))
coordinates_acc_only <- c()
for(i in c(1:length(accessions_name))){
  coordinates_acc_only <- rbind(coordinates_acc_only,c(ACCESSION[which(ACCESSION[,2]==accessions_name[i])[1],c(1,2,3,4)], colMeans(ACCESSION[which(ACCESSION[,2]==accessions_name[i]),c(5:14)])))
}

plotting_acc <- as.data.frame(coordinates_acc_only[-which(is.na(coordinates_acc_only[,3])),])
plotting_pheno <- plotting_acc
plotting_pheno <- as.matrix(plotting_pheno)
for(j in c(5:14)){
  plotting_pheno[which(plotting_pheno[,j]<25),j] <- "green"
  plotting_pheno[which(plotting_pheno[,j]>25),j] <- "red"
}
alt=getData("alt", country="Switzerland", path=tempdir())

# plotting Map_1
png(paste("map_full.png",sep=""), 2000, 1400)
par(cex=2.5, mar=c(5,5,1,6), xpd=T)
plot(alt, col=terrain.colors(20), main="", xlim=c(5.9,10.8), xlab= "longitude coordinates", ylab="latitude coordinates")
points(plotting_acc$Coord_2,plotting_acc$Coord_1, pch=16,cex=3,  col=rgb(0,0,0,alpha = 150,maxColorValue = 255))

points(8.538572, 47.371861, cex=3, col="darkred", pch=16)
points(7.595145,47.558590, cex=3, col="darkred", pch=16)
points(7.451654,46.948312, cex=3, col="darkred", pch=16)
points(6.636060,46.516794, cex=3, col="darkred", pch=16)

text(8.538572, 47.321861,"Zürich", cex=3, col="darkred")
text(7.645145,47.508590, "Basel", cex=3, col="darkred")
text(7.451654,46.898312, "Bern", cex=3, col="darkred")
text(6.908060,46.486794, "Lausanne", cex=3, col="darkred")
dev.off()



#### second map: world map 
# World map is available in the maps package
library(maps)
library(dplyr)

# No margin
par(mar=c(0,0,0,0))

# World map
m <- map('world',
         col="#f2f2f2", fill=TRUE, bg="white", lwd=0.05,
         mar=rep(0,4),border=0, ylim=c(-80,80), plot=FALSE
)

# Get coordinates        
Switzerland <- map('world', regions = "Switzerland", plot=FALSE, exact=TRUE)
Switzerland_long <- mean(Switzerland$range[1:2])
Switzerland_lat <- mean(Switzerland$range[3:4])

Switzerland_1 <- map('world', regions = "Switzerland", plot=FALSE, exact=TRUE)
Switzerland_1_long <- mean(Switzerland_1$range[1:2])
Switzerland_1_lat <- mean(Switzerland_1$range[3:4])

Argentina <- map('world', regions = "Argentina", plot=FALSE)    
Argentina_long <- mean(Argentina$range[1:2])
Argentina_lat <- mean(Argentina$range[3:4])

Poland <- map('world', regions = "Poland", plot=FALSE)    
Poland_long <- mean(Poland$range[1:2])
Poland_lat <- mean(Poland$range[3:4])

Poland_1 <- map('world', regions = "Poland", plot=FALSE)    
Poland_1_long <- mean(Poland_1$range[1:2])
Poland_1_lat <- mean(Poland_1$range[3:4])

Kazakhstan <- map('world', regions = "Kazakhstan", plot=FALSE)    
Kazakhstan_long <- mean(Kazakhstan$range[1:2])
Kazakhstan_lat <- mean(Kazakhstan$range[3:4])

Japan <- map('world', regions = "Japan", plot=FALSE)    
Japan_long <- mean(Japan$range[1:2])
Japan_lat <- mean(Japan$range[3:4])

Iran <- map('world', regions = "Iran", plot=FALSE)    
Iran_long <- mean(Iran$range[1:2])
Iran_lat <- mean(Iran$range[3:4])

UK <- map('world', regions = "UK", plot=FALSE)    
UK_long <- mean(UK$range[1:2])
UK_lat <- mean(UK$range[3:4])

Turkey <- map('world', regions = "Turkey", plot=FALSE)    
Turkey_long <- mean(Turkey$range[1:2])
Turkey_lat <- mean(Turkey$range[3:4])

# Data frame 
data <- data.frame(long = c(Switzerland_long, Argentina_long, Turkey_long, UK_long, Iran_long, Japan_long, Poland_long, Kazakhstan_long),
                   lat= c(Switzerland_lat, Argentina_lat, Turkey_lat, UK_lat, Iran_lat, Japan_lat, Poland_lat, Kazakhstan_lat))

# Plotting map2
pdf(paste("map_isolate.pdf",sep=""), 10, 7)
par(mar=c(0,0,0,0))
map('world',
    col="lightgrey", fill=TRUE, bg="white", lwd=0.05,
    mar=rep(0,4),border=0, ylim=c(-80,80)
)
points(x=data$long, y=data$lat, col="slateblue", cex=1, pch=18)
dev.off()

