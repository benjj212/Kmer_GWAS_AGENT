shinyServer(function(input, output, session) {
library(leaflet)
library(ggplot2)
library(tidyverse)
Switzerland <- maps::map('world', regions = "Switzerland", plot=FALSE, exact=TRUE)
Switzerland_long <- mean(Switzerland$range[1:2])
Switzerland_lat <- mean(Switzerland$range[3:4])

Switzerland_1 <- maps::map('world', regions = "Switzerland", plot=FALSE, exact=TRUE)
Switzerland_1_long <- mean(Switzerland_1$range[1:2])-0.5
Switzerland_1_lat <- mean(Switzerland_1$range[3:4])

Argentina <- maps::map('world', regions = "Argentina", plot=FALSE)    
Argentina_long <- mean(Argentina$range[1:2])
Argentina_lat <- mean(Argentina$range[3:4])

Poland <- maps::map('world', regions = "Poland", plot=FALSE)    
Poland_long <- mean(Poland$range[1:2])
Poland_lat <- mean(Poland$range[3:4])

Poland_1 <- maps::map('world', regions = "Poland", plot=FALSE)    
Poland_1_long <- mean(Poland_1$range[1:2])-0.5
Poland_1_lat <- mean(Poland_1$range[3:4])

Kazakhstan <- maps::map('world', regions = "Kazakhstan", plot=FALSE)    
Kazakhstan_long <- mean(Kazakhstan$range[1:2])
Kazakhstan_lat <- mean(Kazakhstan$range[3:4])

Japan <- maps::map('world', regions = "Japan", plot=FALSE)    
Japan_long <- mean(Japan$range[1:2])
Japan_lat <- mean(Japan$range[3:4])

Ireland <- maps::map('world', regions = "Iran", plot=FALSE)    
Ireland_long <- mean(Ireland$range[1:2])
Ireland_lat <- mean(Ireland$range[3:4])

UK <- maps::map('world', regions = "UK", plot=FALSE)    
UK_long <- mean(UK$range[1:2])
UK_lat <- mean(UK$range[3:4])

Turkey <- maps::map('world', regions = "Turkey", plot=FALSE)    
Turkey_long <- mean(Turkey$range[1:2])
Turkey_lat <- mean(Turkey$range[3:4])

# Data frame 
locations <- data.frame(
  Location=c("Switzerland","Switzerland_1", "Argentina", "Turkey", "United Kingdom", "Iran", "Japan", "Poland", "Poland_1", "Kazakhstan"),
  Latitude= c(Switzerland_lat,Switzerland_1_lat, Argentina_lat, Turkey_lat, UK_lat, Ireland_lat, Japan_lat, Poland_lat,Poland_1_lat, Kazakhstan_lat),
  Longitude = c(Switzerland_long,Switzerland_1_long, Argentina_long, Turkey_long, UK_long, Ireland_long, Japan_long, Poland_long,Poland_1_long, Kazakhstan_long))
## loading phenotype data

data <- read.table("20221024_AGENT_inf_test_final_data.txt", header = T, fill=T,as.is=T)
acc <- levels(as.factor(data[,1]))
iso_name <- matrix(data = colnames(data)[c(3:12)], nrow = 10, ncol = 1)
row.names(iso_name) <- c("Switzerland","Switzerland_1","Iran", "United Kingdom", "Kazakhstan","Poland",  "Argentina", "Japan","Poland_1","Turkey")


all_data_new <- c()
for(i in c(1:length(levels(as.factor(data[,1]))))){
  new_data_frame <- c()
  for(j in c(3:12)){
    new_data_frame <- c(new_data_frame,c(mean(data[which(data[,1]==acc[i]),j]),sd(data[which(data[,1]==acc[i]),j])))
  }
  all_data_new <- rbind(all_data_new, new_data_frame)
}
colnames(all_data_new) <- c(1:20)
colnames(all_data_new) <- c("Switzerland_mean","Switzerland_sd","Switzerland_1_mean","Switzerland_1_sd", "Argentina_mean","Argentina_sd", "Turkey_mean","Turkey_sd", 
                                    "United Kingdom_mean","United Kingdom_sd", "Iran_mean", "Iran_sd",
                                    "Japan_mean", "Japan_sd","Poland_mean","Poland_sd", "Poland_1_mean","Poland_1_sd", "Kazakhstan_mean", "Kazakhstan_sd")

  # Render initial map
  output$map <- renderLeaflet({
    leaflet() %>%
      addTiles() %>%
      addMarkers(data = locations, ~Longitude, ~Latitude, popup = ~Location, layerId = ~Location, ) %>%
      setView(lng = locations$Longitude[1], lat = locations$Latitude[1], zoom = 2)
  })
  # Render additional plot when marker is clicked
  data_order <- c()
  output$additionalPlot <- renderPlot({
    req(input$map_marker_click)
    clicked_location <- locations[locations$Location == input$map_marker_click$id, ]
  
  all_data_new <- all_data_new[order(all_data_new[,paste(input$map_marker_click$id, "_mean", sep="")]),]
  all_data_new <- cbind(c(1:500),all_data_new)
  colnames(all_data_new)[1] <- c("date")
  print(colnames(all_data_new))
  print(c(1,paste(input$map_marker_click$id,"mean", sep="_"),paste(input$map_marker_click$id,"sd", sep="_")))
  for_plot <- all_data_new[,c("date",paste(input$map_marker_click$id,"mean", sep="_"),paste(input$map_marker_click$id,"sd", sep="_"))]
  print(colnames(for_plot))
  colnames(for_plot) <- c("date","mean", "sd")
  eb <- aes(ymax = mean + sd, ymin = mean - sd)
  for_plot <- as.data.frame(for_plot)
  ggplot(data = for_plot, aes(x = date, y = mean)) + 
    geom_line(size = 2) + 
    geom_ribbon(eb, alpha = 0.5) +
    ggtitle(paste("Isolate ",iso_name[input$map_marker_click$id,], sep="")) +
    labs(x = "accessions", y= "infection rate")
})
})





