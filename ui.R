library(leaflet)
shinyUI(fluidPage(
  leafletOutput("map"),
  plotOutput("additionalPlot")
))
