library(shiny)
library(edgebundleR)

shinyUI(fluidPage(
  selectInput("variable", "Variable:",
              c("CAP45.G3" = "CAP45",
                "CAP45.G3M" = "CAP45_301",
                "Du156.12WT" = "Du156",
                "Du156.12M" = "Du156_301")),
    mainPanel(
      edgebundleOutput("distPlot", width = 800, height = 800))
  )
)