library(shiny)
library(shinyIncubator)



shinyUI(pageWithSidebar(
  headerPanel("NINJA: Nematode INdicator Joint Analysis (beta)","NINJA: Nematode INdicator Joint Analysis"),
  sidebarPanel(
    htmlOutput("upload"),
    htmlOutput("adjust"),
    htmlOutput("selectUI"),
    htmlOutput("mailto")
  ),
  mainPanel(
    htmlOutput("disclaimer"),
    tags$link(rel = 'stylesheet', type = 'text/css', href = 'styles.css'),
    progressInit(),
    htmlOutput("tabs")
  )
))


