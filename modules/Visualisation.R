visualisationsidebar <- sidebarPanel(
  selectInput("plottype", "Plot", 
              choices = c("Cell projection", "Violin", "Stacked Barchart")),
  hr(),
  uiOutput("visualisation_listoption"),
  hr(),
  numericInput("width", "Plot width (px)", value = 400),
  numericInput("height", "Plot height (px)", value = 400),
  downloadButton("exportplot", span("Export plot", id = "UpdateAnimateExportPlot", class = ""))
)
visualisationmainpanel <- mainPanel(
  plotOutput("PlotVisualisation")
)




visualisationserver <- function(input, output, session, val){
  output$visualisation_listoption <- renderUI({
    keep <- c()
    for (i in names(val$data@meta.data)){
      if (length(table(val$data@meta.data[,i])) < 50)
        keep <- append(keep, i)}
    
    switch(input$plottype,
           
           "Cell projection" = tagList(
             selectInput("dimreduccellproj", "Plot to project cells", c("UMAP" = "umap", "PCA" = "pca")),
             selectInput("groupbycellproj", "Group cells by", keep),
             sliderInput("ptsizecellproj", "Point size", 0, 10, 1, 0.1),
             sliderInput("labelsizecellproj", "Label size", 0, 20, 5, 0.1),
             checkboxInput("showlegendcellproj", "Show legend", T),
             checkboxInput("labelclustercellproj", "Labels on cells", F),
             checkboxInput("repelcellproj", "Repel labels", F),
             checkboxInput("shufflecellproj", "Shuffle cells", F)
             )
           
           
      )
    
    
  })
  
  
  output$PlotVisualisation <- renderPlot({
    
    switch(input$plottype,
           
           "Cell projection" = {
             if (!input$showlegendcellproj)
              DimPlot(val$data, reduction = input$dimreduccellproj,
                      group.by = input$groupbycellproj, pt.size = input$ptsizecellproj,
                      label.size = input$labelsizecellproj, label = input$labelclustercellproj,
                      shuffle = input$shufflecellproj, repel = input$repelcellproj) +  NoLegend() + 
                      scale_color_manual(values = val$colors[[input$groupbycellproj]])
             else
               DimPlot(val$data, reduction = input$dimreduccellproj,
                       group.by = input$groupbycellproj, pt.size = input$ptsizecellproj,
                       label.size = input$labelsizecellproj, label = input$labelclustercellproj,
                       shuffle = input$shufflecellproj, repel = input$repelcellproj) + 
                       scale_color_manual(values = val$colors[[input$groupbycellproj]])
             }
           
           
           )
  })
  
}