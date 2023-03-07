visualisationsidebar <- sidebarPanel(
  selectInput("plottype", "Plot", 
              choices = c("Cell projection", "Violin Plot", "Stacked Barchart")),
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
    notkeep <- c()
    for (i in names(val$data@meta.data)){
      if (length(table(val$data@meta.data[,i])) < 50)
        keep <- append(keep, i)
      else
        notkeep <- append(notkeep, i)}
    
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
             ),
           
           "Violin Plot" = tagList(
             selectInput("xviolinvisu", "X axis", keep, selected = input$xviolinvisu),
             selectInput("yviolinvisu", "Y axis", notkeep, selected = input$yviolinvisu)
           ),
           
           "Stacked Barchart" = tagList(
             selectInput("xbarvisu", "X axis", keep, selected = input$xbarvisu),
             selectInput("ybarvisu", "Y axis", keep, selected = input$ybarvisu),
             checkboxInput("capbar", "Normalize bars", input$capbar)
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
             },
           
           "Violin Plot" = {
             VlnPlot(val$data, features = input$yviolinvisu, group.by = input$xviolinvisu) + 
               scale_fill_manual(values = as.vector(unlist(val$colors[input$xviolinvisu])))
             },
           
           
           "Stacked Barchart" = {
             if (input$capbar)
               ggplot(val$data@meta.data, 
                      aes(fill=.data[[input$ybarvisu]], x= .data[[input$xbarvisu]])) + 
               geom_bar(position="fill") + 
               scale_fill_manual(values = as.vector(unlist(val$colors[input$ybarvisu]))) + 
               theme_light()
               
             else
               ggplot(val$data@meta.data, 
                      aes(fill=.data[[input$ybarvisu]], x= .data[[input$xbarvisu]])) + 
               geom_bar(position="stack") + 
               scale_fill_manual(values = as.vector(unlist(val$colors[input$ybarvisu]))) + 
               theme_light()
               
           }
           
           
           )
  })
  
}