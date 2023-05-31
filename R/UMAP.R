umap <- tabItem(tabName = "umap", 
                sidebarLayout(
                  sidebarPanel = sidebarPanel(
                    uiOutput("selectassayUMAP"),
                    HelpInput("Which dimensions to use as input features. This is usually determined by examinating the elbow plot in
                              the PCA tab."),
                    numericInput("ndimumap", "Number of dimension",
                                 30, 1, 200),
                    actionButton("doumap", span("Compute UMAP", id = "UpdateAnimateUMAP", class=""), styleclass = "primary")
                  ),
                  
                  mainPanel = mainPanel(
                    uiOutput("groupbyumap"),
                    plotOutput("UMAP")
                  )
                ))


UMAPServer <- function(input, output, session, val){
  
  output$groupbyumap <- renderUI({
    if (is.null(val$data))
      return("")
    
    selectInput("groupbyumap", "group cells by",
                choices = names(val$data@meta.data))
  })
  
  output$selectassayUMAP <- renderUI({
    if (is.null(val$data))
      return("")
    
    selectInput("assayUMAP", "Assay to use", 
                names(val$data@assays))
  })
  
  observeEvent(input$doumap,{
    
    addClass(id = "UpdateAnimateUMAP", class = "loading dots")
    disable("doumap")
    
    tryCatch({
      val$data <- RunUMAP(val$data, dims = 1:input$ndimumap, assay = input$assayUMAP)
    }, error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      enable("doumap")
      removeClass(id = "UpdateAnimateUMAP", class = "loading dots")
      return(0)
    }) 
    enable("doumap")
    removeClass(id = "UpdateAnimateUMAP", class = "loading dots")
  })
  
  
  
  output$UMAP <- renderPlot({
    if (is.null(val$data))
      return(0)
    
    if(!("umap" %in% names(val$data@reductions)))
      return(0)
    
    if(length(table(val$data@meta.data[,input$groupbyumap])) > 50)
      FeaturePlot(val$data, reduction = "umap", features = input$groupbyumap)
    
    else
      DimPlot(val$data, reduction = "umap", group.by = input$groupbyumap) + 
      scale_color_manual(values = as.vector(unlist(val$colors[input$groupbyumap])))
    
  })
  
}