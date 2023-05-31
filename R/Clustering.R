clustering <- tabItem(tabName = "clustering",
                      sidebarLayout(
                        sidebarPanel = sidebarPanel(
                          uiOutput("selectassayCluster"),
                          numericInput("kparam", "k parameter", 
                                       20, 1, 100),
                          HelpInput("Dimensions of PCA to use as input. 
                                    Please select the same value you did to compute the UMAP"),
                          numericInput("ndimscluster", "Number of dimension", 
                                       min = 1, 
                                       max = 200, 
                                       value = 30, 
                                       step = 1),
                          HelpInput("Value of the resolution parameter, use a value above (below) 1.0 
                                    if you want to obtain a larger (smaller) number of communities."),
                          numericInput("resolution", "Resolution", 
                                       0.8, 0, 5, 0.1),
                          selectInput("algocluster", "Algorithm", 
                                      choices = c("Louvain" = "1", 
                                                  "Louvain with multilevel refinement" = "2",
                                                  "SLM" = "3",
                                                  "Leiden" = "4")),
                          actionButton("docluster", span("Compute clustering", id = "UpdateAnimateCluster", class = ""), styleclass = "primary")
                        ),
                        mainPanel = mainPanel(
                          uiOutput("groupbycluster"),
                          plotOutput("UMAPCluster")
                        )
                      ))


ClusteringServer <- function(input, output, session, val){
  output$groupbycluster <- renderUI({
    if (is.null(val$data))
      return("")
    
    selectInput("groupbycluster", "group cells by",
                choices = names(val$data@meta.data),
                selected = "seurat_clusters")
  })
  
  output$selectassayCluster <- renderUI({
    if (is.null(val$data))
      return("")
    
    selectInput("assayCluster", "Assay to use", 
                names(val$data@assays))
  })
  
  output$sliderndimclusters <- renderUI({
    if (is.null(val$data))
      return("")
      numericInput("ndimscluster", "Number of dimension", 
                  min = 1, 
                  max = length(val$data@reductions$pca), 
                  value = 30, 
                  step = 1)
    })
  
  observeEvent(input$docluster, {
    if (length(val$data) == 0){
      alert("No data loaded")
      return(0)
    }
    
    if (input$ndimscluster > length(val$data@reductions$pca)){
      alert(paste("The number of dimensions can't exceed", length(val$data@reductions$pca)))
      return(0)
    }
      
    
    
    addClass(id = "UpdateAnimateCluster", class = "loading dots")
    disable("docluster")
    
    tryCatch({
      
      DefaultAssay(val$data) <- input$assayCluster
      val$data <- FindNeighbors(val$data, k.param = input$kparam, dims = 1:input$ndimscluster)
      val$data <- FindClusters(val$data, resolution = input$resolution, algorithm = as.integer(input$algocluster))
    }, error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      DefaultAssay(val$data) <- "RNA"
      enable("docluster")
      removeClass(id = "UpdateAnimateCluster", class = "loading dots")
      return(0)
    })
    
    
    DefaultAssay(val$data) <- "RNA"
    enable("docluster")
    removeClass(id = "UpdateAnimateCluster", class = "loading dots")
  })
  
  output$UMAPCluster <- renderPlot({
    if (is.null(val$data))
      return(0)
    
    if(!("umap" %in% names(val$data@reductions)))
      return(0)
    
    if(length(table(val$data@meta.data[,input$groupbycluster])) > 50)
      FeaturePlot(val$data, reduction = "umap", features = input$groupbycluster) 
    
    else
      DimPlot(val$data, reduction = "umap", group.by = input$groupbycluster) + 
      scale_color_manual(values = val$colors[[input$groupbycluster]])
  })
  
}