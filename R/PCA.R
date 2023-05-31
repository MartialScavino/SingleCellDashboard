# Dimension reduction
pca <- tabItem(tabName = "pca",
               sidebarLayout(
                 sidebarPanel = sidebarPanel(
                   uiOutput("selectassayPCA"),
                   HelpInput("Features to compute PCA on. 
                             Note that the features must be present in the scaled data. 
                             Any requested features that are not scaled or have 
                             0 variance will be dropped, and the PCA will be run 
                             using the remaining features."),
                   selectInput("featurePCA", "Genes to use",
                               choices = c("Variable genes", "All genes"),
                               selected = "Variable genes"),
                   actionButton("dopca", span("Compute PCA", id="UpdateAnimatePCA", class=""), styleclass = "primary")
                 ),
                 mainPanel = mainPanel(
                   tabsetPanel(
                     tabPanel("Visualisation",
                              uiOutput("groupbypca"),
                              plotOutput("PCA")
                     ),
                     tabPanel("ElbowPlot",
                              sliderInput("ndimElbow", label = "Number of dimension to plot",
                                          10, 50, 30, 5),
                              plotOutput("elbow")),
                     
                     tabPanel("DimHeatmap",
                              sliderInput("ndimheatmap", "Number of dimension to plot",
                                          1, 50, 1, 1),
                              plotOutput("dimheatmap")),
                     
                     tabPanel("PCA Features",
                              uiOutput("selectpcafeatures"),
                              box(verbatimTextOutput("textpcafeaturespos"), title = "Positive"),
                              box(verbatimTextOutput("textpcafeaturesneg"), title = "Negative"))
                   )
                 )
               )
)


PCAServer <- function(input, output, session, val){
  
  ## PCA
  output$groupbypca <- renderUI({
    if (is.null(val$data))
      return("")
    
    selectInput("groupbypca", "group cells by",
                choices = names(val$data@meta.data))
  })
  
  output$selectassayPCA <- renderUI({
    if (is.null(val$data))
      return("")
    
    selectInput("assayPCA", "Assay to use", 
                names(val$data@assays))
  })
  
  
  observeEvent(input$dopca, {
    
    addClass(id = "UpdateAnimatePCA", class = "loading dots")
    disable("dopca")
    
    tryCatch({
      dict <- list("All genes" = rownames(val$data), "Variable genes" = VariableFeatures(val$data))
      val$data <- RunPCA(val$data, features = as.vector(unlist(dict[input$featurePCA])), assay = input$assayPCA)
      
    }, error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      enable("dopca")
      removeClass(id = "UpdateAnimatePCA", class = "loading dots")
      return(0)
    })
    
    enable("dopca")
    removeClass(id = "UpdateAnimatePCA", class = "loading dots")
  })
  
  
  
  output$PCA <- renderPlot({
    if (is.null(val$data))
      return(0)
    
    if (!("pca" %in% names(val$data@reductions)))
      return(0)
    
    if(length(table(val$data@meta.data[,input$groupbypca])) > 50)
      FeaturePlot(val$data, reduction = "pca", features = input$groupbypca)
    
    else
      DimPlot(val$data, reduction = "pca", group.by = input$groupbypca) +
      scale_color_manual(values = as.vector(unlist(val$colors[input$groupbypca])))
  })
  
  output$elbow <- renderPlot({
    if (is.null(val$data))
      return(0)
    
    if (!("pca" %in% names(val$data@reductions)))
      return(0)
    ElbowPlot(val$data, ndims = input$ndimElbow)
  })
  
  output$dimheatmap <- renderPlot({
    if (is.null(val$data))
      return(0)
    
    if (!("pca" %in% names(val$data@reductions)))
      return(0)
    
    DimHeatmap(val$data, dims =  input$ndimheatmap, balanced = TRUE)
  })
  
  
  output$selectpcafeatures <- renderUI({
    if (is.null(val$data))
      return("")
    
    selectInput("pcafeatures", "PC to show", 
                choices = 1:length(val$data@reductions$pca))
  })
  
  output$textpcafeaturespos <- renderText({
    if (is.null(val$data))
      return("")
    
    if (! "pca" %in% names(val$data@reductions))
      return("No pca found")
    
    gene_list <- TopFeatures(object = val$data[["pca"]], dim = as.integer(input$pcafeatures), nfeatures = 20, balanced = T)
    
    return(paste0(gene_list$positive, collapse = "\n"))
    
  })
  
  output$textpcafeaturesneg <- renderText({
    if (is.null(val$data))
      return("")
    
    if (! "pca" %in% names(val$data@reductions))
      return("No pca found")
    
    gene_list <- TopFeatures(object = val$data[["pca"]], dim = as.integer(input$pcafeatures), nfeatures = 20, balanced = T)
    
    return(paste0(gene_list$negative, collapse = "\n"))
    
  })
  
  
}