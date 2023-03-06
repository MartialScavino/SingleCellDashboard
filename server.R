server <- function(input, output, session) {
  
  # Reactive value for dataset
  val <- reactiveValues(colors = list(), 
                        markers = data.frame(list()),
                        degs = data.frame())
  
  # Loading dataset
  observeEvent(input$data,{
    # Check for file type
    if (file_ext(input$data$datapath) != "rds"){ 
      alert("Enter a rds file")
      return(0)
    }
    
    df <- readRDS(input$data$datapath)

    if (!("percent.mt" %in% names(df@meta.data))){
      df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
    }  
    
    val$data <- df
    
    # Show Metadata
    output$dataset <- renderDataTable(val$data@meta.data, extensions = 'Buttons', 
                                      options = list(dom = 'Bfrtip', fixedColumns = TRUE,
                                                     buttons = c('copy', 'csv', 'excel')))
    })
 
  # Load data with 10X files
  observeEvent(input$filtered, {

    addClass(id = "UpdateAnimateLoad", class = "loading dots")
    disable("filtered")

    counts <- ReadMtx(input$matrixfile$datapath,
            cells = input$barcodesfile$datapath,
            features = input$featuresfile$datapath,
            feature.column = 1)

    df <- CreateSeuratObject(counts = counts)

    if (!("percent.mt" %in% names(df@meta.data)))
      df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")

    val$data <- df

    # Show Metadata
    output$dataset <- renderDataTable(val$data@meta.data, extensions = 'Buttons',
                                      options = list(dom = 'Bfrtip', fixedColumns = TRUE,
                                                     buttons = c('copy', 'csv', 'excel')))

    enable("filtered")
    removeClass(id = "UpdateAnimateLoad", class = "loading dots")
  })
  
  
  # QC
  ## Done in modules/QC.R
  ## Plots and inputs to trim data
  qcserver(input, output, session, val)
  
  
  # Preprocessing
  observeEvent(input$LaunchPreprocessing, {
    
    addClass(id = "UpdateAnimatePreprocessing", class = "loading dots")
    disable("LaunchPreprocessing")
    
    val$data <- NormalizeData(val$data, normalization.method = input$NormMethod, 
                              scale.factor = input$ScaleFactor)
    val$data <- FindVariableFeatures(val$data, selection.method = input$VariableMethod,
                                     nfeatures = input$nfeatures)
    
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    
    # Convert S genes
    s.genes_mus.musculus <- orthologs(genes = s.genes, species = "mouse")
    
    # Convert G2M genes
    g2m.genes_mus.musculus <- orthologs(genes = g2m.genes, species = "mouse")
    
    if (input$species == "Human")
      val$data <- CellCycleScoring(val$data, s.features = s.genes, g2m.features = g2m.genes)
    
    else
      val$data <- CellCycleScoring(val$data, s.features = s.genes_mus.musculus$symbol, g2m.features = g2m.genes_mus.musculus$symbol)
    
    
    dict <- list("All genes" = rownames(val$data), "Variable genes" = VariableFeatures(val$data))
    
    if (input$regressing)
      val$data <- ScaleData(val$data, features = as.vector(unlist(dict[input$FeatureScale])), 
                            vars.to.regress = c("S.Score", "G2M.Score"))

    else
      val$data <- ScaleData(val$data, features = as.vector(unlist(dict[input$FeatureScale])))
    
    
    enable("LaunchPreprocessing")
    removeClass(id = "UpdateAnimatePreprocessing", class = "loading dots")
    
  })
  
  
  output$variablefeatureplot <- renderPlot({
    
    if (length(VariableFeatures(val$data)) == 0)
      return(0)
    
    top20 <- head(VariableFeatures(val$data), 20)
    p1 <- VariableFeaturePlot(val$data)
    LabelPoints(plot = p1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)
    
    })
  
  
  # Dimension reduction
  ## PCA
  output$groupbypca <- renderUI(selectInput("groupbypca", "group cells by",
                                         choices = names(val$data@meta.data)))
  
  observeEvent(input$dopca, {
    
    addClass(id = "UpdateAnimatePCA", class = "loading dots")
    disable("dopca")
    
    dict <- list("All genes" = rownames(val$data), "Variable genes" = VariableFeatures(val$data))
    
    val$data <- RunPCA(val$data, features = as.vector(unlist(dict[input$featurePCA])))
    
    enable("dopca")
    removeClass(id = "UpdateAnimatePCA", class = "loading dots")
  })
  
  
  
  output$PCA <- renderPlot({
    if (!("pca" %in% names(val$data@reductions)))
      return(0)
    
    if(length(table(val$data@meta.data[,input$groupbypca])) > 50)
      FeaturePlot(val$data, reduction = "pca", features = input$groupbypca)
    
    else
      DimPlot(val$data, reduction = "pca", group.by = input$groupbypca) +
      scale_color_manual(values = as.vector(unlist(val$colors[input$groupbypca])))
    })
  
  output$elbow <- renderPlot({
    if (!("pca" %in% names(val$data@reductions)))
      return(0)
    ElbowPlot(val$data, ndims = input$ndimElbow)
    })
  
  output$dimheatmap <- renderPlot({
    if (!("pca" %in% names(val$data@reductions)))
      return(0)
    
    DimHeatmap(val$data, dims =  input$ndimheatmap, balanced = TRUE, nfeatures = 20)
    })
  
  
  
  ## UMAP
  output$groupbyumap <- renderUI(selectInput("groupbyumap", "group cells by",
                                             choices = names(val$data@meta.data)))
  
  observeEvent(input$doumap,{
    
    addClass(id = "UpdateAnimateUMAP", class = "loading dots")
    disable("doumap")
      
               val$data <- RunUMAP(val$data, dims = 1:input$ndimumap)
               
    enable("doumap")
    removeClass(id = "UpdateAnimateUMAP", class = "loading dots")
               })
  
  output$UMAP <- renderPlot({
    if(!("umap" %in% names(val$data@reductions)))
      return(0)
    
    if(length(table(val$data@meta.data[,input$groupbyumap])) > 50)
      FeaturePlot(val$data, reduction = "umap", features = input$groupbyumap)
    
    else
      DimPlot(val$data, reduction = "umap", group.by = input$groupbyumap) + 
      scale_color_manual(values = as.vector(unlist(val$colors[input$groupbyumap])))
    
  })
  
  
  ## Clustering
  output$groupbycluster <- renderUI(selectInput("groupbycluster", "group cells by",
                                             choices = names(val$data@meta.data),
                                             selected = "seurat_clusters"))
  
  output$sliderndimclusters <- renderUI(isolate(
    sliderInput("ndimscluster", "Number of dimension", 
                min = 1, 
                max = length(val$data@reductions$pca), 
                value = 30, 
                step = 1)
  ))
  
  observeEvent(input$docluster, {
    addClass(id = "UpdateAnimateCluster", class = "loading dots")
    disable("docluster")
    
    val$data <- FindNeighbors(val$data, k.param = input$kparam, dims = 1:input$ndimscluster)
    val$data <- FindClusters(val$data, resolution = input$resolution, algorithm = as.integer(input$algocluster))
    
    enable("docluster")
    removeClass(id = "UpdateAnimateCluster", class = "loading dots")
  })
    
  output$UMAPCluster <- renderPlot({
    if(length(table(val$data@meta.data[,input$groupbycluster])) > 50)
      FeaturePlot(val$data, reduction = "umap", features = input$groupbycluster) 
    
    else
      DimPlot(val$data, reduction = "umap", group.by = input$groupbycluster) + 
      scale_color_manual(values = val$colors[[input$groupbycluster]])
  })
  
  
  # Button to export data as rds
  output$savedata <- downloadHandler(
    filename = "seurat_object.rds",
    content = function(file) {
      addClass(id = "UpdateAnimateSave", class = "loading dots")
      disable("savedata")
      
      saveRDS(val$data, file)
      
      enable("savedata")
      removeClass(id = "UpdateAnimateSave", class = "loading dots")
    }
  )
  
  
  ## FeaturePlot
  # Done in modules/FeaturePlot.R
  featureplotserver(input, output, session, val)
  
  
  ## Options
  # Done in mdoules/Options.R
  optionsserver(input, output, session, val)
  
  
  ## Markers by ident
  # Done in modules/AllMarkers.R
  allmarkersserver(input, output, session, val)
  
  
  ## Markers between 2 idents
  # Done in modules/DifferentialExpression.R
  degserver(input, output, session, val)
  
  ## Enrichment analysis
  # Done in modules/Enrichment.R
  enrichmentserver(input, output, session, val)
  
  
  ## Make your own plot
  # Done in modules/Visualisation.R
  visualisationserver(input, output, session, val)
   
}