server <- function(input, output, session) {
  
  # Reactive value for dataset
  val <- reactiveValues(colors = list(), 
                        markers = data.frame(list()),
                        degs = data.frame(),
                        enrichment = data.frame(),
                        plot = character(0))
  
  # Loading dataset
  observeEvent(input$data,{
    # Check for file type
    if (file_ext(input$data$datapath) != "rds"){ 
      alert("Enter a rds file")
      return(0)
    }
    
    tryCatch({
    df <- readRDS(input$data$datapath)
    
    if (!("percent.mt" %in% names(df@meta.data))){
      
      if (length(grep( "^mt-", rownames(df), value = T) > 0))
        df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^mt-")
      
      else if (length(grep( "^MT-", rownames(df), value = T) > 0))
        df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
    }  
    
    val$data <- df
    
    }, error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      return(0)
    }) 
    
    # Show Metadata
    output$dataset <- renderDataTable(val$data@meta.data, 
                                      options = list(dom = 'Bfrtip', 
                                                     fixedColumns = TRUE, 
                                                     scrollX = T))
    })
 
  # Load data with 10X files
  observeEvent(input$filtered, {

    addClass(id = "UpdateAnimateLoad", class = "loading dots")
    disable("filtered")

    tryCatch({
      
    counts <- ReadMtx(input$matrixfile$datapath,
            cells = input$barcodesfile$datapath,
            features = input$featuresfile$datapath,
            feature.column = 1)

    df <- CreateSeuratObject(counts = counts)

    if (!("percent.mt" %in% names(df@meta.data))){
      
      if (length(grep( "^mt-", rownames(df), value = T) > 0))
        df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^mt-")
      
      else if (length(grep( "^MT-", rownames(df), value = T) > 0))
        df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
    }

    val$data <- df
    
    }, error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      enable("filtered")
      removeClass(id = "UpdateAnimateLoad", class = "loading dots")
      return(0)
    })

    # Show Metadata
    output$dataset <- renderDataTable(val$data@meta.data, 
                                      options = list(dom = 'frtip', 
                                                     fixedColumns = TRUE, 
                                                     scrollX = T))

    enable("filtered")
    removeClass(id = "UpdateAnimateLoad", class = "loading dots")
  })
  
  
  # QC
  ## Done in modules/QC.R
  ## Plots and inputs to trim data
  qcserver(input, output, session, val)
  
  
  output$selectassayPreprocessing <- renderUI(selectInput("assayPreprocessing", "Assay to use", names(val$data@assays)))
  
  # Preprocessing
  observeEvent(input$LaunchPreprocessing, {
    
    addClass(id = "UpdateAnimatePreprocessing", class = "loading dots")
    disable("LaunchPreprocessing")
  
  tryCatch({
    
    val$data <- NormalizeData(val$data, normalization.method = input$NormMethod, 
                              scale.factor = input$ScaleFactor, assay = "RNA")
    val$data <- FindVariableFeatures(val$data, selection.method = input$VariableMethod,
                                     nfeatures = input$nfeatures, assay = "RNA")
    
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
                            vars.to.regress = c("S.Score", "G2M.Score"), assay = "RNA")

    else
      val$data <- ScaleData(val$data, features = as.vector(unlist(dict[input$FeatureScale])), assay = "RNA")
  }, error = function(e){
    alert("There has been an error (printed in R console)")
    print(e)
    enable("LaunchPreprocessing")
    removeClass(id = "UpdateAnimatePreprocessing", class = "loading dots")
    return(0)
  })
    
    enable("LaunchPreprocessing")
    removeClass(id = "UpdateAnimatePreprocessing", class = "loading dots")
    
  })
  
  
  output$variablefeatureplot <- renderPlot({
    if (is.null(val$data))
      return(0)
    
    if (length(VariableFeatures(val$data)) == 0)
      return(0)
    
    top20 <- head(VariableFeatures(val$data), 20)
    tryCatch({
    p1 <- VariableFeaturePlot(val$data, assay = "RNA")
    LabelPoints(plot = p1, points = top20, repel = TRUE, xnudge = 0, ynudge = 0)
    }, error = function(e){
      
      if (grepl("'FindVariableFeatures'", e) & "integrated" %in% names(val$data@assays))
        return(ggplot() + 
               ggtitle("Can't show plot because the variable feature analysis was performed on the integrated assay") + 
               theme_minimal())
      else{
        alert("There has been an error (printed in R console)")
        print(e)
        return(0)
      }
    })
    
    })
  
  
  # Dimension reduction
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
    dict <- list("All genes" = rownames(val$data), "Variable genes" = VariableFeatures(val$data))
    
    tryCatch({
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
    
    gene_list <- TopFeatures(object = val$data[["pca"]], dim = as.integer(input$pcafeatures), nfeatures = 20, balanced = T)
    
    return(paste0(gene_list$positive, collapse = "\n"))
    
    })
  
  output$textpcafeaturesneg <- renderText({
    if (is.null(val$data))
      return("")
    
    gene_list <- TopFeatures(object = val$data[["pca"]], dim = as.integer(input$pcafeatures), nfeatures = 20, balanced = T)
    
    return(paste0(gene_list$negative, collapse = "\n"))
    
  })
  
  
  ## UMAP
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
  
  
  ## Clustering
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
    
    isolate(
      sliderInput("ndimscluster", "Number of dimension", 
                min = 1, 
                max = length(val$data@reductions$pca), 
                value = 30, 
                step = 1)
  )})
  
  observeEvent(input$docluster, {
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