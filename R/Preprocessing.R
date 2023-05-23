preprocess <- tabItem(tabName = "preprocessing",
                      sidebarLayout(
                        sidebarPanel = sidebarPanel(
                          radioButtons("species", "Species", c("Human", "Mouse"),"Human"),
                          hr(),
                          HelpInput(markdown("
                          * LogNormalize: Feature counts for 
                          each cell are divided by the total counts for that cell
                          and multiplied by the scale.factor. 
                          This is then natural-log transformed using log1p.
                          
                          * CLR: Applies a centered log ratio transformation
                          
                          * RC: Relative counts. Feature counts for each cell are 
                                             divided by the total counts for that 
                                             cell and multiplied by the scale.factor. 
                                             No log-transformation is applied. 
                                             For counts per million (CPM) 
                                             set scale.factor = 1e6
                                             
                          In practice we almost never change this parameter so it is advised to keep the default value.")),
                          selectInput("NormMethod", "Normalization method",
                                      choices = c("LogNormalize" = "LogNormalize", 
                                                  "Centered log ratio" = "CLR",
                                                  "Relative counts" = "RC"), 
                                      selected = "LogNormalize"),
                          HelpInput("Sets the scale factor for cell-level normalization.
                                    
                                    Instead of normalizing all counts to 1, it normalize all counts to the scale factor.
                                    It makes expression values more readable but it doesn't change further analyses.
                                    
                                    It is advised to keep the default value"),
                          numericInput("ScaleFactor", "Scale factor",
                                       value = 10000,
                                       min = 1000, max = 1000000, step = 1000),
                          hr(),
                          HelpInput(markdown("How to choose top variable features. Choose one of :
                          
                          * vst: First, fits a line to the relationship of log(variance) and 
                          log(mean) using local polynomial regression (loess). 
                          Then standardizes the feature values using the observed 
                          mean and expected variance (given by the fitted line). 
                          Feature variance is then calculated on the standardized 
                          values after clipping to a maximum (see clip.max parameter).
                          
                          * dispersion (disp): selects the genes with the highest dispersion values")),
                          selectInput("VariableMethod", "Variable features selection method",
                                      choices = c("Vst" = "vst", "Dispersion" = "dispersion"),
                                      selected = "Vst"),
                          HelpInput("Number of features to select as top variable features"),
                          numericInput("nfeatures", "Number of variable features",
                                       value = 2000, 
                                       min = 1000, max = 10000, step = 500),
                          hr(),
                          # uiOutput("selectassayPreprocessing"),
                          HelpInput("Variables to regress out. Note that selecting all genes will be more time-consuming and may generate noise
                                    but can also allow a deeper analysis"),
                          selectInput("FeatureScale", "Set of genes to use for scaling",
                                      choices = c("Variable genes", "All genes")),
                          HelpInput("Regressing the cell cycle will be more time-consuming but it can be a good idea if you notice that your cells clusterize
                                    by cell cycle phase"),
                          checkboxInput("regressing", "Regress cell cycle"),
                          hr(),
                          actionButton("LaunchPreprocessing", span("Compute preprocessing", id="UpdateAnimatePreprocessing", class=""), styleclass = "primary")
                        ),
                        mainPanel = mainPanel(
                          
                          plotOutput("variablefeatureplot")
                        )
                      ))


PreprocessingServer <- function(input, output, session, val){
  
  
  
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
      
      
      dict <- list("All genes" = rownames(GetAssay(val$data), "RNA"), "Variable genes" = VariableFeatures(val$data))
      
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
  
  
  
}

