sidebarDEG <- sidebarPanel(
  uiOutput("selectidentdeg"),
  uiOutput("selectident1deg"),
  uiOutput("selectident2deg"),
  hr(),
  p(markdown("#### Subset")),
  HelpInput("This part lets you test DEGs between 2 conditions within a subset of cell. 
            For example, if you wish to find the DEGs between the tissu 1 and 2 within a celltype A, 
            just select the tissu column with the above buttons and select the celltype column with the following buttons"),
  uiOutput("selectgroupbyDE"),
  uiOutput("selectsubsetident"),
  hr(),
  HelpInput(text_placement = "right", trigger = "click", content = markdown("Denotes which test to use. Available options are:
  
  * wilcox : Identifies differentially expressed genes between two groups of cells
  using a Wilcoxon Rank Sum test (default)
  
  * bimod : Likelihood-ratio test for single cell gene expression, 
  (McDavid et al., Bioinformatics, 2013)
  
  * roc : Identifies 'markers' of gene expression using ROC analysis. 
  For each gene, evaluates (using AUC) a classifier built on that gene alone, 
  to classify between two groups of cells. An AUC value of 1 means that expression
  values for this gene alone can perfectly classify the two groupings 
  (i.e. Each of the cells in cells.1 exhibit a higher level than each of the cells in cells.2).
  An AUC value of 0 also means there is perfect classification, but in the other direction.
  A value of 0.5 implies that the gene has no predictive power to classify the two groups. 
  Returns a 'predictive power' abs(AUC-0.5)*2 ranked matrix of putative differentially expressed genes.
  
  * t : Identify differentially expressed genes between two groups of cells 
  using the Student's t-test.
  
  * negbinom : Identifies differentially expressed genes between two groups of cells 
  using a negative binomial generalized linear model. Use only for UMI-based datasets
  
  * poisson : Identifies differentially expressed genes between two groups of cells 
  using a poisson generalized linear model. Use only for UMI-based datasets
  
  * LR : Uses a logistic regression framework to determine differentially expressed genes. 
  Constructs a logistic regression model predicting group membership based on each feature individually and 
  compares this to a null model with a likelihood ratio test.
  
  * MAST : Identifies differentially expressed genes between two groups of cells 
  using a hurdle model tailored to scRNA-seq data. Utilizes the MAST package to run the DE testing.
  
  * DESeq2 : Identifies differentially expressed genes between two groups of cells 
                     based on a model using DESeq2 which uses a negative binomial distribution 
                     (Love et al, Genome Biology, 2014).This test does not support pre-filtering of 
                     genes based on average difference (or percent detection rate) between cell groups. 
                     However, genes may be pre-filtered based on their minimum detection rate (min.pct) across 
                     both cell groups. To use this method, please install DESeq2, using the instructions at 
                     https://bioconductor.org/packages/release/bioc/html/DESeq2.html")),
  selectInput("TestToUseDEG", "Test to use", 
              choices = c("Wilcoxon" = "wilcox", "Likelihood Ratio" = "bimod", "ROC" = "roc",
                          "T test" = "t", "Negative Binomial" = "negbinom",
                          "Poisson" = "poisson", "Logistic Regression" = "LR",
                          "MAST" = "MAST", "DESeq2" = "DESeq2")),
  HelpInput("only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. 
            Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1"),
  sliderInput("MinimumPercentDEG", "Minimum percent expressed",
              0, 1, 0.1, 0.01),
  HelpInput("Limit testing to genes which show, on average, at least X-fold difference (log-scale) 
            between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, 
            but can miss weaker signals."),
  sliderInput("LogFCThresholdDEG", "LogFC threshold", 
              0, 3, 0.25, 0.05),
  sliderInput("PValueThresholdDEG", "Adjusted p-value threshold",
              0, 1, 0.01, 0.01),
  checkboxInput("OnlyPosDEG", "Only upregulated genes", value = T),
  actionButton("dofindDEG", span("Compute", id = "UpdateAnimateDEG", class=""), styleclass = "primary")
)


degtable <- tabPanel("DEGs",
                     uiOutput("copygenes"),
                     dataTableOutput("DEGs"))

volcanodeg <- tabPanel("Volcano plot",
                       plotlyOutput("VolcanoDEGS"))



degserver <- function(input, output, session, val){
  
  output$selectidentdeg <- renderUI({
    if (is.null(val$data))
      return("")
    
    keep <- c()
    for (i in names(val$data@meta.data)){
      if (length(table(val$data@meta.data[,i])) < 50)
        keep <- append(keep, i)
    }
    selectInput("choiceidentdeg", "Select a variable", choices = keep, selected = input$choiceidentdeg)
  })
  
  output$selectident1deg <- renderUI({
    if (is.null(val$data))
      return("")
    
    selectInput("ident1deg", "First group to compare", 
                choices = names(table(val$data@meta.data[,input$choiceidentdeg])), 
                selected = input$ident1deg)
    })
  output$selectident2deg <- renderUI({
    if (is.null(val$data))
      return("")
    
    selectInput("ident2deg", "Second group to compare", 
                choices = c("All others", names(table(val$data@meta.data[,input$choiceidentdeg]))),
                selected = input$ident2deg)})
  
  output$selectgroupbyDE <- renderUI({
    if (is.null(val$data))
      return("")
    
    keep <- c()
    for (i in names(val$data@meta.data)){
      if (length(table(val$data@meta.data[,i])) < 50)
        keep <- append(keep, i)
    }
    selectInput("choiceidentsubset", "Select a group to subset", choices = c("None", keep), selected = input$choiceidentsubset)
  })
  
  output$selectsubsetident <- renderUI({
    if (is.null(val$data))
      return("")
    
    if (input$choiceidentsubset != 'None')
      selectInput("subsetident", "Choose the cells to keep", 
                  choices = names(table(val$data@meta.data[,input$choiceidentsubset])),
                  selected = input$subsetident)
    
  })
  
  
  observeEvent(input$dofindDEG,{
    
    addClass(id = "UpdateAnimateDEG", class = "loading dots")
    disable("dofindDEG")
    
    tryCatch({
    Idents(val$data) <- input$choiceidentdeg
    
    if (input$ident2deg == "All others"){
      if (input$choiceidentsubset == "None"){
      val$degs <- FindMarkers(object = val$data, logfc.threshold = input$LogFCThresholdDEG,
                                    min.pct = input$MinimumPercentDEG, test.use = input$TestToUseDEG,
                                    only.pos = input$OnlyPosDEG, ident.1 = input$ident1deg, assay = "RNA")
      }
      
      else{
        Idents(val$data) <- input$choiceidentsubset
        val$degs <- FindMarkers(val$data, logfc.threshold = input$LogFCThresholdDEG,
                                min.pct = input$MinimumPercentDEG, test.use = input$TestToUseDEG,
                                only.pos = input$OnlyPosDEG, subset.ident = input$subsetident,
                                group.by = input$choiceidentdeg, ident.1 = input$ident1deg, assay = "RNA")
      }
      
    }
    
    else{
      if (input$choiceidentsubset == "None"){
      val$degs <- FindMarkers(val$data, logfc.threshold = input$LogFCThresholdDEG,
                                 min.pct = input$MinimumPercentDEG, test.use = input$TestToUseDEG,
                                 only.pos = input$OnlyPosDEG, ident.1 = input$ident1deg, ident.2 = input$indent2deg, assay = "RNA")
      }
      
      else{
        
        Idents(val$data) <- input$choiceidentsubset
        val$degs <- FindMarkers(val$data, logfc.threshold = input$LogFCThresholdDEG,
                                min.pct = input$MinimumPercentDEG, test.use = input$TestToUseDEG,
                                only.pos = input$OnlyPosDEG, subset.ident = input$subsetident,
                                group.by = input$choiceidentdeg, ident.1 = input$ident1deg, ident.2 = input$indent2deg, assay = "RNA")
        
      }
      
    }
    
    val$degs <- val$degs[which(val$degs$p_val_adj < input$PValueThreshold),]
    val$degs[["gene"]] <- rownames(val$degs)
    
    }, error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      enable("dofindDEG")
      removeClass(id = "UpdateAnimateDEG", class = "loading dots")
      return(0)
    })
    
    enable("dofindDEG")
    removeClass(id = "UpdateAnimateDEG", class = "loading dots")
    
  })
  
  output$DEGs <- renderDataTable(val$degs, extensions = 'Buttons', server = F,
                                                options = list(dom = 'Bfrtip', fixedColumns = TRUE, scrollX = TRUE,
                                                               buttons = c('csv', 'excel')))
  
  output$copygenes <- renderUI({
    if (is.null(val$data))
      return("")
    
    if (length(val$degs) != 0)
      rclipButton("copybtm", "Copy all genes", 
                  paste(rownames(val$degs), collapse = "\n"), modal = T, styleclass = "primary", icon = icon("clipboard"))
    
  })
  
  
  output$VolcanoDEGS <- renderPlotly({
    if (is.null(val$data))
      return(ggplotly(ggplot() + theme_minimal()))
    
    up <- val$degs[which(val$degs$avg_log2FC > input$LogFCThreshold),]
    down <- val$degs[which(val$degs$avg_log2FC < -input$LogFCThreshold),]
    
    top10up <- up %>% slice_max(n = 10, order_by = avg_log2FC)
    top10down <- down %>% slice_max(n = 10, order_by = -avg_log2FC)
    
    ggplotly(
      ggplot(val$degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(data = up, color = "red") + 
        geom_point(data = down, color = "blue") + 
        geom_text(data = top10up, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene), nudge_y = 0.05) +
        geom_text(data = top10down, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene), nudge_y = 0.05) + theme_light()
      
    )
    
  })
  
}

#test <- readRDS("../../Downloads/test_tissue.rds")
#head(test)
#Idents(test)
## je veux séparer les cellules du cluster 1 entre ADK et AT
#FindMarkers(test, ident.1 = "ADK", ident.2 = "AT", subset.ident = 1, group.by = "tissue")
