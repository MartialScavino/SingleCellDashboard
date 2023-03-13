sidebarDEG <- sidebarPanel(
  uiOutput("selectidentdeg"),
  uiOutput("selectident1deg"),
  uiOutput("selectident2deg"),
  hr(),
  p(markdown("#### Subset")),
  uiOutput("selectgroupbyDE"),
  uiOutput("selectsubsetident"),
  hr(),
  selectInput("TestToUseDEG", "Test to use", 
              choices = c("Wilcoxon" = "wilcox", "Likelihood Ratio" = "bimod", "ROC" = "roc",
                          "T test" = "t", "Negative Binomial" = "negbinom",
                          "Poisson" = "poisson", "Logistic Regression" = "LR",
                          "MAST" = "MAST", "DESeq2" = "DESeq2")),
  sliderInput("MinimumPercentDEG", "Minimum percent expressed",
              0, 1, 0.1, 0.01),
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
    keep <- c()
    for (i in names(val$data@meta.data)){
      if (length(table(val$data@meta.data[,i])) < 50)
        keep <- append(keep, i)
    }
    selectInput("choiceidentdeg", "Select a variable", choices = keep, selected = input$choiceidentdeg)
  })
  
  output$selectident1deg <- renderUI(selectInput("ident1deg", "First group to compare", 
                                                 choices = names(table(val$data@meta.data[,input$choiceidentdeg])), 
                                     selected = input$ident1deg))
  output$selectident2deg <- renderUI(selectInput("ident2deg", "Second group to compare", 
                                                 choices = c("All others", names(table(val$data@meta.data[,input$choiceidentdeg]))),
                                                 selected = input$ident2deg))
  
  output$selectgroupbyDE <- renderUI({
    keep <- c()
    for (i in names(val$data@meta.data)){
      if (length(table(val$data@meta.data[,i])) < 50)
        keep <- append(keep, i)
    }
    selectInput("choiceidentsubset", "Select a group to subset", choices = c("None", keep), selected = input$choiceidentsubset)
  })
  
  output$selectsubsetident <- renderUI({
    if (input$choiceidentsubset != 'None')
      selectInput("subsetident", "Choose the cells to keep", choices = names(table(val$data@meta.data[,input$choiceidentsubset])))
    
  })
  
  
  observeEvent(input$dofindDEG,{
    
    addClass(id = "UpdateAnimateDEG", class = "loading dots")
    disable("dofindDEG")
    
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
    
    enable("dofindDEG")
    removeClass(id = "UpdateAnimateDEG", class = "loading dots")
    
  })
  
  output$DEGs <- renderDataTable(val$degs, extensions = 'Buttons', server = F,
                                                options = list(dom = 'Bfrtip', fixedColumns = TRUE,
                                                               buttons = c('csv', 'excel')))
  
  output$copygenes <- renderUI({
    if (length(val$degs) != 0)
      rclipButton("copybtm", "Copy all genes", 
                  paste(rownames(val$degs), collapse = "\n"), modal = T, styleclass = "primary", icon = icon("clipboard"))
    
  })
  
  
  output$VolcanoDEGS <- renderPlotly({
    
    up <- val$degs[which(val$degs$avg_log2FC > input$LogFCThreshold),]
    down <- val$degs[which(val$degs$avg_log2FC < -input$LogFCThreshold),]
    
    top10up <- up %>% slice_max(n = 10, order_by = avg_log2FC)
    top10down <- down %>% slice_max(n = 10, order_by = -avg_log2FC)
    
    ggplotly(
      ggplot(val$degs, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(data = up, color = "red") + 
        geom_point(data = down, color = "blue") + 
        geom_text(data = top10up, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene), nudge_y = 5) +
        geom_text(data = top10down, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene), nudge_y = 5) + theme_light()
      
    )
    
  })
  
}

#test <- readRDS("../../Downloads/test_tissue.rds")
#head(test)
#Idents(test)
## je veux séparer les cellules du cluster 1 entre ADK et AT
#FindMarkers(test, ident.1 = "ADK", ident.2 = "AT", subset.ident = 1, group.by = "tissue")
