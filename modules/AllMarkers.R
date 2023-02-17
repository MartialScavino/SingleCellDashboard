sidebarAllMarkers <- sidebarPanel(
  uiOutput("SelectIdentAllMarkers"),
  selectInput("TestToUse", "Test to use", 
              choices = c("Wilcoxon" = "wilcox", "Likelihood Ratio" = "bimod", "ROC" = "roc",
                          "T test" = "t", "Negative Binomial" = "negbinom",
                          "Poisson" = "poisson", "Logistic Regression" = "LR",
                          "MAST" = "MAST", "DESeq2" = "DESeq2")),
  sliderInput("MinimumPercent", "Minimum percent expressed",
              0, 1, 0.1, 0.01),
  sliderInput("LogFCThreshold", "LogFC threshold", 
              0, 3, 0.25, 0.05),
  sliderInput("PValueThreshold", "Adjusted p-value threshold",
              0, 1, 0.01, 0.01),
  checkboxInput("OnlyPos", "Only upregulated genes", value = T),
  actionButton("dofindallmarkers", span("Compute", id = "UpdateAnimateAllMarkers", class=""), styleclass = "primary"),
  hr(),
  p(markdown("**Import marker table**")),
  fileInput("importtablemarkers", label = "Enter a .csv or a excel file", accept = ".csv")
)



markers <- tabPanel("Markers", 
                    dataTableOutput("AllMarkersDataTable"))
heatmap <- tabPanel("Heatmap", 
                    plotOutput("HeatmapAllMarkers"))
dotplot <- tabPanel("Dot plot",
                    plotOutput("DotplotAllMarkers"))
volcano <- tabPanel("Volcano plot",
                    uiOutput("selectVolcanoPlotAllMarkers"),
                    plotlyOutput("VolcanoPlotAllMarkers"))





allmarkersserver <- function(input, output, session, val){
  
  output$SelectIdentAllMarkers <- renderUI({
    keep <- c()
    for (i in names(val$data@meta.data)){
      if (length(table(val$data@meta.data[,i])) < 50)
        keep <- append(keep, i)
    }
    selectInput("identallmarkers", "Select Ident", choices = keep, selected = input$identallmarkers)
    
  })
  
  observeEvent(input$dofindallmarkers, {
    addClass(id = "UpdateAnimateAllMarkers", class = "loading dots")
    disable("dofindallmarkers")
    
    Idents(val$data) <- input$identallmarkers
    val$markers <- FindAllMarkers(val$data, logfc.threshold = input$LogFCThreshold,
                                  min.pct = input$MinimumPercent, test.use = input$TestToUse,
                                  only.pos = input$OnlyPos)
    
    val$markers <- val$markers[which(val$markers$p_val_adj < input$PValueThreshold),]
    
    enable("dofindallmarkers")
    removeClass(id = "UpdateAnimateAllMarkers", class = "loading dots")
    
  })
  
  
  output$AllMarkersDataTable <- renderDataTable(val$markers, extensions = 'Buttons', 
                                                                  options = list(dom = 'Bfrtip', fixedColumns = TRUE,
                                                                                 buttons = c('copy', 'csv', 'excel')))
  
  
  output$HeatmapAllMarkers <- renderPlot({
    
    top10 <- val$markers %>% 
      group_by(cluster) %>% 
      slice_max(order_by = avg_log2FC, n = 10)
    
    top10 <- top10 %>% 
      arrange(factor(cluster, levels = levels(val$data@meta.data[,input$identallmarkers])))
    
    DoHeatmap(val$data, features = top10$gene, 
              group.colors = as.vector(unlist(val$colors[input$identallmarkers])),
              group.by = input$identallmarkers)
  }, height = 500)
  
  
  
  output$DotplotAllMarkers <- renderPlot({
    
    top10 <- val$markers %>% 
      group_by(cluster) %>% 
      slice_max(order_by = avg_log2FC, n = 10)
    
    top10 <- top10 %>% 
      arrange(factor(cluster, levels = levels(val$data@meta.data[,input$identallmarkers])))
    
    DotPlot(val$data, features = unique(top10$gene), group.by = input$identallmarkers) +
       theme(axis.text.x = element_text(angle = 60, size = 8, vjust = 0.85))
    
  })
  
  # Select Input pour le volcano plot
  output$selectVolcanoPlotAllMarkers <- renderUI(
    selectInput("choiceVolcanoPlotAllMarkers", "Choose a variable",
                choices = names(table(val$data@meta.data[,input$identallmarkers])))
  )
  
  output$VolcanoPlotAllMarkers <- renderPlotly({
    
    subdf <- val$markers[which(val$markers$cluster == input$choiceVolcanoPlotAllMarkers),]
    subdf_up <- subdf[which(subdf$avg_log2FC > input$LogFCThreshold),]
    subdf_down <- subdf[which(subdf$avg_log2FC < -input$LogFCThreshold),]
    # subdf_notde <- subdf[which(subdf$avg_log2FC > -input$LogFCThreshold & 
    #                              subdf$avg_log2FC < input$LogFCThreshold),]
    
    top10up <- subdf_up %>% slice_max(n = 10, order_by = avg_log2FC)
    top10down <- subdf_down %>% slice_max(n = 10, order_by = -avg_log2FC)
    
    ggplotly(
      ggplot(subdf, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(data = subdf_up, color = "red") + 
        geom_point(data = subdf_down, color = "blue") + 
        geom_text(data = top10up, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene), nudge_y = 5) +
        geom_text(data = top10down, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene), nudge_y = 5) + theme_light()
      
    )
    
  })
  
  
  observeEvent(input$importtablemarkers, {
    
    if (file_ext(input$importtablemarkers$datapath) == "xlsx")
      val$markers <- readxl::read_xlsx(input$importtablemarkers$datapath, col_names = T)
    
    else if (file_ext(input$importtablemarkers$datapath) == "csv")
      val$markers <- read.csv(input$importtablemarkers$datapath, header = T, sep = ',', row.names = 1)
    
    else{
      alert("Please enter a csv or a xlsx file")
      return(0)
    }
      
  })
  

}


# testdf <- FindMarkers(test, ident.1 = 1, group.by = "seurat_clusters", logfc.threshold = 2, min.diff.pct = 0.1)
# testdf$avg_log2FC
# testdf$pct.2 - testdf$pct.1
