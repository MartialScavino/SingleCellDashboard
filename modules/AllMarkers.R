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
  fileInput("importtablemarkers", label = "Enter a .csv file", accept = ".csv"),
  hr(),
  p(markdown("**Export marker table**")),
  downloadButton("exporttablemarkers", span("Export", id = "UpdateAnimateSaveMarkers", class = ""))
)



markers <- tabPanel("Markers", 
                    dataTableOutput("AllMarkersDataTable"))
heatmap <- tabPanel("Heatmap", 
                    plotOutput("HeatmapAllMarkers"))
dotplot <- tabPanel("Dot plot",
                    plotOutput("DotplotAllMarkers"))
volcano <- tabPanel("Volcano plot")





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
  
  
}
