sidebarAllMarkers <- sidebarPanel(
  uiOutput("SelectIdentAllMarkers"),
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
  selectInput("TestToUse", "Test to use", 
              choices = c("Wilcoxon" = "wilcox", "Likelihood Ratio" = "bimod", "ROC" = "roc",
                          "T test" = "t", "Negative Binomial" = "negbinom",
                          "Poisson" = "poisson", "Logistic Regression" = "LR",
                          "MAST" = "MAST", "DESeq2" = "DESeq2")),
  HelpInput("only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. 
            Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1"),
  sliderInput("MinimumPercent", "Minimum percent expressed",
              0, 1, 0.1, 0.01),
  HelpInput("Limit testing to genes which show, on average, at least X-fold difference (log-scale) 
            between the two groups of cells. Default is 0.25 Increasing logfc.threshold speeds up the function, 
            but can miss weaker signals."),
  sliderInput("LogFCThreshold", "LogFC threshold", 
              0, 3, 0.25, 0.05),
  sliderInput("PValueThreshold", "Adjusted p-value threshold",
              0, 1, 0.01, 0.01),
  checkboxInput("OnlyPos", "Only upregulated genes", value = T),
  actionButton("dofindallmarkers", span("Compute", id = "UpdateAnimateAllMarkers", class=""), styleclass = "primary"),
  hr(),
  p(markdown("**Import marker table**")),
  HelpInput("Computing marker genes can take some time and could need a huge amount of resources. 
            If you've already done this analysis, you can download the table with the csv button and reimport it with the following button."),
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
    if (is.null(val$data))
      return("")
    
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
    
    tryCatch({
    Idents(val$data) <- input$identallmarkers
    val$markers <- FindAllMarkers(val$data, logfc.threshold = input$LogFCThreshold,
                                  min.pct = input$MinimumPercent, test.use = input$TestToUse,
                                  only.pos = input$OnlyPos, assay = "RNA")
    
    val$markers <- val$markers[which(val$markers$p_val_adj < input$PValueThreshold),]
    
    }, error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      enable("dofindallmarkers")
      removeClass(id = "UpdateAnimateAllMarkers", class = "loading dots")
      return(0)
      
    })
    
    enable("dofindallmarkers")
    removeClass(id = "UpdateAnimateAllMarkers", class = "loading dots")
    
  })
  
  
  output$AllMarkersDataTable <- renderDataTable(val$markers, extensions = 'Buttons', server = F,
                                                                  options = list(dom = 'Bfrtip', fixedColumns = TRUE,
                                                                                 buttons = c('copy', 'csv', 'excel'), scrollX = T))
  
  
  output$HeatmapAllMarkers <- renderPlot({
    if (is.null(val$data))
      return(0)
    
    top10 <- val$markers %>% 
      group_by(cluster) %>% 
      slice_max(order_by = avg_log2FC, n = 10)
    
    top10 <- top10 %>% 
      arrange(factor(cluster, levels = levels(val$data@meta.data[,input$identallmarkers])))
    
    DoHeatmap(val$data, features = top10$gene, 
              group.colors = as.vector(unlist(val$colors[input$identallmarkers])),
              group.by = input$identallmarkers, )
  }, height = 800)
  
  
  
  output$DotplotAllMarkers <- renderPlot({
    if (is.null(val$data))
      return(0)
    
    top10 <- val$markers %>% 
      group_by(cluster) %>% 
      slice_max(order_by = avg_log2FC, n = 10)
    
    top10 <- top10 %>% 
      arrange(factor(cluster, levels = levels(val$data@meta.data[,input$identallmarkers])))
    
    DotPlot(val$data, features = unique(top10$gene), group.by = input$identallmarkers) +
       theme(axis.text.x = element_text(angle = 60, size = 6, vjust = 0.85))
    
  })
  
  # Select Input pour le volcano plot
  output$selectVolcanoPlotAllMarkers <- renderUI({
    if (is.null(val$data))
      return("")
    
    selectInput("choiceVolcanoPlotAllMarkers", "Choose a variable",
                choices = names(table(val$data@meta.data[,input$identallmarkers])))
  })
  
  output$VolcanoPlotAllMarkers <- renderPlotly({
    if (is.null(val$data))
      return(ggplotly(ggplot() + theme_minimal()))
    
    subdf <- val$markers[which(val$markers$cluster == input$choiceVolcanoPlotAllMarkers),]
    subdf_up <- subdf[which(subdf$avg_log2FC > input$LogFCThreshold),]
    subdf_down <- subdf[which(subdf$avg_log2FC < -input$LogFCThreshold),]
    # subdf_notde <- subdf[which(subdf$avg_log2FC > -input$LogFCThreshold & 
    #                              subdf$avg_log2FC < input$LogFCThreshold),]
    
    top10up <- subdf_up %>% slice_max(n = 10, order_by = avg_log2FC)
    top10down <- subdf_down %>% slice_max(n = 10, order_by = -avg_log2FC)
    
    ggplotly(
      ggplot(subdf, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene)) +
        geom_point(data = subdf_up, color = "red") + 
        geom_point(data = subdf_down, color = "blue") + 
        geom_text(data = top10up, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene), nudge_y = 2) +
        geom_text(data = top10down, aes(x = avg_log2FC, y = -log10(p_val_adj), label = gene), nudge_y = 2) + theme_light()
      
    )
    
  })
  
  
  observeEvent(input$importtablemarkers, {
    
    tryCatch({
    
    if (file_ext(input$importtablemarkers$datapath) == "xlsx")
      val$markers <- readxl::read_xlsx(input$importtablemarkers$datapath, col_names = T)
    
    else if (file_ext(input$importtablemarkers$datapath) == "csv")
      val$markers <- read.csv(input$importtablemarkers$datapath, header = T, sep = ',')
    
    else{
      alert("Please enter a csv or a xlsx file")
      return(0)
    }
    
    }, error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      return(0)
    })
      
  })
  

}
