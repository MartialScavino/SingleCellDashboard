# QC panel
sidebar_QC <- sidebarPanel(
  sliderInput("mt", "Mamimum percentage of mitochondrial genes",
              min = 0,
              max = 100,
              value = 20,
              step = 1),
  uiOutput("sliderFeature"),
  textOutput("texte"),
  actionButton("Trim", span("Compute trimming", id="UpdateAnimateTrim", class=""), styleclass = "primary"))


qc <- tabItem(tabName = "qc",
              sidebarLayout(sidebarPanel = sidebar_QC, 
                            mainPanel = mainPanel(
                              tabsetPanel(
                                tabPanel("Scatter", 
                                         plotOutput("scatter_QC_MT"),
                                         plotOutput("scatter_QC_Feature")
                                ),
                                
                                tabPanel("Hist", 
                                         plotOutput("hist_QC_MT"),
                                         plotOutput("hist_QC_Feature")),
                                
                                tabPanel("Violin", 
                                         uiOutput("select_violin_QC"),
                                         plotOutput("violin_QC_MT"),
                                         plotOutput("violin_QC_Features")),
                                tabPanel("Stacked",
                                         fluidRow(
                                           uiOutput("select_stacked_1"),
                                           uiOutput("select_stacked_2")),
                                         fluidRow(
                                           plotlyOutput("stacked_barplot")))
                              )
                            )
              )
)

QC <- qcserver <- function(input, output, session, val){
  
  
  # Slider featre slide bar
  output$sliderFeature <- renderUI(sliderInput("features",
                                               "Number of features",
                                               min = 0,
                                               max = max(val$data$nFeature_RNA),
                                               value = c(600, max(val$data$nFeature_RNA)), 
                                               step = 50))
  
  # Number of cell removed text
  output$texte <- renderText({
    
    nb_cell_bad <- length(rownames(val$data@meta.data[which(val$data$percent.mt > input$mt | 
                                                              val$data$nFeature_RNA < input$features[1] | 
                                                              val$data$nFeature_RNA > input$features[2]),]))
    nb_cell_tot <- length(rownames(val$data@meta.data))
    
    paste0("Removing ", nb_cell_bad, " cells from the dataset (", round(100*nb_cell_bad/nb_cell_tot, 2), "% of total cells)")
    
  })
  
  output$scatter_QC_MT <- renderPlot({
    
    
    ggplot(val$data@meta.data, aes(nCount_RNA, percent.mt, color = nFeature_RNA)) +
      geom_point(data = val$data@meta.data[val$data$percent.mt > input$mt,], alpha = 0.2) + 
      geom_point(data = val$data@meta.data[val$data$percent.mt <= input$mt,]) + 
      geom_hline(yintercept = input$mt, col = "#FEB078") + 
      scale_color_viridis(option = "B") + theme_light()
    
    
  })
  
  output$scatter_QC_Feature <- renderPlot({
    
    ggplot(val$data@meta.data, aes(nCount_RNA, nFeature_RNA, color = percent.mt)) +
      geom_point(data = val$data@meta.data[val$data$nFeature_RNA < input$features[1],], alpha = 0.2) + 
      geom_point(data = val$data@meta.data[val$data$nFeature_RNA > input$features[2],], alpha = 0.2) + 
      geom_point(data = val$data@meta.data[val$data$nFeature_RNA <= input$features[2] & val$data$nFeature_RNA >= input$features[1],]) + 
      scale_color_viridis() +
      geom_hline(yintercept = input$features, col = "#832681") + theme_light()
    
  })
  
  
  output$hist_QC_MT <- renderPlot({
    
    ggplot(val$data@meta.data, aes(x = percent.mt)) + 
      geom_histogram(aes(y = after_stat(density)), bins = 200, fill = "#FEB078") + 
      geom_density(color = "#f8765c") + 
      geom_vline(xintercept = input$mt, col = "#800000") + theme_light()
    
  })
  
  output$hist_QC_Feature <- renderPlot({
    
    ggplot(val$data@meta.data, aes(x = nFeature_RNA)) + 
      geom_histogram(aes(y = after_stat(density)), bins = 200, fill = "#DDA0DD") + 
      geom_density(color = "#8B0000") + 
      geom_vline(xintercept = input$features, col = "#832681") + theme_light()
    
  })
  
  
  output$select_violin_QC <- renderUI({
    keep <- c()
    for (i in names(val$data@meta.data)){
      if (length(table(val$data@meta.data[,i])) < 50)
        keep <- append(keep, i)
    }
    selectInput("choice_violin_QC", "Select a variable", choices = keep)
    
  })
  
  
  output$violin_QC_MT <- renderPlot({
    VlnPlot(val$data, "nFeature_RNA", group.by = input$choice_violin_QC) + 
      geom_hline(yintercept = input$features, col = "#832681") + NoLegend()
  })
  
  output$violin_QC_Features <- renderPlot({
    VlnPlot(val$data, "percent.mt", group.by = input$choice_violin_QC) + 
      geom_hline(yintercept = input$mt, col = '#FEB078') + NoLegend()
    
  })
  
  
  output$select_stacked_1 <- renderUI({
    keep <- c()
    for (i in names(val$data@meta.data)){
      if (length(table(val$data@meta.data[,i])) < 50)
        keep <- append(keep, i)
    }
    selectInput("choice_stacked_x", "Split by", choices = keep)
    
  })
  
  
  output$select_stacked_2 <- renderUI({
    keep <- c()
    for (i in names(val$data@meta.data)){
      if (length(table(val$data@meta.data[,i])) < 50)
        keep <- append(keep, i)
    }
    selectInput("choice_stacked_y", "Color by", choices = keep)
    
  })
  
  
  output$stacked_barplot <- renderPlotly({
    
    ggplotly(ggplot(val$data@meta.data, aes(fill=.data[[input$choice_stacked_y]], 
                                            x= .data[[input$choice_stacked_x]])) + 
               geom_bar(position="fill"))
    
    
  })

}
