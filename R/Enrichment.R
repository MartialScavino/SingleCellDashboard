enrichmentsidebar <- sidebarPanel(
  selectInput("database", "Choose a database", 
              c("Gene Ontology Biological Process" = "GO_Biological_Process_2021",
                "Gene Ontology Cellular Component" = "GO_Cellular_Component_2021",
                "Gene Ontology Molecular Function" = "GO_Molecular_Function_2021",
                "KEGG" = "KEGG_2021_Human", "Reactome" = "Reactome_2022")),
  sliderInput("pvaladj_enrichment", label = "Adjusted P-value threshold", min = 0, max = 1, value = 0.01, step = 0.01),
  textAreaInput("listgeneenrichment", "Input a list of genes (separated by a linebreak)", height = '300px'),
  actionButton("doenrichment", span("Compute", id = "UpdateAnimateEnrichment", class=""), styleclass = "primary")
  )

enrichmentmainpanel <- mainPanel(tabsetPanel(
  tabPanel("Dataframe",
           dataTableOutput("Enrichmentdt")),
  tabPanel("Dotplot",
           uiOutput("sliderinputenrichment"),
           plotOutput("Dotplot"))
))


enrichmentserver <- function(input, output, session, val){
  output$Enrichmentdt <- renderDataTable(val$enrichment, extensions = 'Buttons', server = F, 
                                         options = list(dom = 'Bfrtip', fixedColumns = TRUE, scrollX = TRUE,
                                                        buttons = c("copy", 'csv', 'excel')))
  
  
  observeEvent(input$doenrichment,  {
    
    if (input$listgeneenrichment == ""){
      alert("Please enter a list of gene")
      return(0)
    }
    
    addClass(id = "UpdateAnimateEnrichment", class = "loading dots")
    disable("doenrichment")
    
    tryCatch({
    string <- input$listgeneenrichment
    string <- sub(" ", "", string)
    liste_gene_enrich <- strsplit(string, "\n")[[1]]
    
    
    enrichmenttemp <- enrichr(liste_gene_enrich, input$database)[[1]]
    val$enrichment <- enrichmenttemp[which(enrichmenttemp$Adjusted.P.value < input$pvaladj_enrichment),]
    
    percentage <- lapply(val$enrichment$Overlap, function(i) {
      return(as.integer(unlist(str_split(i , "/")[[1]][1])) /as.integer(unlist(str_split(i , "/")[[1]][2])) )
    })
    val$enrichment[["Overlap_Percentage"]] <- unlist(percentage)
    
    val$enrichment <- val$enrichment %>% arrange(desc(Combined.Score))
    
    }, error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      enable("doenrichment")
      removeClass(id = "UpdateAnimateEnrichment", class = "loading dots")
      return(0)
      
    })
    enable("doenrichment")
    removeClass(id = "UpdateAnimateEnrichment", class = "loading dots")
    
  })
  
  
  output$sliderinputenrichment <- renderUI({
    if (is.null(val$data))
      return("")
    
    sliderInput("termtoplot", "Number of term to plot", 
                1, length(val$enrichment$Term), value = 30, step = 1)
  })
  
  output$Dotplot <- renderPlot({
    if (is.null(val$data))
      return(0)
    
    term_order <- val$enrichment$Term[order(val$enrichment$Combined.Score)]
    val$enrichment$Term <- factor(val$enrichment$Term, levels = term_order)
    
    ggplot(val$enrichment[1:input$termtoplot,], aes(x = Combined.Score, y = Term, size = Overlap_Percentage, color = Adjusted.P.value)) + 
      geom_point() + theme_light()
    
  })
  
}
