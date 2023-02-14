featureplotsidebar <- sidebarPanel(
  textInput("signaturename", "Signature name", value = ''),
  p(markdown('---')),
  p("Signature genes"),
  p("(Click on a button to remove it)"),
  uiOutput("list_button"),
  p(markdown('---')),
  p(HTML("<br><br><br>")),
  actionButton('dofeature', span('Compute plot', id="UpdateAnimateFeature", class=""), styleclass = "primary"))

featureplotmain <-mainPanel(
  uiOutput("typeaheadFeature"),
  actionButton("reset", "Reset signature", styleclass = "danger"),
  plotOutput("featuresplot"))


## Feature Plot
featureplotserver <- function(input, output, session, val){
  
output$typeaheadFeature <- renderUI({
  
  genes <- data.frame(list(gene = rownames(val$data)))
  
  textInput.typeahead("gene", placeholder = "Enter a gene",
                      local = genes, 
                      valueKey = "gene",
                      tokens =  c(1:length(genes$gene)),
                      template = HTML("<p class='repo-language'>{{info}}</p> <p class='repo-gene'>{{gene}}</p>")
  )
  
})

liste = reactiveValues(df = data.frame(list(gene_name = character(0), show_button = character(0))))

# Ajoute un gène à la liste lorsqu'on écrit un nouveau gène
observeEvent(input$gene,{
  
  if (input$gene != "" & input$gene %in% liste$df$gene_name == FALSE){
    liste$df[nrow(liste$df) + 1,] <- c(input$gene, TRUE)
  }
  
  else if (input$gene != "" & input$gene %in% liste$df$gene_name == TRUE){
    if (liste$df[which(liste$df$gene_name == input$gene), "show_button"] == FALSE){
      
      liste$df[which(liste$df$gene_name == input$gene), "show_button"] <- TRUE
      
    }
  }
})


# Créer un bouton pour chaque élément dans la liste de gène signature
observeEvent(liste$df$gene_name, output$list_button <- renderUI({
  
  if (length(liste$df$gene_name) > 0){
    tagList(
      lapply(1:length(liste$df$gene_name), function(i){
        
        if (liste$df$show_button[i] == TRUE){
          
          id1 <- paste0('slider_',liste$df$gene_name[i])
          actionButton(id1, liste$df$gene_name[i], "info", 
                       onclick = "Shiny.onInputChange('myclick', {id : this.id, val : this})")
        }
      })
    )
  }
  
  else{
    
    output$list_button <- renderUI(actionButton('bla', "No gene selected", "inverse"))
    
  }
}))


# Reset list quand on appuie sur le bouton
observeEvent(input$reset, {
  
  liste$df <- data.frame(list(gene_name = character(0), show_button = character(0)))
  
})

# Make button removable when clicking on it
observeEvent(input$myclick, {
  
  id <- strsplit(input$myclick$id, "_")[[1]][2]
  liste$df[which(liste$df$gene_name == id), "show_button"] <- FALSE
  
})


# Calcule le plot quand le bouton est sélectionné 
observeEvent(input$dofeature,{
  output$featuresplot <- renderPlot({
    
    isolate({
      
      liste_gene_plot <- liste$df[which(liste$df$show_button == TRUE), "gene_name"]
      
      if (length(liste_gene_plot) == 1){
        
        p <- FeaturePlot(val$data, features = liste_gene_plot)
        return(p)
      }
      
      else if (length(liste_gene_plot) > 1){
        
        addClass(id = "UpdateAnimateFeature", class = "loading dots")
        disable("dofeature")
        
        if (input$signaturename == ""){
          val$data <- AddModuleScore_UCell(val$data, features = list(Signature = liste_gene_plot))
          p <- FeaturePlot(object = val$data, features = "Signature_UCell") + ggtitle("Signature")
        }
        
        else{
          name <- sub(" ", "", input$signaturename)
          val$data <- AddModuleScore_UCell(val$data ,features = list(Signature_ = liste_gene_plot), name = paste0(name, "_UCell"))
          p <- FeaturePlot(val$data, features = paste0("Signature_", name, "_UCell")) + ggtitle(input$signaturename)
        }
        
        enable("dofeature")
        removeClass(id = "UpdateAnimateFeature", class = "loading dots")
        return(p)
        
      }
      
      else { return() }
      
    })
  })
})

}