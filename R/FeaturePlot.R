singlegenesidebar <- sidebarPanel(
  uiOutput("typeaheadFeature"),
  hr(),
  uiOutput("selectsinglegene"))

singlegenemainpanel <- mainPanel(plotOutput("featureplotsinglegene"),
                                 plotOutput("violinsinglegene"))


signaturesidebar <- sidebarPanel(
  textInput("signaturename", "Signature name", value = ''),
  textAreaInput("signaturelist", "Gene list (separate each gene with a line break)", height = '300px'),
  textOutput("genenotpresent"),
  actionButton('dosignatureplot', span('Compute plot', id="UpdateAnimateFeature", class=""), styleclass = "primary"),
  hr(),
  uiOutput("selectsignature"))

signaturemainpanel <- mainPanel(plotOutput("featureplotsignature"),
                                plotOutput("violinsignature"))


## Feature Plot
featureplotserver <- function(input, output, session, val){
  
  
  ## SINGLE GENE
output$typeaheadFeature <- renderUI({
  
  genes <- data.frame(list(gene = rownames(GetAssay(val$data, "RNA"))))
  
  textInput.typeahead("gene", placeholder = "Enter a gene",
                      local = genes, 
                      valueKey = "gene",
                      tokens =  c(1:length(genes$gene)),
                      template = HTML("<p class='repo-language'>{{info}}</p> <p class='repo-gene'>{{gene}}</p>")
  )
  
})

output$selectsinglegene <- renderUI({
  if (is.null(val$data))
    return("")
  
  keep <- c()
  for (i in names(val$data@meta.data)){
    if (length(table(val$data@meta.data[,i])) < 50)
      keep <- append(keep, i)
  }
  selectInput("choice_singlegene", "Color by", choices = keep, selected = input$choice_singlegene)
  
})


output$featureplotsinglegene <- renderPlot({
  if (input$gene != "")
      FeaturePlot(val$data, features = input$gene)
})


output$violinsinglegene <- renderPlot({
  if (input$gene != "")
    VlnPlot(val$data, features = input$gene, group.by = input$choice_singlegene, assay = "RNA") + 
    scale_fill_manual(values = val$colors[[input$choice_singlegene]]) + 
    theme_light()
})

## SIGNATURE
output$selectsignature <- renderUI({
  if (is.null(val$data))
    return("")
  
  keep <- c()
  for (i in names(val$data@meta.data)){
    if (length(table(val$data@meta.data[,i])) < 50)
      keep <- append(keep, i)
  }
  selectInput("choice_signature", "Color by", choices = keep, selected = input$choice_singlegene)
  
})

output$genenotpresent <- renderText({
  
  string <- input$signaturelist
  string <- gsub(" ", "", string)
  genes <- strsplit(string, "\n")[[1]]
  
  if (length(genes) == 0)
    return("")
  
  genesnotpresent <- genes[which(!genes %in% rownames(GetAssay(val$data, "RNA")))]
  if (length(genesnotpresent) == 0)
    return("")
  
  return(paste0("The following genes are not present with this name in the dataset : ",toString(genesnotpresent)))
  
})



# Calcule le plot quand le bouton est sélectionné 
observeEvent(input$dosignatureplot,{
  
  if (input$signaturelist == ""){
    alert("Please enter a list of gene")
    return(0)
  }
  
  addClass(id = "UpdateAnimateFeature", class = "loading dots")
  disable("dosignatureplot")
  
  
  tryCatch({
    
  string <- input$signaturelist
  string <- sub(" ", "", string)
  liste_gene_plot <- strsplit(string, "\n")[[1]] 
  
  if (input$signaturename == ""){
    val$data <- AddModuleScore_UCell(val$data, features = list(Signature = liste_gene_plot), assay = "RNA")
    title <- "Signature"
    p <- FeaturePlot(object = val$data, features = "Signature_UCell") + ggtitle("Signature")
    
    output$violinsignature <- renderPlot({
      VlnPlot(val$data, features = "Signature_UCell", group.by = input$choice_signature, assay = "RNA") + 
        scale_fill_manual(values = val$colors[[input$choice_signature]]) + 
        theme_light() + ggtitle(title)
    })
    
  }
  
  else{
    name <- sub(" ", "", input$signaturename)
    val$data <- AddModuleScore_UCell(val$data ,features = list(Signature_ = liste_gene_plot), name = paste0(name, "_UCell"), assay = "RNA")
    title <- input$signaturename
    p <- FeaturePlot(val$data, features = paste0("Signature_", name, "_UCell")) + ggtitle(input$signaturename)
    
    output$violinsignature <- renderPlot({
      VlnPlot(val$data, features = paste0("Signature_", name, "_UCell"), group.by = input$choice_signature, assay = "RNA") + 
        scale_fill_manual(values = val$colors[[input$choice_signature]]) + 
        theme_light() + ggtitle(title)
    })
    
  }
  
  
  output$featureplotsignature <- renderPlot(p)
  
  enable("dosignatureplot")
  removeClass(id = "UpdateAnimateFeature", class = "loading dots")
  
  }, error = function(e){
    alert("There has been an error (printed in R console)")
    print(e)
    enable("dosignatureplot")
    removeClass(id = "UpdateAnimateFeature", class = "loading dots")
    return(0)
  })
  
})

}