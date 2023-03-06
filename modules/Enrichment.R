enrichmentsidebar <- sidebarPanel(
  selectInput("database", "Choose a database", 
              c("Gene Ontology Biological Process" = "GO_Biological_Process_2021",
                "Gene Ontology Cellular Component" = "GO_Cellular_Component_2021",
                "Gene Ontology Molecular Function" = "GO_Molecular_Function_2021",
                "KEGG" = "KEGG_2021_Human", "Reactome" = "Reactome_2022")),
  textAreaInput("listgeneenrichment", "Input a list of genes (separated by a linebreak)", height = '300px'),
  actionButton("doenrichment", span("Compute", id = "UpdateAnimateEnrichment", class=""), styleclass = "primary")
  )

enrichmentmainpanel <- mainPanel()


enrichmentserver <- function(input, output, session, val){
  
# dbs[order(dbs$libraryName), "libraryName"]
  
  
}
