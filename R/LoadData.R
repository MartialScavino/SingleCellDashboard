# Loading data
load_data <- tabItem(tabName = "load_data",
                     fluidRow(
                       
                       box(title = "RDS", fileInput(inputId = "data", label = "Enter a RDS file", accept = "")),
                       
                       box(title = "filtered_feature_bc_matrix", 
                           fileInput(inputId = "barcodesfile", "barcodes.tsv.gz"),
                           fileInput(inputId = "featuresfile", "features.tsv.gz"),
                           fileInput(inputId = "matrixfile", "matrix.mtx.gz"),
                           actionButton("filtered", span("Load Data", id = "UpdateAnimateLoad", class = ""), styleclass = "primary"))),
                     fluidRow(
                       box(title = "Metadata", width = 12,
                           dataTableOutput("dataset")
                       )
                     )
)



LoadDataServer <- function(input, output, session, val){
  
  # Loading dataset
  observeEvent(input$data,{
    # Check for file type
    if (file_ext(input$data$datapath) != "rds"){ 
      alert("Enter a rds file")
      return(0)
    }
    
    tryCatch({
      df <- readRDS(input$data$datapath)
      
      if (!("percent.mt" %in% names(df@meta.data))){
        
        if (length(grep( "^mt-", rownames(df), value = T) > 0))
          df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^mt-")
        
        else if (length(grep( "^MT-", rownames(df), value = T) > 0))
          df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
      }  
      
      val$data <- df
      
    }, error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      return(0)
    }) 
    
    # Show Metadata
    output$dataset <- renderDataTable({
      if (is.null(val$data))
        return(data.table(list()))
      
      val$data@meta.data
    } , options = list(dom = 'frtip', 
                       fixedColumns = TRUE, 
                       scrollX = T))})
  
  # Load data with 10X files
  observeEvent(input$filtered, {
    
    addClass(id = "UpdateAnimateLoad", class = "loading dots")
    disable("filtered")
    
    tryCatch({
      
      counts <- ReadMtx(input$matrixfile$datapath,
                        cells = input$barcodesfile$datapath,
                        features = input$featuresfile$datapath,
                        feature.column = 1)
      
      df <- CreateSeuratObject(counts = counts)
      
      if (!("percent.mt" %in% names(df@meta.data))){
        
        if (length(grep( "^mt-", rownames(df), value = T) > 0))
          df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^mt-")
        
        else if (length(grep( "^MT-", rownames(df), value = T) > 0))
          df[["percent.mt"]] <- PercentageFeatureSet(df, pattern = "^MT-")
      }
      
      val$data <- df
      
    }, error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      enable("filtered")
      removeClass(id = "UpdateAnimateLoad", class = "loading dots")
      return(0)
    })
    
    # Show Metadata
    output$dataset <- renderDataTable({
      if (is.null(val$data))
        return(data.frame(list()))
      
      val$data@meta.data
    } , options = list(dom = 'frtip', 
                       fixedColumns = TRUE, 
                       scrollX = T))
    
    enable("filtered")
    removeClass(id = "UpdateAnimateLoad", class = "loading dots")
  })
  
  
  
}