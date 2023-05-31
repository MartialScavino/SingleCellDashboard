options(shiny.maxRequestSize = 30000 * 1024^2)
set.seed(123)

Animation <- tags$head(tags$style(type="text/css", '
            .loading {
                display: inline-block;
                overflow: hidden;
                height: 1.3em;
                margin-top: -0.3em;
                line-height: 1.5em;
                vertical-align: text-bottom;
                box-sizing: border-box;
            }
            .loading.dots::after {
                text-rendering: geometricPrecision;
                content: "⠋\\A⠙\\A⠹\\A⠸\\A⠼\\A⠴\\A⠦\\A⠧\\A⠇\\A⠏";
                animation: spin10 1s steps(10) infinite;
                animation-duration: 1s;
                animation-timing-function: steps(10);
                animation-delay: 0s;
                animation-iteration-count: infinite;
                animation-direction: normal;
                animation-fill-mode: none;
                animation-play-state: running;
                animation-name: spin10;
            }
            .loading::after {
                display: inline-table;
                white-space: pre;
                text-align: left;
            }
            @keyframes spin10 { to { transform: translateY(-15.0em); }}
            '))




## UI
changelabel <- tabPanel("Change labels",
                        fluidRow(
                          box(uiOutput("SelectLabelInf50")),
                          box(textInput("newcolumn", value = "", "Set the new column name (it can be an existing one)"),
                              actionButton("dochangelabel", "Apply changes", styleclass = "primary"))),
                        fluidRow(
                          
                          box(title = "Old label",width = 6,
                              uiOutput("listOldLabel")),
                          box(title = "New label", width = 6,
                              uiOutput("listNewLabel"))
                          
                        )
                        
)


removecolumn <- tabPanel("Remove column",
                         fluidRow(
                           box(uiOutput("colnames")),
                           box(actionButton("doremove", "Remove column", styleclass = "primary"))
                         ),
                         fluidRow(dataTableOutput("headmetadata"))
)


changecolors <- tabPanel("Change colors",
                         fluidRow(
                           box(uiOutput("SelectLabelInf50_2"),
                               actionButton("resetcolors", "Reset colors of selected column", styleclass = "danger")),
                           box(actionButton("dochangecolors", "Apply changes", styleclass = "primary"))
                           
                         ),
                         fluidRow(
                           
                           box(title = "Label",width = 6,
                               uiOutput("listLabel")),
                           box(title = "Color", width = 6,
                               uiOutput("listColor"))
                         ))


changeorder <- tabPanel("Change order",
                        fluidRow(
                          box(uiOutput("SelectLabelInf50_3")),
                          box(actionButton("dochangeorder", "Apply changes", styleclass = "primary"))
                        ),
                        fluidRow(
                          box(title = "Label", width = 6,
                              uiOutput("listLabelOrder")),
                          box(title = "Order", width = 6,
                              uiOutput("listOrder"))
                        ))


## SERVER

optionsserver <- function(input, output, session, val){
  
  output$SelectLabelInf50 <- renderUI({
    if (is.null(val$data))
      return("")
    
    keep <- c()
    for (i in names(val$data@meta.data)){
      if (length(table(val$data@meta.data[,i])) < 50)
        keep <- append(keep, i)
    }
    selectInput("changeLabel", "Select a variable", choices = keep)
    
  })
  
  
  
  observeEvent(input$changeLabel, {
    output$listOldLabel <- renderUI({
      if (is.null(val$data))
        return("")
      
      tryCatch({
      
      tagList(
        lapply(1:length(table(val$data@meta.data[,input$changeLabel])), 
               function(i){
                 
                 id <- paste0('oldLabel_', i)
                 disabled(
                   textInput(id, "", value = names(table(val$data@meta.data[,input$changeLabel]))[i])
                 )
                 
               })
      )
      }, error = function(e){
        alert("There has been an error (printed in R console)")
        print(e)
        return(0)
      })
    })
    
    output$listNewLabel <- renderUI({
      if (is.null(val$data))
        return("")
      
      tagList(
        lapply(1:length(table(val$data@meta.data[,input$changeLabel])), 
               function(i){
                 
                 id <- paste0('newLabel_', i)
                 textInput(id, "", value = names(table(val$data@meta.data[,input$changeLabel]))[i])
                 
               })
      )
    })
  })
  
  
  
  ## Apply changes
  observeEvent(input$dochangelabel, {
    
    alreadyexist = T
    
    if (input$newcolumn == ""){
      alert("Enter a column name")
      return(0)
    }
    
    if (!(input$newcolumn %in% colnames(val$data@meta.data))){
      alreadyexist = F
      val$data@meta.data[[input$newcolumn]] <- val$data@meta.data[, input$changeLabel]
    }
    
    lapply(1:length(table(val$data@meta.data[,input$changeLabel])),
           function(i){
             oldid  <- paste0('oldLabel_', i)
             newid <- paste0('newLabel_', i)
             
             if (alreadyexist)
               val$data@meta.data[, input$newcolumn] <- gsub(input[[oldid]], 
                                                             input[[newid]], 
                                                             val$data@meta.data[, input$changeLabel])
             
             else{
               
               val$data@meta.data[, input$newcolumn] <- gsub(input[[oldid]], 
                                                             input[[newid]], 
                                                             val$data@meta.data[, input$newcolumn])
             }
           } 
    )
  })
  
  
  ## Remove column
  output$colnames <- renderUI({
    if (is.null(val$data))
      return("")
    
    selectInput("coltoremove", "Column to remove",
                choices = colnames(val$data@meta.data))
    })
  
  output$headmetadata <- renderDataTable(
    if (!(is.null(val$data)))
      val$data@meta.data, 
    options = list(pageLength = 10, width="100%", scrollX = TRUE))
  
  observeEvent(input$doremove, 
               val$data@meta.data <- val$data@meta.data[, -which(colnames(val$data@meta.data) == input$coltoremove)])
  
  
  
  observeEvent(val$data,{
    
    sapply(colnames(val$data@meta.data), function(i){
      
      if (length(names(table(val$data@meta.data[,i]))) < 50 & 
          length(names(table(val$data@meta.data[,i]))) != 0){
        
        if (!(i %in% names(val$colors)))
          val$colors[[i]] <- hue_pal()(length(names(table(val$data@meta.data[,i]))))
        
        else{
          
          if (length(val$colors[[i]]) != length(names(table(val$data@meta.data[,i]))))
            val$colors[[i]] <- hue_pal()(length(names(table(val$data@meta.data[,i]))))
          
        }
      }
    })
  })
  
  
  output$SelectLabelInf50_2 <- renderUI({
    if (is.null(val$data))
      return("")
    
    keep <- c()
    for (i in names(val$data@meta.data)){
      if (length(table(val$data@meta.data[,i])) < 50)
        keep <- append(keep, i)
    }
    selectInput("selectchangecolor", "Select a variable", choices = keep)
    
  })
  
  
  observeEvent(input$selectchangecolor, {
    output$listLabel <- renderUI({
      if (is.null(val$data))
        return("")
      
      tagList(
        lapply(1:length(table(val$data@meta.data[,input$selectchangecolor])), 
               function(i){
                 
                 id <- paste0('Label_', i)
                 disabled(
                   textInput(id, "", value = names(table(val$data@meta.data[,input$selectchangecolor]))[i])
                 )
                 
               })
      )
    })
    
    
    output$listColor <- renderUI({
      if (is.null(val$data))
        return("")
      
      tagList(
        lapply(1:length(table(val$data@meta.data[,input$selectchangecolor])), 
               function(i){
                 
                 id <- paste0('Label_', i)
                 textInput(id, "", value = unlist(val$colors[input$selectchangecolor])[i])
                 
               })
      )
    })
    
    
    
  })
  
  
  observeEvent(input$resetcolors,{
    tryCatch({
      
      val$colors[[input$selectchangecolor]] <- hue_pal()(length(names(table(val$data@meta.data[,input$selectchangecolor]))))
      },
      error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      return(0)
  } 
  )})
  
  
  observeEvent(input$dochangecolors, {
    
    tryCatch({
    lapply(1:length(table(val$data@meta.data[, input$selectchangecolor])),
           function(i){
             id <- paste0('Label_', i)
             
             val$colors[[input$selectchangecolor]][i] <- input[[id]]
           }
    )
    }, error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      return(0)
    })
  })
  
 
  output$SelectLabelInf50_3 <- renderUI({
    if (is.null(val$data))
      return("")
    
    keep <- c()
    for (i in names(val$data@meta.data)){
      if (length(table(val$data@meta.data[,i])) < 50)
        keep <- append(keep, i)
    }
    selectInput("selectchangeorder", "Select a variable", choices = keep, selected = input$selectchangeorder)
    
  })
  
  
  observeEvent(input$selectchangeorder, {
    output$listLabelOrder <- renderUI({
      if (is.null(val$data))
        return("")
      
      tryCatch({
      
      tagList(
        lapply(1:length(table(val$data@meta.data[,input$selectchangeorder])), 
               function(i){
                 
                 id <- paste0('LabelOrder_', i)
                 disabled(
                   textInput(id, "", value = names(table(val$data@meta.data[,input$selectchangeorder]))[i])
                 )
                 
               })
      )
      }, error = function(e){
        alert("There has been an error (printed in R console)")
        print(e)
        return(0)
      })
    })
    
    
    output$listOrder <- renderUI({
      if (is.null(val$data))
        return("")
      
      tagList(
        lapply(1:length(table(val$data@meta.data[,input$selectchangeorder])), 
               function(i){
                 
                 id <- paste0('Order_', i)
                 textInput(id, "", 
                           value = which(names(table(val$data@meta.data[,input$selectchangeorder])) == names(table(val$data@meta.data[,input$selectchangeorder]))[i]))
                 
               })
      )
    })
    
  })
  
  
  ## Faire les 3 boutons 
  
  observeEvent(input$dochangeorder, {
    
    tryCatch({
    
    df <- data.frame(list(label = "", order = ""))
    df_colors <- data.frame(list(color = '', order = ''))
    
    for (i in 1:length(names(table(val$data@meta.data[,input$selectchangeorder])))){
      
      id_labelorder <- paste0('LabelOrder_', i)
      id_order <- paste0("Order_", i)
      
      df[i,] <- c(input[[id_labelorder]], input[[id_order]])
      df_colors[i,] <- c(val$colors[[input$selectchangeorder]][i], input[[id_order]])
      
    }
    
    df <- df %>% arrange(order)
    df_colors <- df_colors %>% arrange(order)
    
    levels_test <- as.vector(df$label)
    
    val$data@meta.data[,input$selectchangeorder] <- factor(val$data@meta.data[,input$selectchangeorder], levels = levels_test)
    
    if (all(val$colors[[input$selectchangeorder]] != hue_pal()(length(names(table(val$data@meta.data[, input$selectchangeorder]))))))
      val$colors[[input$selectchangeorder]] <- as.vector(df_colors$color)
    }, error = function(e){
      alert("There has been an error (printed in R console)")
      print(e)
      return(0)
    })
  })
  
  
   
}
