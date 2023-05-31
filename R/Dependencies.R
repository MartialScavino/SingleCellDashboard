if (!require("pacman")) install.packages("pacman")

packages <- c("shinydashboard", "shiny", "shinyWidgets", "shinyjs", "shinyFiles", "shinyBS", "Seurat",
              "ggplot2", "viridis", "cowplot", "DT", "babelgene", "plotly", "stringr", "scales", "tidyverse", "enrichR"
              , "rclipboard")

if (!require("shinysky")) devtools::install_github("AnalytixWare/ShinySky")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

pacman::p_load(char = packages)
df <- installed.packages()

if (!"UCell" %in% rownames(df))
  BiocManager::install("UCell", force = T, update = F)


if (!"DESeq2" %in% rownames(df))
  BiocManager::install("DESeq2", force = T, update = F)

library(shinysky)
library(dplyr)
library(tools)
library(UCell)
library(DESeq2)

poubelle <- sapply(c(packages, "UCell", "shinysky", "dplyr", "tools", "DESeq2"), function(x){
  
  if (!x %in% rownames(df)){
    
    print(paste0("Le Package ", x, " n'a pas été installé correctement"))
    return(0) 
  }
  
  return(1)
  
})

rm(poubelle)
rm(df)
rm(packages)


## Defining functions
helpPopup <- function(title, content,
                      placement=c('right', 'top', 'left', 'bottom'),
                      trigger=c('click', 'hover', 'focus', 'manual'),
                      width){
  if (trigger == 'hover'){ 
  
  tagList(
    singleton(
      tags$head(
        tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })"),
        tags$style(type = "text/css",
                   paste0(".popover{
                                    max-width:",width,";
                                  }"))
      )
    ),
    tags$i(
      href = "#", class = "fa-regular fa-circle-question", `data-toggle` = "popover",
      title = title, `data-content` = content, `data-animation` = TRUE,
      `data-placement` = match.arg(placement, several.ok=TRUE)[1],
      `data-trigger` = match.arg(trigger, several.ok=TRUE)[1],
      `data-html` = TRUE, `data-container` = "body"
    )
  )
  }
  
  else{
    tagList(
      singleton(
        tags$head(
          tags$script("$(function() { $(\"[data-toggle='popover']\").popover(); })"),
          tags$style(type = "text/css",
                     paste0(".popover{
                                    max-width:",width,";
                                  }"))
        )
      ),
      tags$a(
        href = "#", class = "fa-regular fa-circle-question", `data-toggle` = "popover",
        title = title, `data-content` = content, `data-animation` = TRUE,
        `data-placement` = match.arg(placement, several.ok=TRUE)[1],
        `data-trigger` = match.arg(trigger, several.ok=TRUE)[1],
        `data-html` = TRUE, `data-container` = "body"
      )
    )
    
  }
}

HelpInput <- function(content, text_placement = "bottom", icon_placement = "right", trigger = "hover", width = "500px"){
  
  style <- paste0("text-align: ", icon_placement)
  
  return(tags$div(helpPopup(title = NULL, content = content, placement = text_placement, trigger = trigger, width = width),
                  style= style))
}






