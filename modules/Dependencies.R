if (!require("pacman")) install.packages("pacman")

packages <- c("shinydashboard", "shiny", "shinyWidgets", "shinyjs", "shinyFiles", "shinyBS", "Seurat",
              "ggplot2", "viridis", "cowplot", "DT", "babelgene", "plotly", "stringr", "scales", "tidyverse")

if (!require("shinysky")) devtools::install_github("AnalytixWare/ShinySky")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

if (require("UCell"))
  BiocManager::install("UCell", force = T, update = F)

pacman::p_load(char = packages)

library(shinysky)
library(dplyr)
library(tools)
library(UCell)

df <- installed.packages()
poubelle <- sapply(c(packages, "UCell", "shinysky", "dplyr", "tools"), function(x){
  
  if (!x %in% rownames(df)){
    
    print(paste0("Le Package ", x, " n'a pas été installé correctement"))
    return(0) 
  }
  
  return(1)
  
})

rm(poubelle)
rm(df)
rm(packages)