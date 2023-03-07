if (!require("pacman")) install.packages("pacman")

packages <- c("shinydashboard", "shiny", "shinyWidgets", "shinyjs", "shinyFiles", "shinyBS", "Seurat",
              "ggplot2", "viridis", "cowplot", "DT", "babelgene", "plotly", "stringr", "scales", "tidyverse", "enrichR", "rclipboard")

if (!require("shinysky")) devtools::install_github("AnalytixWare/ShinySky")

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")

pacman::p_load(char = packages)
df <- installed.packages()

if (!"UCell" %in% rownames(df))
  BiocManager::install("UCell", force = T, update = F)

library(shinysky)
library(dplyr)
library(tools)
library(UCell)

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