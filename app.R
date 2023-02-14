if (!require("rstudioapi")) install.packages("rstudioapi")

# Setting working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(shinydashboard)
library(shiny)
library(shinysky)
library(shinyWidgets)
library(shinyjs)
library(shinyFiles)
library(shinyBS)
library(Seurat)
library(ggplot2)
library(viridis)
library(cowplot)
library(DT)
library(babelgene)
library(UCell)
library(plotly)
library(stringr)
library(scales)
library(dplyr)


runApp('.')


# test <- readRDS("../../test.rds")
# Idents(test) <- "orig.ident"
# markers <- FindAllMarkers(test, only.pos = T)
# 
# markers <- markers[which(markers$p_val_adj < 0.01),]
# top10test <- markers %>% 
#   group_by(cluster) %>% 
#   slice_max(n = 5, order_by = avg_log2FC)
# 
# 
# ### MARCHE
# factors <- c("ADK_02032", "ADK_02036", "ADK_02043",
#                                                 "Polype_02032", "Polype_02036", "Polype_02043",
#                                                 "TS_02032", "TS_02036", "TS_02043")
# 
# test$orig.ident <- factor(test$orig.ident, levels = factors)
# 
# top10test$cluster <- factor(top10test$cluster, levels = factors)
# 
# top10test <- top10test %>% arrange(factor(cluster, levels = factors))
# 
# DotPlot(test, features = unique(top10test$gene), group.by = "orig.ident") +
#   theme(axis.text.x = element_text(angle = 60, size = 8, vjust = 0.85))

#levels(test$orig.ident)                                    
