if (!require("rstudioapi")) install.packages("rstudioapi")

# Setting working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Install packages
source("R/Dependencies.R")

runApp('.')