if (!require("rstudioapi")) install.packages("rstudioapi")

# Setting working directory
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("modules/Dependencies.R")

runApp('.')
