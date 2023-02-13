# Options done in modules/Options.R
## Set the max size to import data
## Describe the spinner when clicking on a button

# Sourcing modules
source("modules/Options.R")
source("modules/FeaturePlot.R")
source("modules/QC.R")


head <- dashboardHeader(title = "Single cell analysis", tags$li(class = "dropdown", 
                                                                downloadButton("savedata", span("Save data", id = "UpdateAnimateSave", class = ""))))

side <- dashboardSidebar(
  sidebarMenu(
    menuItem("Load data", tabName ="load_data", icon = icon("upload")),
    menuItem("QC", tabName = 'qc', icon = icon("circle-check")),
    menuItem("Preprocessing", tabName = "preprocessing"),
    menuItem("Dimension reduction", tabName = "dimreduc", 
             menuSubItem("PCA", tabName = "pca"),
             menuSubItem("UMAP", tabName = "umap")),
    menuItem("Clustering", tabName = "clustering"),
    menuItem("Options", tabName = "options", icon = icon("cog")),
    menuItem("FeaturePlot", tabName = "featureplot")

  )
)

# Loading data
load_data <- tabItem(tabName = "load_data",
                     fluidRow(
                       
                       box(title = "RDS", fileInput(inputId = "data", label = "Select the path to your Rds", accept = "")),
                       
                       box(title = "filtered_feature_bc_matrix", 
                           fileInput(inputId = "barcodesfile", "barcodes.tsv.gz"),
                           fileInput(inputId = "featuresfile", "features.tsv.gz"),
                           fileInput(inputId = "matrixfile", "matrix.mtx.gz"),
                           actionButton("filtered", span("Load Data", id = "UpdateAnimateLoad", class = ""), styleclass = "primary"))),
                     fluidRow(dataTableOutput("dataset"))
)
 

preprocess <- tabItem(tabName = "preprocessing",
                      sidebarLayout(
                        sidebarPanel = sidebarPanel(
                          radioButtons("species", "Species", c("Human", "Mouse"),"Human"),
                          hr(),
                          selectInput("NormMethod", "Normalization method",
                                      choices = c("LogNormalize" = "LogNormalize", 
                                                  "Centered log ratio" = "CLR",
                                                  "Relative counts" = "RC"), 
                                      selected = "LogNormalize"),
                          numericInput("ScaleFactor", "Scale factor",
                                       value = 10000,
                                       min = 1000, max = 1000000, step = 1000),
                          hr(),
                          selectInput("VariableMethod", "Variable features selection method",
                                      choices = c("Vst" = "vst", "Dispersion" = "dispersion"),
                                      selected = "Vst"),
                          numericInput("nfeatures", "Number of variable features",
                                       value = 2000, 
                                       min = 1000, max = 10000, step = 500),
                          hr(),
                          selectInput("FeatureScale", "Set of genes to use for scaling",
                                      choices = c("Variable genes", "All genes")),
                          checkboxInput("regressing", "Regress cell cycle"),
                          hr(),
                          actionButton("LaunchPreprocessing", span("Compute preprocessing", id="UpdateAnimatePreprocessing", class=""), styleclass = "primary")
                        ),
                        mainPanel = mainPanel(
                          
                          plotOutput("variablefeatureplot")
                        )
                      ))


# Dimension reduction
pca <- tabItem(tabName = "pca",
        sidebarLayout(
          sidebarPanel = sidebarPanel(
            selectInput("featurePCA", "Genes to use",
                        choices = c("Variable genes", "All genes"),
                        selected = "Variable genes"),
            actionButton("dopca", span("Compute PCA", id="UpdateAnimatePCA", class=""), styleclass = "primary")
          ),
          mainPanel = mainPanel(
            tabsetPanel(
              tabPanel("Visualisation",
                       uiOutput("groupbypca"),
                       plotOutput("PCA")
                       ),
              tabPanel("ElbowPlot",
            sliderInput("ndimElbow", label = "Number of dimension to plot",
                         10, 50, 30, 5),
            plotOutput("elbow")),
            
              tabPanel("DimHeatmap",
            sliderInput("ndimheatmap", "Number of dimension to plot",
                        10, 30, 15, 5),
            plotOutput("dimheatmap"))
            )
          )
        )
)

umap <- tabItem(tabName = "umap", 
                sidebarLayout(
                  sidebarPanel = sidebarPanel(
                    numericInput("ndimumap", "Number of dimension",
                                 30, 1, 200),
                    actionButton("doumap", span("Compute UMAP", id = "UpdateAnimateUMAP", class=""), styleclass = "primary")
                  ),
                  
                  mainPanel = mainPanel(
                    uiOutput("groupbyumap"),
                    plotOutput("UMAP")
                  )
                ))


clustering <- tabItem(tabName = "clustering",
                      sidebarLayout(
                        sidebarPanel = sidebarPanel(
                          numericInput("kparam", "k parameter", 
                                       20, 1, 100),
                          uiOutput("sliderndimclusters"),
                          numericInput("resolution", "Resolution", 
                                       0.8, 0, 5, 0.1),
                          selectInput("algocluster", "Algorithm", 
                                      choices = c("Louvain" = "1", 
                                                  "Louvain with multilevel refinement" = "2",
                                                  "SLM" = "3",
                                                  "Leiden" = "4")),
                          actionButton("docluster", span("Compute clustering", id = "UpdateAnimateCluster", class = ""), styleclass = "primary")
                        ),
                        mainPanel = mainPanel(
                          uiOutput("groupbycluster"),
                          plotOutput("UMAPCluster")
                        )
                      ))



# Donz in modules/options.R
option <- tabItem(tabName = "options",
                  tabsetPanel(
                    changelabel,
                    removecolumn,
                    changecolors
                  ))

featureplot <- tabItem(tabName = "featureplot", 
                       ## Done in modules/FeaturePlot.R
                       featureplotsidebar, 
                       featureplotmain)

bod <- dashboardBody(
  #chooseSliderSkin(skin = "Shiny"),
  #setSliderColor(color = c("#FEB078", "#832681"),sliderId =  c(1, 2)),
  useShinyjs(), 
  
  ## Done in modules/Options.R
  Animation,
  tabItems(
    load_data,
    qc,
    preprocess,
    pca,
    umap,
    clustering,
    option,
    featureplot
  )
)



dashboardPage(
  skin = "black",
  header = head,
  sidebar = side,
  body = bod
)

