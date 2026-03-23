
library(shiny)
library(ggplot2)
library(dplyr)
library(viridisLite)
library(ggpubr)

#----# #----# #----# #----# #----# #----# #----# #----# 
#----# #----# #----#  HOW TO RUN   #----# #----# #----# 
#                                                     #
#   load needed libraries                             #
#   click on `> Run App` in the command bar           #
#                                                     # 
#----# #----# #----# #----# #----# #----# #----# #----# 



addResourcePath(prefix = "neuron_maps_pdf", directoryPath = "output_images/all_maps_neurons_interpolation")
# Determine total neurons available based on the loaded data

# -------------------------------------------------------------------------

source("00_auxiliary_functions_DONT_RUN.R")
data <- readRDS("../../../Data/data_binary_position.RDS")
WIND <- readRDS("../../02_data_analysis/02_bSCDC_neuronal_response_to_position/output_RDS/Windows_list.RDS")
idx <- readRDS("../../02_data_analysis/02_bSCDC_neuronal_response_to_position/output_RDS/indices.RDS")
loc_neurons <- readRDS("../../../Data/M3424F_loc_neurons.RDS")
total_neurons_available <- nrow(loc_neurons)

est_cluster <- NULL
for(i in 1:length(idx)){
  n_window <- idx[i]
  filename <- paste0("../01_bSCDC_individual_trials/output_RDS/est_cluster_neurons_win", n_window, ".RDS")
  est_cl <- readRDS(file=filename)
  est_cluster <-rbind(est_cluster, est_cl)
}
rm(est_cl)

cocluster_neurons <- matrix(1,ncol(est_cluster),ncol(est_cluster))
for(i in 2:nrow(cocluster_neurons)){
  for(j in 1:(i-1)) {
    cocluster_neurons[i,j] = cocluster_neurons[j,i] = sum(apply(est_cluster[,c(i,j)], 1, function(x) x[1]==x[2]))
  }
}

# -------------------------------------------------------------------------


# add sequential window number to the calcium data
str(data)
data$pos_binary[7867:(7867+256)] = 3
data$pos_binary[(7867+256):8378] = 4
idcirc = (data$pos1^2 + data$pos2^2)<1.05
win = unlist(apply(
  cbind(
    1:length(which(diff(data$pos_binary)!=0)),
    c(0,(which(diff(data$pos_binary)!=0)))[1:138],
    (which(diff(data$pos_binary)!=0))), 1, function(x) rep(x[1], x[3]-x[2])))
win = c(win, rep(139, 10870-10862))
data$win = win
data = data[,c(1:5,331,6:330)]



# ==============================================================================
# 2. SHINY UI
# ==============================================================================
ui <- fluidPage(
  theme = bslib::bs_theme(bootswatch = "flatly"), # Cleaner theme
  
  titlePanel("Neuron Activity & Location Viewer"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Select Neurons"),
      p("Choose two indices to compare spatial firing patterns and anatomical locations."),
      
      # Numeric inputs bounded by available data
      numericInput("neu1_input", "Neuron 1 Index:", value = 1, min = 1, max = total_neurons_available, step = 1),
      numericInput("neu2_input", "Neuron 2 Index:", value = 2, min = 1, max = total_neurons_available, step = 1),
      
      hr(),
      helpText(paste("Total neurons available:", total_neurons_available)),
      helpText("Data must be loaded in the global environment.")
    ),
    
    mainPanel(
      width = 9,
      
      # Top Row: Spatial Density Maps
      fluidRow(
        column(width = 6,
               uiOutput("spatial_plot1") # <-- NEW
        ),
        # Use uiOutput to conditionally show the second plot only if neu1 != neu2
        column(width = 6,
               uiOutput("spatial_plot2")
        )
      ),
      
      hr(),
      h4("Anatomical Locations", align = "center"),
      
      # Bottom Row: Location Map centered
      fluidRow(
        column(width = 8, offset = 2,
               plotOutput("location_plot", height = "500px")
        )
      )
    )
  )
)


# ==============================================================================
# 3. SHINY SERVER
# ==============================================================================
server <- function(input, output, session) {
  
  # --- Reactive Inputs with Validation ---
  v_neu1 <- reactive({
    req(input$neu1_input)
    validate(need(input$neu1_input > 0 && input$neu1_input <= total_neurons_available, 
                  paste("Index must be between 1 and", total_neurons_available)))
    input$neu1_input
  })
  
  v_neu2 <- reactive({
    req(input$neu2_input)
    validate(need(input$neu2_input > 0 && input$neu2_input <= total_neurons_available, 
                  paste("Index must be between 1 and", total_neurons_available)))
    input$neu2_input
  })
  
  # --- Spatial Plot 1 Output ---
  output$spatial_plot1 <- renderUI({
    req(v_neu1())
    idx <- v_neu1()
    
    pdf_filename <- paste0("interpolation_neuron", idx, ".pdf")
    filesystem_path <- file.path("output_images/all_maps_neurons_interpolation", pdf_filename)
    
    if (file.exists(filesystem_path)) {
      web_path <- paste0("neuron_maps_pdf/", pdf_filename)
      tags$iframe(
        src = web_path, 
        width = "100%",
        height = "900px",
        style = "border: 1px solid gray;",
        scrolling = "auto"
      )
    } else {
      wellPanel(
        h4("PDF Not Found"),
        p(paste("Expected location:", filesystem_path))
      )
    }
  })
  
  # --- Spatial Plot 2 UI Container ---
  # Only generate the second plot structure if the indices are different
  output$spatial_plot2_container <- renderUI({
    if(v_neu1() == v_neu2()) {
      return(wellPanel(h5("Select a different second neuron to compare spatial maps.", align="center", style="color:gray;")))
    } else {
      plotOutput("spatial_plot2", height = "400px")
    }
  })
  
  # --- Spatial Plot 2 Output ---
  output$spatial_plot2 <- renderUI({
    req(v_neu2())
    idx <- v_neu2()
    
    pdf_filename <- paste0("interpolation_neuron", idx, ".pdf")
    filesystem_path <- file.path("output_images/all_maps_neurons_interpolation", pdf_filename)
    if (file.exists(filesystem_path)) {
      web_path <- paste0("neuron_maps_pdf/", pdf_filename)
      tags$iframe(
        src = web_path,
        width = "100%",
        height = "900px",
        style = "border: 1px solid gray;",
        scrolling = "auto"
      )
    } else {
      wellPanel(
        h4("PDF Not Found"),
        p(paste("Expected location:", filesystem_path))
      )
    }
  })
  
  # --- Location Plot Output ---
  output$location_plot <- renderPlot({
    plot_neuron_locations(neu1 = v_neu1(), neu2 = v_neu2())
  })
}


# ==============================================================================
# 4. RUN APP
# ==============================================================================
shinyApp(ui = ui, server = server)
