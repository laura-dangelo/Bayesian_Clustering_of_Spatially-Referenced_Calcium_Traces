
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
idx <- readRDS("../../../Data/Time_windows/indices.RDS")
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
win = c(win, rep(139, nrow(data)-length(win)))
data$win = win
data = data[,c(1:5,331,6:330)]


pdf_dir <- normalizePath(
  file.path(
    "output_images",
    "all_maps_neurons_interpolation"
  ),
  winslash = "/",
  mustWork = TRUE
)

addResourcePath(
  prefix = "neuron_maps_pdf",
  directoryPath = pdf_dir
)
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
      width = 12,
      
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

  v_neu1 <- reactive({
    req(input$neu1_input)

    validate(
      need(
        input$neu1_input >= 1 &&
          input$neu1_input <= total_neurons_available,
        paste("Index must be between 1 and", total_neurons_available)
      )
    )

    as.integer(input$neu1_input)
  })

  v_neu2 <- reactive({
    req(input$neu2_input)

    validate(
      need(
        input$neu2_input >= 1 &&
          input$neu2_input <= total_neurons_available,
        paste("Index must be between 1 and", total_neurons_available)
      )
    )

    as.integer(input$neu2_input)
  })

  make_pdf_viewer <- function(neuron_index) {
    
    pdf_filename <- sprintf(
      "interpolation_neuron%d.pdf",
      neuron_index
    )
    
    filesystem_path <- file.path(
      pdf_dir,
      pdf_filename
    )
    
    if (!file.exists(filesystem_path)) {
      return(
        wellPanel(
          h4("PDF not found"),
          tags$code(filesystem_path)
        )
      )
    }
    
    web_path <- paste0(
      "/neuron_maps_pdf/",
      URLencode(pdf_filename, reserved = TRUE)
    )
    
    tagList(
      tags$object(
        data = web_path,
        type = "application/pdf",
        width = "100%",
        height = "400px",
        tags$p(
          "Inline PDF preview unavailable.",
          tags$a(
            href = web_path,
            target = "_blank",
            "Open PDF"
          )
        )
      ),
      tags$p(
        tags$a(
          href = web_path,
          target = "_blank",
          "Open PDF in a new tab"
        )
      )
    )
  }
  
  output$spatial_plot1 <- renderUI({
    make_pdf_viewer(v_neu1())
  })

  output$spatial_plot2 <- renderUI({
    if (v_neu1() == v_neu2()) {
      return(
        wellPanel(
          h5(
            "Select a different second neuron to compare spatial maps.",
            align = "center",
            style = "color: gray;"
          )
        )
      )
    }

    make_pdf_viewer(v_neu2())
  })

  output$location_plot <- renderPlot({
    plot_neuron_locations(
      neu1 = v_neu1(),
      neu2 = v_neu2(),
      loc_neurons = loc_neurons,
      shapes = c(21, 24)
    )
  })
}

# ==============================================================================
# 2. THEME
# ==============================================================================

app_theme <- bslib::bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#0f766e",
  secondary = "#64748b",
  bg = "#f3f6f9",
  fg = "#172033"
)


# ==============================================================================
# 3. SHINY UI
# ==============================================================================

ui <- bslib::page_sidebar(
  
  theme = app_theme,
  fillable = FALSE,
  window_title = "Neuron Explorer",
  
  title = tags$div(
    class = "app-title",
    
    tags$div(
      class = "title-mark",
      "N"
    ),
    
    tags$div(
      tags$div(
        class = "title-main",
        "Neuron Explorer"
      ),
      tags$div(
        class = "title-subtitle",
        "Spatial activity and anatomical comparison"
      )
    )
  ),
  
  sidebar = bslib::sidebar(
    
    width = "300px",
    open = "desktop",
    resizable = FALSE,
    bg = "#ffffff",
    padding = "1.4rem",
    
    title = tags$span(
      class = "sidebar-heading",
      "Comparison controls"
    ),
    
    tags$div(
      class = "sidebar-description",
      "Select two neurons to compare their spatial response maps and anatomical locations."
    ),
    
    tags$div(
      class = "control-label",
      tags$span(class = "neuron-dot neuron-dot-1"),
      "First neuron"
    ),
    
    numericInput(
      inputId = "neu1_input",
      label = NULL,
      value = 1,
      min = 1,
      max = total_neurons_available,
      step = 1,
      width = "100%"
    ),
    
    tags$div(
      class = "control-label second-control",
      tags$span(class = "neuron-dot neuron-dot-2"),
      "Second neuron"
    ),
    
    numericInput(
      inputId = "neu2_input",
      label = NULL,
      value = 2,
      min = 1,
      max = total_neurons_available,
      step = 1,
      width = "100%"
    ),
    
    actionButton(
      inputId = "swap_neurons",
      label = "Swap neurons",
      class = "btn-outline-primary w-100 swap-button"
    ),
    
    tags$hr(class = "sidebar-divider"),
    
    uiOutput("selected_pair"),
    
    tags$div(
      class = "dataset-summary",
      
      tags$div(
        class = "dataset-summary-label",
        "Dataset"
      ),
      
      tags$div(
        class = "dataset-summary-value",
        format(total_neurons_available, big.mark = ","),
        tags$span(" neurons")
      )
    )
  ),
  
  # ---------------------------------------------------------------------------
  # Spatial maps
  # ---------------------------------------------------------------------------
  
  tags$div(
    class = "section-heading",
    
    tags$div(
      tags$h2("Spatial activity maps"),
      tags$p(
        "Interpolated neuronal response across the experimental environment."
      )
    ),
    
    tags$span(
      class = "section-badge",
      "PDF maps"
    )
  ),
  
  bslib::layout_column_wrap(
    
    width = "450px",
    heights_equal = "row",
    fill = FALSE,
    
    bslib::card(
      class = "map-card",
      fill = FALSE,
      full_screen = TRUE,
      
      bslib::card_header(
        class = "map-card-header",
        
        tags$div(
          class = "card-heading-content",
          
          tags$div(
            class = "card-neuron-icon neuron-icon-1",
            "1"
          ),
          
          tags$div(
            tags$div(
              class = "card-title-text",
              "Spatial response"
            ),
            uiOutput("neuron_badge1")
          )
        )
      ),
      
      bslib::card_body(
        class = "p-0",
        uiOutput("spatial_plot1")
      )
    ),
    
    bslib::card(
      class = "map-card",
      fill = FALSE,
      full_screen = TRUE,
      
      bslib::card_header(
        class = "map-card-header",
        
        tags$div(
          class = "card-heading-content",
          
          tags$div(
            class = "card-neuron-icon neuron-icon-2",
            "2"
          ),
          
          tags$div(
            tags$div(
              class = "card-title-text",
              "Spatial response"
            ),
            uiOutput("neuron_badge2")
          )
        )
      ),
      
      bslib::card_body(
        class = "p-0",
        uiOutput("spatial_plot2")
      )
    )
  ),
  
  # ---------------------------------------------------------------------------
  # Anatomical locations
  # ---------------------------------------------------------------------------
  
  tags$div(
    class = "section-heading location-heading",
    
    tags$div(
      tags$h2("Anatomical locations"),
      tags$p(
        "Position of the selected neurons in the imaging field."
      )
    ),
    
    tags$span(
      class = "section-badge",
      "Location map"
    )
  ),
  
  bslib::card(
    class = "location-card",
    fill = FALSE,
    full_screen = TRUE,
    
    bslib::card_header(
      class = "location-card-header",
      
      tags$div(
        class = "location-legend",
        
        tags$span(
          class = "legend-item",
          tags$span(class = "neuron-dot neuron-dot-1"),
          uiOutput("location_label1")
        ),
        
        tags$span(
          class = "legend-item",
          tags$span(class = "neuron-dot neuron-dot-2"),
          uiOutput("location_label2")
        )
      )
    ),
    
    bslib::card_body(
      class = "location-card-body",
      
      plotOutput(
        outputId = "location_plot",
        height = "560px"
      )
    )
  ),
  
  # ---------------------------------------------------------------------------
  # CSS
  # ---------------------------------------------------------------------------
  
  tags$head(
    tags$style(
      HTML(
        "
        /* ---------------------------------------------------------------
           Global appearance
        --------------------------------------------------------------- */

        body {
          background:
            radial-gradient(
              circle at top right,
              rgba(15, 118, 110, 0.08),
              transparent 28rem
            ),
            #f3f6f9;
        }

        .bslib-page-title {
          background: #ffffff;
          border-bottom: 1px solid #e6ebef;
          padding-top: 0.85rem;
          padding-bottom: 0.85rem;
        }

        .app-title {
          display: flex;
          align-items: center;
          gap: 0.85rem;
        }

        .title-mark {
          width: 44px;
          height: 44px;
          display: flex;
          align-items: center;
          justify-content: center;
          border-radius: 13px;
          background: linear-gradient(135deg, #0f766e, #14b8a6);
          color: white;
          font-size: 1.3rem;
          font-weight: 800;
          box-shadow: 0 6px 16px rgba(15, 118, 110, 0.25);
        }

        .title-main {
          color: #172033;
          font-size: 1.25rem;
          font-weight: 750;
          letter-spacing: -0.02em;
          line-height: 1.1;
        }

        .title-subtitle {
          color: #728096;
          font-size: 0.82rem;
          font-weight: 400;
          margin-top: 0.15rem;
        }


        /* ---------------------------------------------------------------
           Sidebar
        --------------------------------------------------------------- */

        .bslib-sidebar-layout > .sidebar {
          border-right: 1px solid #e6ebef;
          box-shadow: 6px 0 24px rgba(30, 41, 59, 0.035);
        }

        .sidebar-heading {
          color: #172033;
          font-size: 1rem;
          font-weight: 750;
          letter-spacing: -0.01em;
        }

        .sidebar-description {
          color: #718096;
          font-size: 0.88rem;
          line-height: 1.55;
          margin-bottom: 1.35rem;
        }

        .control-label {
          display: flex;
          align-items: center;
          gap: 0.55rem;
          margin-bottom: 0.45rem;
          color: #354052;
          font-size: 0.84rem;
          font-weight: 700;
        }

        .second-control {
          margin-top: 1rem;
        }

        .neuron-dot {
          display: inline-block;
          width: 10px;
          height: 10px;
          flex: 0 0 10px;
          border-radius: 50%;
        }

        .neuron-dot-1 {
          background: #0f766e;
          box-shadow: 0 0 0 3px rgba(15, 118, 110, 0.13);
        }

        .neuron-dot-2 {
          background: #f59e0b;
          box-shadow: 0 0 0 3px rgba(245, 158, 11, 0.15);
        }

        .sidebar .form-control {
          min-height: 45px;
          border: 1px solid #dce3e9;
          border-radius: 10px;
          background: #f8fafc;
          font-size: 0.95rem;
          font-weight: 600;
          transition: all 0.18s ease;
        }

        .sidebar .form-control:focus {
          border-color: #14b8a6;
          background: #ffffff;
          box-shadow: 0 0 0 0.2rem rgba(20, 184, 166, 0.12);
        }

        .swap-button {
          min-height: 42px;
          margin-top: 0.65rem;
          border-radius: 10px;
          font-weight: 650;
        }

        .sidebar-divider {
          border-color: #e5eaee;
          margin: 1.4rem 0;
        }

        .selection-summary {
          padding: 0.9rem;
          border: 1px solid #e1e7eb;
          border-radius: 12px;
          background: #f8fafc;
        }

        .selection-summary-title {
          margin-bottom: 0.65rem;
          color: #788497;
          font-size: 0.72rem;
          font-weight: 750;
          letter-spacing: 0.08em;
          text-transform: uppercase;
        }

        .selected-neuron {
          display: flex;
          align-items: center;
          justify-content: space-between;
          padding: 0.35rem 0;
          color: #3c4858;
          font-size: 0.87rem;
        }

        .selected-neuron strong {
          color: #172033;
          font-size: 0.92rem;
        }

        .dataset-summary {
          margin-top: 1rem;
          padding: 0.9rem;
          border-radius: 12px;
          background: linear-gradient(135deg, #0f766e, #0d9488);
          color: white;
        }

        .dataset-summary-label {
          margin-bottom: 0.25rem;
          font-size: 0.7rem;
          font-weight: 700;
          letter-spacing: 0.09em;
          opacity: 0.75;
          text-transform: uppercase;
        }

        .dataset-summary-value {
          font-size: 1.35rem;
          font-weight: 750;
        }

        .dataset-summary-value span {
          font-size: 0.82rem;
          font-weight: 500;
          opacity: 0.82;
        }


        /* ---------------------------------------------------------------
           Section headings
        --------------------------------------------------------------- */

        .section-heading {
          display: flex;
          align-items: center;
          justify-content: space-between;
          gap: 1rem;
          margin: 0.4rem 0 1rem;
        }

        .section-heading h2 {
          margin: 0;
          color: #172033;
          font-size: 1.2rem;
          font-weight: 750;
          letter-spacing: -0.02em;
        }

        .section-heading p {
          margin: 0.2rem 0 0;
          color: #7a8798;
          font-size: 0.86rem;
        }

        .section-badge {
          padding: 0.4rem 0.7rem;
          border: 1px solid #dce5e8;
          border-radius: 999px;
          background: rgba(255, 255, 255, 0.85);
          color: #607080;
          font-size: 0.72rem;
          font-weight: 700;
          white-space: nowrap;
        }

        .location-heading {
          margin-top: 1.8rem;
        }


        /* ---------------------------------------------------------------
           Cards
        --------------------------------------------------------------- */

        .map-card,
        .location-card {
          overflow: hidden;
          border: 1px solid rgba(218, 226, 232, 0.9);
          border-radius: 16px;
          background: #ffffff;
          box-shadow:
            0 2px 3px rgba(30, 41, 59, 0.025),
            0 10px 30px rgba(30, 41, 59, 0.065);
          transition:
            transform 0.18s ease,
            box-shadow 0.18s ease;
        }

        .map-card:hover {
          transform: translateY(-2px);
          box-shadow:
            0 3px 5px rgba(30, 41, 59, 0.035),
            0 16px 38px rgba(30, 41, 59, 0.09);
        }

        .map-card-header,
        .location-card-header {
          min-height: 68px;
          display: flex;
          align-items: center;
          padding: 0.9rem 1.1rem;
          border-bottom: 1px solid #edf0f2;
          background: #ffffff;
        }

        .card-heading-content {
          display: flex;
          align-items: center;
          gap: 0.8rem;
        }

        .card-neuron-icon {
          width: 38px;
          height: 38px;
          display: flex;
          align-items: center;
          justify-content: center;
          border-radius: 11px;
          color: white;
          font-size: 0.9rem;
          font-weight: 800;
        }

        .neuron-icon-1 {
          background: linear-gradient(135deg, #0f766e, #14b8a6);
        }

        .neuron-icon-2 {
          background: linear-gradient(135deg, #d97706, #f59e0b);
        }

        .card-title-text {
          color: #263142;
          font-size: 0.92rem;
          font-weight: 750;
        }

        .neuron-badge {
          display: inline-block;
          margin-top: 0.15rem;
          color: #788597;
          font-size: 0.78rem;
          font-weight: 500;
        }


        /* ---------------------------------------------------------------
           PDF viewer
        --------------------------------------------------------------- */

        .pdf-viewer {
          background: #eef2f4;
        }

        .pdf-frame {
          display: block;
          width: 100%;
          height: 430px;
          border: 0;
          background: white;
        }

        .pdf-toolbar {
          min-height: 48px;
          display: flex;
          align-items: center;
          justify-content: space-between;
          gap: 0.8rem;
          padding: 0.55rem 0.85rem;
          border-top: 1px solid #e2e7ea;
          background: #ffffff;
        }

        .pdf-filename {
          overflow: hidden;
          color: #7b8795;
          font-size: 0.73rem;
          text-overflow: ellipsis;
          white-space: nowrap;
        }

        .pdf-open-button {
          flex: 0 0 auto;
          padding: 0.32rem 0.7rem;
          border-radius: 8px;
          font-size: 0.73rem;
          font-weight: 700;
          text-decoration: none;
        }

        .pdf-empty-state {
          min-height: 478px;
          display: flex;
          flex-direction: column;
          align-items: center;
          justify-content: center;
          padding: 2rem;
          background: #f8fafc;
          color: #6f7d8e;
          text-align: center;
        }

        .pdf-empty-state strong {
          margin-bottom: 0.4rem;
          color: #364152;
          font-size: 1rem;
        }

        .pdf-empty-state code {
          max-width: 100%;
          margin-top: 0.6rem;
          padding: 0.4rem 0.6rem;
          border-radius: 6px;
          background: #edf1f3;
          color: #506070;
          font-size: 0.72rem;
          word-break: break-all;
        }


        /* ---------------------------------------------------------------
           Anatomical plot
        --------------------------------------------------------------- */

        .location-card-body {
          padding: 0.8rem 1rem 1rem;
          background: #ffffff;
        }

        .location-legend {
          display: flex;
          align-items: center;
          gap: 1.4rem;
        }

        .legend-item {
          display: flex;
          align-items: center;
          gap: 0.55rem;
          color: #4c596a;
          font-size: 0.82rem;
          font-weight: 650;
        }


        /* ---------------------------------------------------------------
           Responsive layout
        --------------------------------------------------------------- */

        @media (max-width: 768px) {
          .title-subtitle {
            display: none;
          }

          .section-heading {
            align-items: flex-start;
          }

          .section-badge {
            display: none;
          }

          .pdf-frame {
            height: 370px;
          }

          .location-legend {
            align-items: flex-start;
            flex-direction: column;
            gap: 0.5rem;
          }
        }
        "
      )
    )
  )
)
# ==============================================================================
# 4. RUN APP
# ==============================================================================
shinyApp(ui = ui, server = server)
