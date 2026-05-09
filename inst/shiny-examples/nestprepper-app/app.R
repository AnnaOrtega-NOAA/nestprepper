library(shiny)
library(nestprepper)
library(dplyr)
library(tidyr)
library(ggplot2)
library(shinythemes)
library(shinyBS)

ui <- fluidPage(
  theme = shinythemes::shinytheme("flatly"),
  titlePanel("nestprepper: Sea Turtle Population Dashboard"),

  sidebarLayout(
    sidebarPanel(
      h4("1. Data Entry"),
      # Allows uploading JM.csv and W.csv simultaneously
      fileInput("file1", "Upload Nesting CSV(s)", accept = ".csv", multiple = TRUE),

      uiOutput("mapper_ui"),
      uiOutput("year_filter_ui"),
      hr(),

      h4("2. Biological Parameters"),
      numericInput("clutch_freq", "Clutch Frequency:", 5.5, min = 1, max = 15, step = 0.1),
      bsTooltip("clutch_freq", "Average nests laid per female per season.", "right"),

      numericInput("remig_int", "Remigration Interval:", 3.0, min = 1, max = 10, step = 0.1),
      bsTooltip("remig_int", "Average years between nesting seasons.", "right"),

      hr(),
      h4("3. Model Settings"),
      sliderInput("iterations", "Model Precision (MCMC):",
                  min = 5000, max = 150000, value = 100000, step = 5000),

      actionButton("run_model", "Run Bayesian Model",
                   class = "btn-success", style="width: 100%; font-weight: bold; height: 50px;")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Mapping & Preview",
                 br(),
                 # Shows the raw stacked data before modeling
                 tableOutput("preview")),

        tabPanel("Results & Interpretation",
                 br(),
                 uiOutput("summary_stats"),
                 hr(),
                 plotOutput("ghost_plot", height = "550px"),
                 br(),
                 wellPanel(
                   h4("Understanding Your Results:"),
                   tags$ul(
                     tags$li(strong("Combined Trend:"), "If multiple files were uploaded, the blue line represents the total regional population."),
                     tags$li(strong("Monthly Handling:"), "If monthly columns were selected, data was aggregated into annual totals using state-space imputation."),
                     tags$li(strong("Precision:"), "Gray bands show 95% Credible Intervals (The 'Window of Truth').")
                   )
                 )
        )
      )
    )
  )
)

server <- function(input, output, session) {

  # --- MULTI-FILE LOADING ---
  raw_df <- reactive({
    req(input$file1)
    # Stacks all uploaded CSVs into one dataframe
    lapply(input$file1$datapath, read.csv, stringsAsFactors = FALSE) %>%
      dplyr::bind_rows()
  })

  # --- DYNAMIC MAPPING ---
  output$mapper_ui <- renderUI({
    cols <- names(raw_df())
    tagList(
      selectInput("year_col", "Year Column (e.g. Season):", choices = c("", cols)),

      # Detects if a 'Month' column exists to trigger aggregation logic
      if(any(grepl("month", cols, ignore.case = TRUE))) {
        selectInput("month_col", "Month Column (found monthly data):",
                    choices = c("", cols))
      },

      selectizeInput("count_cols", "Count Column(s):", choices = cols, multiple = TRUE)
    )
  })

  output$year_filter_ui <- renderUI({
    req(raw_df(), input$year_col)
    yrs <- sort(unique(as.numeric(raw_df()[[input$year_col]])))
    sliderInput("selected_years", "Analysis Year Range:",
                min = min(yrs, na.rm=T), max = max(yrs, na.rm=T),
                value = c(min(yrs, na.rm=T), max(yrs, na.rm=T)), sep = "")
  })

  # --- DATA STANDARDIZATION & AGGREGATION ---
  mapped_data <- reactive({
    req(input$year_col, input$count_cols)

    # 1. Pivot to long format
    df <- raw_df() %>%
      dplyr::rename(Year = !!input$year_col)

    # Map month if selected
    if(!is.null(input$month_col) && input$month_col != "") {
      df <- df %>% dplyr::rename(Month = !!input$month_col)
    }

    df <- df %>%
      tidyr::pivot_longer(cols = dplyr::all_of(input$count_cols),
                          names_to = "Site", values_to = "Count") %>%
      dplyr::mutate(Count = as.numeric(Count))

    # 2. RUN AGGREGATION (Handles monthly -> annual sum)
    # Ensure nestprepper::aggregate_monthly_to_annual is in your R/data_prep.R
    df <- aggregate_monthly_to_annual(df)

    return(df)
  })

  filtered_data <- reactive({
    req(mapped_data(), input$selected_years)
    mapped_data() %>%
      dplyr::filter(Year >= input$selected_years[1], Year <= input$selected_years[2])
  })

  final_results <- reactiveVal(NULL)

  observeEvent(input$run_model, {
    req(filtered_data())
    withProgress(message = 'Modeling...', detail = "Aggregating sites and running MCMC", value = 0.5, {

      # 1. Apply biology (Clutch/Remig)
      abund <- nestprepper::calculate_abundance(filtered_data(),
                                                input$clutch_freq,
                                                input$remig_int,
                                                quiet = TRUE)

      # 2. Run the JAGS model
      res <- nestprepper::run_turtle_model(abund,
                                           iter = input$iterations,
                                           parallel = FALSE)

      final_results(list(res = res, abund = abund))
    })
  })

  # --- OUTPUTS ---
  output$preview <- renderTable({ req(filtered_data()); head(filtered_data(), 20) })

  output$summary_stats <- renderUI({
    req(final_results())
    fit <- final_results()$res$fit
    years <- final_results()$res$years

    growth <- round(fit$mean$U * 100, 1)
    last_idx <- length(years)

    # Calculate Regional Estimate for the final year
    est_total <- round(sum(exp(fit$mean$X[last_idx] + fit$mean$A)))

    wellPanel(
      style = "background: #ffffff; border-left: 10px solid #27ae60;",
      h3("Integrated Population Summary"),
      h4(paste("Estimated Nesters in", years[last_idx], ":", format(est_total, big.mark=","))),
      p("Regional Growth Rate (U): ",
        span(style = ifelse(growth > 0, "color: green;", "color: red;"),
             strong(growth, "% per year")))
    )
  })

  output$ghost_plot <- renderPlot({
    req(final_results())
    nestprepper::plot_turtle_status(res = final_results()$res,
                                    abund = final_results()$abund)
  })
}

shinyApp(ui, server)
