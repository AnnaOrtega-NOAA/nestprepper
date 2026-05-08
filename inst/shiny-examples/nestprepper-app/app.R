library(shiny)
library(nestprepper)
library(dplyr)
library(tidyr)
library(ggplot2)
library(shinythemes)
library(zoo)
library(shinyBS)

ui <- fluidPage(
  theme = shinythemes::shinytheme("flatly"),
  titlePanel("nestprepper: Sea Turtle Population Dashboard"),

  sidebarLayout(
    sidebarPanel(
      h4("1. Data Entry"),
      fileInput("file1", "Upload Nesting CSV", accept = ".csv"),
      uiOutput("mapper_ui"),

      hr(),
      h4("3. Biological Parameters"),

      # Attach directly to 'clutch_freq'
      numericInput("clutch_freq", "Clutch Frequency:", 5.5, min = 1, max = 15, step = 0.1),
      bsTooltip("clutch_freq", "Average number of nests one female turtle lays in a single season.",
                "right", options = list(container = "body")),

      # Attach directly to 'remig_int'
      numericInput("remig_int", "Remigration Interval:", 3.0, min = 1, max = 10, step = 0.1),
      bsTooltip("remig_int", "Average number of years a female turtle takes between nesting seasons.",
                "right", options = list(container = "body")),

      hr(),
      h4("4. Model Settings"),

      # Attach directly to 'iterations'
      sliderInput("iterations", "Model Precision (MCMC):", 500, 20000, 2000, step = 500),
      bsTooltip("iterations", "Higher values lead to more stable results but take longer to run.",
                "right", options = list(container = "body")),

      actionButton("run_model", "Run Bayesian Model",
                   class = "btn-success", style="width: 100%; font-weight: bold; height: 50px;")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Mapping & Preview",
                 br(),
                 uiOutput("outlier_alert"),
                 tableOutput("outlier_table"),
                 hr(),
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
                     tags$li(strong("The Estimate:"), "Our best guess for the population size based on the data provided."),
                     tags$li(strong("Confidence Range:"), "The window where the true population likely falls."),
                     tags$li(strong("Trend Line:"), "Upward means recovery; downward means decline.")
                   )
                 ))
      )
    )
  )
)

server <- function(input, output, session) {

  # --- 1. DATA LOADING ---
  raw_df <- reactive({
    req(input$file1)
    read.csv(input$file1$datapath, stringsAsFactors = FALSE)
  })

  # --- 2. DYNAMIC MAPPING ---
  output$mapper_ui <- renderUI({
    cols <- names(raw_df())
    tagList(
      h4("2. Map Your Columns"),
      selectInput("year_col", "Year Column:", choices = c("", cols)),
      selectizeInput("count_cols", "Count Column(s):",
                     choices = cols, multiple = TRUE,
                     options = list(placeholder = 'Select one or more columns'))
    )
  })

  # --- 3. AUTO-DATA PROCESSING ---
  mapped_data <- reactive({
    req(input$year_col, input$count_cols)

    # 1. Start with the raw data and rename the Year
    df <- raw_df() %>%
      dplyr::rename(Year = !!input$year_col)

    # 2. Identify if we have multiple count columns
    selected_counts <- input$count_cols

    # 3. PIVOT: This is the key. It takes all selected columns and stacks them.
    # If you select Beach A and Beach B, it creates a 'Site' column with those names.
    df_long <- df %>%
      tidyr::pivot_longer(
        cols = dplyr::all_of(selected_counts),
        names_to = "Site",
        values_to = "Count"
      )

    # 4. Clean up: ensure numbers are numeric and filter NAs in Year
    df_long %>%
      dplyr::mutate(Count = as.numeric(Count)) %>%
      dplyr::select(Year, Site, Count) %>%
      dplyr::filter(!is.na(Year))
  })

  # --- 4. QAQC ---
  outliers <- reactive({
    req(mapped_data())
    mapped_data() %>% group_by(Site) %>%
      filter(Count > (mean(Count, na.rm=TRUE) + 3*sd(Count, na.rm=TRUE)))
  })

  output$outlier_alert <- renderUI({
    req(nrow(outliers()) > 0)
    div(class = "alert alert-warning", icon("triangle-exclamation"),
        "Biological Spikes Detected: Review high counts below for potential data entry errors.")
  })

  output$outlier_table <- renderTable({ req(nrow(outliers()) > 0); outliers() })

  # --- 5. EXECUTION ---
  final_results <- reactiveVal(NULL)

  observeEvent(input$run_model, {
    req(mapped_data())
    if (any(is.na(mapped_data()$Count))) {
      showModal(modalDialog(
        title = "Data Gaps Detected",
        "How should the model handle missing years?",
        footer = tagList(
          actionButton("na_zero", "Assume 0 nests", class = "btn-info"),
          actionButton("na_interp", "Fill Gaps (Interpolate)", class = "btn-info"),
          modalButton("Cancel")
        )
      ))
    } else { execute_model("none") }
  })

  observeEvent(input$na_zero, { removeModal(); execute_model("zero") })
  observeEvent(input$na_interp, { removeModal(); execute_model("interpolate") })

  execute_model <- function(na_handle) {
    withProgress(message = 'Estimating Population...', value = 0.2, {
      dat <- mapped_data()
      if(na_handle == "zero") dat$Count[is.na(dat$Count)] <- 0
      if(na_handle == "interpolate") {
        dat <- dat %>% group_by(Site) %>%
          mutate(Count = round(zoo::na.approx(Count, na.rm = FALSE))) %>% ungroup()
      }

      results <- nestprepper::run_nestprepper_workflow(
        df = dat, clutch_freq = input$clutch_freq,
        remig_int = input$remig_int, iter = input$iterations, quiet = TRUE
      )
      final_results(results)
    })
  }

  # --- 6. SUMMARY UI ---
  output$summary_stats <- renderUI({
    req(final_results())

    fit <- final_results()$res$fit
    years <- final_results()$res$years

    est_log <- fit$summary[grep("X", rownames(fit$summary)), ]
    est_nat <- exp(est_log[, c("2.5%", "50%", "97.5%")])

    last_idx <- nrow(est_nat)
    final_year <- years[last_idx]
    final_est  <- round(est_nat[last_idx, "50%"])
    low_ci     <- round(est_nat[last_idx, "2.5%"])
    high_ci    <- round(est_nat[last_idx, "97.5%"])
    growth     <- round(fit$mean$U * 100, 1)

    wellPanel(
      style = "background: #ffffff; border-left: 10px solid #27ae60; box-shadow: 0 4px 6px rgba(0,0,0,0.1);",
      h3(paste("Snapshot for", final_year)),
      p(style = "font-size: 20px; color: #2c3e50;",
        "Estimated nesting females: ", strong(final_est)),
      p("Confidence Range: ", strong(low_ci), " to ", strong(high_ci), "."),
      p(icon("chart-line"), "Annual growth: ",
        span(style = ifelse(growth > 0, "color: green;", "color: red;"), strong(growth, "%")))
    )
  })

  output$preview <- renderTable({ req(mapped_data()); head(mapped_data(), 10) })

  output$ghost_plot <- renderPlot({
    req(final_results())
    plot_turtle_status(final_results()$res, final_results()$abund, "Population Trend Estimate")
  })
}

shinyApp(ui, server)
