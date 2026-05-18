library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(jagsUI)
library(bsicons)

ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "minty"),
  titlePanel("NestPrepper: NOAA Technical Memo Standard"),

  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose CSV File(s)", multiple = TRUE, accept = ".csv"),
      uiOutput("select_year"),
      uiOutput("year_range_slider"),
      uiOutput("select_month"),
      uiOutput("select_beaches"),
      hr(),
      h4("Biological Parameters"),
      numericInput("clutch_freq", "Clutch Frequency", value = 5.5, step = 0.1),
      numericInput("remig_int", "Remigration Interval", value = 3.06, step = 0.1),
      hr(),
      h4("Projection Horizon"),
      numericInput("proj_years", "Years to Project:", value = 50, min = 10, max = 100),
      actionButton("run_model", "Run Standard Analysis (150k Iterations)",
                   class = "btn-primary w-100", style = "font-weight: bold;")
    ),

    mainPanel(
      tabsetPanel(
        tabPanel("Historical Trend",
                 br(), uiOutput("summary_stats"), br(),
                 card(plotOutput("trend_plot", height = "500px"))),
        tabPanel("Future Projections",
                 br(), card(tableOutput("risk_table")),
                 card(plotOutput("projection_plot", height = "400px")))
      )
    )
  )
)

server <- function(input, output, session) {
  vault <- reactiveValues(nesters = NULL, total = NULL, trend = NULL, year = NULL, fit = NULL, data = NULL)

  # Data Logic
  df_raw <- reactive({ req(input$file1); purrr::map_df(input$file1$datapath, ~read.csv(.x)) })
  output$select_year <- renderUI({ req(df_raw()); selectInput("year_col", "Year Column", choices = names(df_raw())) })
  output$year_range_slider <- renderUI({
    req(df_raw(), input$year_col); yrs <- sort(unique(df_raw()[[input$year_col]]))
    sliderInput("year_range", "Year Range:", min = min(yrs), max = max(yrs), value = c(min(yrs), max(yrs)), sep = "")
  })
  output$select_beaches <- renderUI({ req(df_raw()); selectizeInput("beach_cols", "Beaches", choices = names(df_raw()), multiple = TRUE) })

  processed_data <- reactive({
    req(input$year_col, input$beach_cols, input$year_range)
    d <- df_raw() %>% filter(!!sym(input$year_col) >= input$year_range[1], !!sym(input$year_col) <= input$year_range[2])
    d_long <- d %>% select(all_of(c(input$year_col, input$beach_cols))) %>%
      rename(Year = !!input$year_col) %>%
      pivot_longer(cols = all_of(input$beach_cols), names_to = "Site", values_to = "Count")
    calculate_abundance(d_long, clutch_freq = input$clutch_freq, remig_int = input$remig_int)
  })

  # THE ENGINE (Defensive Version)
  observeEvent(input$run_model, {
    req(processed_data())

    withProgress(message = 'Martin et al. (2020) Workflow', value = 0, {

      # Step 1: Check if the function actually supports the new arguments
      has_args <- all(c("burnin", "thin") %in% names(formals(run_turtle_model)))

      if(!has_args){
        showNotification("Error: run_turtle_model function in models.R is outdated. Please update function arguments.", type = "error")
        return(NULL)
      }

      # Step 2: Run with Error Handling
      result <- tryCatch({
        setProgress(0.2, detail = "Sampling 150k MCMC Iterations...")
        run_turtle_model(
          processed_data(),
          iter = 150000,
          burnin = 50000,
          thin = 10
        )
      }, error = function(e) {
        showNotification(paste("JAGS Error:", e$message), type = "error")
        return(NULL)
      })

      # Step 3: Vault results if successful
      if(!is.null(result)){
        setProgress(0.8, detail = "Vaulting results...")
        final_idx <- which(result$years == max(result$years))
        vault$year <- max(result$years)
        vault$nesters <- round(exp(result$fit$mean$X[final_idx]))
        vault$total <- round(vault$nesters * input$remig_int)
        vault$trend <- round((exp(result$fit$mean$U) - 1) * 100, 2)
        vault$fit <- result
        vault$data <- processed_data()
        setProgress(1, detail = "Complete!")
      }
    })
  })

  # Renderers
  output$summary_stats <- renderUI({
    validate(need(vault$nesters, "Run analysis to see results."))
    layout_column_wrap(width = 1/3,
                       value_box(title = paste("Nesters", vault$year), value = vault$nesters, theme = "primary"),
                       value_box(title = "Total Females", value = vault$total, theme = "secondary"),
                       value_box(title = "Annual Trend", value = paste0(vault$trend, "%"), theme = if(vault$trend >=0) "success" else "danger")
    )
  })

  output$trend_plot <- renderPlot({ req(vault$fit); plot_turtle_status(vault$fit, vault$data) })

  output$risk_table <- renderTable({
    req(vault$fit); u <- vault$fit$fit$sims.list$U; h <- input$proj_years
    lapply(c(0.5, 0.25, 0.125), function(tr) {
      yrs <- log(tr) / u; yrs[yrs < 0] <- Inf
      data.frame("Status" = paste0(tr*100, "% Decline"),
                 "Risk" = paste0(round(mean(yrs <= h) * 100, 1), "%"))
    }) %>% bind_rows()
  })

  output$projection_plot <- renderPlot({
    req(vault$fit); u <- vault$fit$fit$sims.list$U
    n_final <- vault$fit$fit$mean$X[length(vault$fit$years)]
    yrs <- 0:input$proj_years
    sims <- replicate(500, exp(n_final + (sample(u, 1) * yrs)))
    matplot(yrs, sims, type = 'l', col = rgb(0,0,0,0.02), lty = 1, xlab = "Years", ylab = "Abundance")
  })
}

shinyApp(ui, server)
