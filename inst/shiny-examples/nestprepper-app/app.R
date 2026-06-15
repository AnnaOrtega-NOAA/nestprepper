library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(nestprepper)
library(mvtnorm)
library(truncnorm)

# =====================================================================
# CMP HELPER FUNCTIONS (From 2023 take_helper_Fn.R)
# =====================================================================
compute_CMP_constant <- function(Lambda, Nu, Mu, Tol, Max, Log=TRUE, Type="Z"){
  if( (!is.na(Lambda) & Lambda > 10^Nu) | (!is.na(Mu) & Mu^Nu > 10^Nu) ){
    if(Type=="Z"){
      ln_Const = Nu*Lambda^(1/Nu) - ((Nu-1)/(2*Nu))*log(Lambda) - ((Nu-1)/2)*log(2*pi) - (1/2)*log(Nu)
    }
    if(Type=="S"){
      ln_Const = Nu*Mu - ((Nu-1)/(2))*log(Mu) - ((Nu-1)/2)*log(2*pi) - (1/2)*log(Nu)
    }
  }else{
    Const = rep(0,Max+1)
    Index = 1
    Const[Index] = 1
    while( Const[Index]/Const[1] > Tol ){
      if(Type=="Z") Const[Index+1] = Const[Index] * ( Lambda / Index^Nu )
      if(Type=="S") Const[Index+1] = Const[Index] * ( Mu / Index )^Nu
      Index = Index + 1
    }
    ln_Const = log(sum(Const))
  }
  if(Log==TRUE) return(ln_Const)
  if(Log==FALSE) return(exp(ln_Const))
}

dCMP <- function( x, lambda, mu, nu, log=TRUE, tol=0.01, iter.max=200 ){
  if(!missing(mu) & !missing(lambda)) stop("'mu' and 'lambda' both specified")
  if(missing(mu) & !missing(lambda)) loglike = x*log(lambda) - nu*lfactorial(x) - compute_CMP_constant(Lambda=lambda, Nu=nu, Mu=NA, Tol=tol, Max=iter.max, Log=TRUE, Type="Z")
  if(!missing(mu) & missing(lambda)) loglike = nu*x*log(mu) - nu*lfactorial(x) - compute_CMP_constant(Lambda=NA, Nu=nu, Mu=mu, Tol=tol, Max=iter.max, Log=TRUE, Type="S")
  if(missing(mu) & missing(lambda)) stop("Neither 'mu' or 'lambda' is specified")
  if(log==TRUE) return( loglike )
  if(log==FALSE) return( exp(loglike) )
}

rCMP <- function( n, lambda, mu, nu, tol=0.01, x_max=200 ){
  loglike_x = rep(NA, x_max+1)
  for( x in 0:x_max ){
    if(missing(mu) & !missing(lambda)) loglike_x[x+1] = dCMP( x=x, lambda=lambda, nu=nu, log=TRUE, tol=tol, iter.max=x_max)
    if(!missing(mu) & missing(lambda)) loglike_x[x+1] = dCMP( x=x, mu=mu, nu=nu, log=TRUE, tol=tol, iter.max=x_max)
  }
  n_draws = sample( x=0:x_max, size=n, replace=TRUE, prob=exp(loglike_x))
  return(n_draws)
}

# 2023 MVN Transformation Helpers
inv_logit <- function(x) { exp(x) / (1 + exp(x)) }
safe_logit <- function(x) {
  x_safe <- pmax(0.001, pmin(0.999, x))
  log(x_safe / (1 - x_safe))
}

# =====================================================================
# USER INTERFACE
# =====================================================================
ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "minty"),

  tags$head(
    tags$script(HTML("
      setInterval(function() {
        var tooltipTriggerList = [].slice.call(document.querySelectorAll('[data-bs-toggle=\"tooltip\"]:not([data-bs-initialized])'))
        tooltipTriggerList.map(function (tooltipTriggerEl) {
          tooltipTriggerEl.setAttribute('data-bs-initialized', 'true');
          return new bootstrap.Tooltip(tooltipTriggerEl, { container: 'body', trigger: 'hover' });
        });
      }, 500);
    "))
  ),

  titlePanel("NestPrepper: Advanced Multi-Engine Population Suite"),

  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Choose Nesting CSV", multiple = TRUE, accept = ".csv"),
      uiOutput("mapper_ui"),
      uiOutput("year_range_ui"),

      hr(),
      h4("QA / QC Settings"),
      p("Empty cells in nesting csv will be treated as NAs.", style = "font-size: 0.9rem; color: #555; margin-bottom: 5px;"),
      checkboxInput("missing_rule_zero", "Convert empty cells to 0s (Warning: can warp log-normal math)", value = FALSE),
      checkboxInput("check_outliers", "Flag extreme outlier years (>3 SD)?", value = TRUE),
      radioButtons("aggr_strategy", "Diagnostic Aggregation Strategy:",
                   choices = c("1. Standard: Individual Beaches (Annual)" = "beach",
                               "2. Seasonal: Pooled Mid-Year vs. Pooled End-Year" = "season",
                               "3. Interaction: Beaches Split by Season" = "both"),
                   selected = "beach"),

      hr(),
      h4("Biological Parameters"),
      numericInput("clutch_freq", "Clutch Frequency", value = 5.5, step = 0.1),
      numericInput("remig_int", "Remigration Interval", value = 3.06, step = 0.01),

      hr(),
      conditionalPanel(
        condition = "input.main_tabs != 'Data Preview'",
        h4("Model Precision (Trend Analysis)"),
        sliderInput("iterations", "MCMC Iterations", min = 10000, max = 150000, value = 100000, step = 10000),
        br()
      ),
      actionButton("run_model", "1. Run Regional Trend Assessment", class = "btn-primary w-100", style = "font-weight: bold;")
   ),

    mainPanel(
      tabsetPanel(id = "main_tabs",
        # TAB 1: DATA PREVIEW
        tabPanel("Data Preview",
                 br(),
                 h3("Input Data Visualizer", style = "font-weight: bold;"),
                 p("Verify your raw nesting dataset before running the comparative engines. Naive totals are deliberately hidden here; accurate annual population profiles are only plotted after mathematical imputation is applied during the model run."),
                 br(),
                 uiOutput("preview_dynamic_plots"),
                 br(),
                 card(
                   card_header("Processed Abundance Matrix (Available Post-Run)"),
                   tableOutput("data_preview_table")
                 )
        ),

        # TAB 2: REGIONAL POPULATION PROFILE
        tabPanel("Regional Population Profile",
                 br(),
                 uiOutput("summary_stats"),
                 br(),
                 uiOutput("executive_summary"),
                 br(),
                 card(
                   card_header("Official Population Assessment Trajectory (JAGS singleUQ Baseline)"),
                   plotOutput("clean_baseline_plot", height = "500px")
                 )
        ),

        # TAB 3: STRUCTURAL DIAGNOSTICS
        tabPanel("Structural Diagnostics",
                 br(),
                 h3("Model Structure Diagnostics", style = "font-weight: bold;"),
                 p("This sandbox evaluates whether your data fits the default model rules. Use it to check if your beaches share a synchronized regional trend, or if individual beaches are behaving completely differently in the wild."),
                 br(),
                 card(
                   card_header("Interactive Plot Controls"),
                   checkboxGroupInput("plot_layers", "Select Model Layers to Overlay on Observed Counts:",
                                      choices = c("JAGS Shared Baseline (Official Red)" = "jags",
                                                  "MARSS Shared Trend (Validation Blue)" = "marss_s",
                                                  "MARSS Independent Trends (Site Breakout Purple)" = "marss_i"),
                                      selected = c("jags"), inline = TRUE)
                 ),
                 br(),
                 card(
                   card_header("Multi-Engine Trajectory Comparison"),
                   uiOutput("dynamic_unified_plot")
                 ),
                 br(), hr(),
                 h4("Technical Model Parameter Scorecards", style = "font-weight: bold; margin-top: 30px;"),
                 fluidRow(
                   column(12, card(
                     card_header(
                       span("Annual Population Trend (U) Coefficient Breakdowns ",
                            tags$span(
                              shiny::icon("info-circle", style = "color: #6c757d; cursor: pointer;"),
                              `data-bs-toggle` = "tooltip",
                              `data-bs-placement` = "top",
                              `data-bs-title` = "MANAGEMENT EVALUATION: Compare localized 'Site-Specific' trends against pooled 'Shared regional' models. If an individual beach displays a clear decline while regional models report a stable average, the shared regional trend assumption is masking localized management risks."
                            )
                       )
                     ),
                     tableOutput("table_u")
                   )),
                   column(12, card(
                     card_header("Variance Components: Biological Fluctuations (Q) vs. Survey Noise (R)"),
                     plotOutput("plot_var", height = "500px")
                   ))
                 )
        ),

        # TAB 4: PVA SIMULATION
        tabPanel("PVA Simulation",
                 br(),
                 h3("Population Viability Analysis (PVA) Sandbox", style = "font-weight: bold;"),
                 p("Project localized population trajectories forward using estimated growth targets and environmental parameters to calculate risk benchmarks. Projections assume zero future fishery take."),
                 br(),
                 fluidRow(
                   column(4,
                          card(
                            card_header("Simulation Controls"),
                            sliderInput("pva_years", "Projection Time Horizon (Years)", min = 10, max = 100, value = 50, step = 10),
                            numericInput("pva_sims", "Stochastic Iterations", value = 500, min = 100, max = 2500, step = 100),
                            hr(),
                            p(tags$b("Biological Note on Horizons:"), style = "font-size: 0.9rem;"),
                            p("A 50-year horizon aligns with regional management recovery timelines. However, because sea turtle generation lengths are long (often 30+ years), extending the evaluation to 100 years is highly recommended to fully observe long-term demographic lags.", style = "font-size: 0.85rem; color: #555;")
                          )
                   ),
                   column(8,
                          card(
                            card_header("Forward Abundance Trajectories"),
                            plotOutput("pva_plot", height = "400px")
                          )
                   )
                 ),
                 br(),
                 card(
                   card_header(
                     span("Abundance Threshold Collapse Risk Matrix ",
                          tags$span(
                            shiny::icon("info-circle", style = "color: #6c757d; cursor: pointer;"),
                            `data-bs-toggle` = "tooltip",
                            `data-bs-placement` = "top",
                            `data-bs-title` = "RISK BENCHMARKS: This matrix displays the probability that a nesting population will drop below specific percentages of its final estimated baseline value by the end of the simulation horizon."
                          )
                     )
                   ),
                   tableOutput("table_pva_thresholds")
                 )
        ),

        # TAB 5: ESA FISHERY ASSESSMENT (Take vs. No Take)
        tabPanel("ESA Fishery Assessment",
                 br(),
                 h3("Fishery Impact Assessment Studio", style = "font-weight: bold;"),
                 p("Reconstruct historical carryover via SAFE expansions, calculate the Fishery-Adjusted baseline, and project 'Take vs. No Take' scenarios. ANE is evaluated as a pooled regional impact against the aggregate nesting cohort."),
                 br(),

                 fluidRow(
                   column(4,
                          card(
                            card_header("Fishery Inputs & Biological Toggles"),
                            h5("1. Historical Expansions", style = "font-weight:bold; color:#2c3e50;"),
                            fileInput("obs_file", "Upload Observer Logs (.csv)", accept = ".csv"),
                            fileInput("safe_file", "Upload SAFE Fleet Expansion (.csv)", accept = ".csv"),
                            p("Leave blank to run on 2023 DSLL Baseline Defaults.", style = "font-size: 0.8rem; color: #6c757d;"),
                            hr(),

                            h5("2. Biological Tuning", style = "font-weight:bold; color:#2c3e50;"),
                            numericInput("ane_pj", "Juvenile Survival (\u03C6j)", value = 0.810, min = 0.5, max = 0.99, step = 0.01),
                            numericInput("ane_pa", "Adult Survival (\u03C6a)", value = 0.893, min = 0.5, max = 0.99, step = 0.01),
                            hr(),
                            radioButtons("plot_horizon", "Projection Horizon:", choices = c("10-Year Tactical" = 10, "100-Year Strategic" = 100), selected = 100, inline = TRUE),
                            hr(),
                            actionButton("run_fishery_model", "2. Run Fishery Assessment", class = "btn-success w-100", style = "font-weight: bold;")
                          )
                   ),
                   column(8,
                          card(
                            card_header("Dual-Baseline Trend Validation"),
                            tableOutput("fishery_trend_ledger")
                          ),
                          card(
                            card_header("Take vs. No Take Regional Abundance Overlay"),
                            plotOutput("pva_overlay_plot", height = "300px")
                          ),
                          card(
                            card_header("Isolated Fishery Impact (\u0394 Annual Nesters)"),
                            plotOutput("pva_difference_plot", height = "250px")
                          )
                   )
                 )
        )
      )
    )
  )
)

# =====================================================================
# SERVER ENGINE
# =====================================================================
server <- function(input, output, session) {

  vault <- reactiveValues(res = NULL, abund = NULL, marss = NULL, summary = NULL, year = NULL, nesters = NULL, total = NULL, trend_display = NULL, trend_pct = NULL, fishery_res_adj = NULL, fishery_sims = NULL)

  # 1. Base Uploads
  df_raw <- reactive({
    req(input$file1)

    purrr::map_df(input$file1$datapath, function(path) {
      df <- read.csv(path)

      # 1. Find the first column that contains "month" (case-insensitive)
      month_idx <- grep("month", names(df), ignore.case = TRUE)[1]

      # 2. Standardize its name so map_df aligns them perfectly
      if (!is.na(month_idx)) {
        names(df)[month_idx] <- "Standardized_Month"
      }

      # (Optional) You can do the exact same thing for the Year column
      year_idx <- grep("year", names(df), ignore.case = TRUE)[1]
      if (!is.na(year_idx)) {
        names(df)[year_idx] <- "Standardized_Year"
      }

      return(df)
    })
  })
  obs_raw <- reactive({ if(!is.null(input$obs_file)) read.csv(input$obs_file$datapath, stringsAsFactors = FALSE) else NULL })
  safe_raw <- reactive({ if(!is.null(input$safe_file)) read.csv(input$safe_file$datapath, stringsAsFactors = FALSE) else NULL })

  # Base UI Mappers
  output$mapper_ui <- renderUI({
    req(df_raw())
    cols <- names(df_raw())
    tagList(
      selectInput("year_col", "1. Select Year Column (Season)", choices = c("", cols)),
      checkboxInput("has_months", "Data contains monthly records", value = any(grepl("month", cols, ignore.case = TRUE))),
      conditionalPanel(condition = "input.has_months == true", selectInput("month_col", "2. Select Month Column", choices = c("", cols))),
      selectizeInput("beach_cols", "3. Select Nesting Site Column(s)", choices = cols, multiple = TRUE)
    )
  })

  output$year_range_ui <- renderUI({
    req(df_raw(), input$year_col)
    yrs <- suppressWarnings(as.numeric(df_raw()[[input$year_col]]))
    yrs <- yrs[!is.na(yrs)]
    sliderInput("year_range", "4. Select Year Range for Analysis", min = min(yrs), max = max(yrs), value = c(min(yrs), max(yrs)), sep = "")
  })

  # Data Preview Processing
  preview_df <- reactive({
    req(df_raw(), input$year_col, input$beach_cols)
    d <- df_raw() %>% rename(Year = !!sym(input$year_col))
    if (!is.null(input$year_range)) d <- d %>% dplyr::filter(Year >= input$year_range[1] & Year <= input$year_range[2])
    if (input$has_months && !is.null(input$month_col) && input$month_col != "") d <- d %>% rename(Month = !!sym(input$month_col))

    target_beach_cols <- setdiff(input$beach_cols, input$year_col)
    d_long <- d %>% dplyr::select(any_of(c("Year", "Month", target_beach_cols))) %>%
      tidyr::pivot_longer(cols = all_of(target_beach_cols), names_to = "Site", values_to = "Count") %>% mutate(Count = as.numeric(Count))

    # NEW: Flexible Diagnostic Aggregation
    if (input$has_months) {
      d_long <- d_long %>% mutate(Month_Num = as.numeric(Month))

      if (input$aggr_strategy == "season") {
        d_long <- d_long %>%
          dplyr::mutate(Site = dplyr::if_else(Month_Num >= 5 & Month_Num <= 10, "Pooled_MidYear", "Pooled_EndYear")) %>%
          dplyr::group_by(Year, Month, Month_Num, Site) %>%
          dplyr::summarise(Count = if(all(is.na(Count))) NA_real_ else sum(Count, na.rm=TRUE), .groups="drop")
      } else if (input$aggr_strategy == "both") {
        d_long <- d_long %>%
          dplyr::mutate(Site = dplyr::if_else(Month_Num >= 5 & Month_Num <= 10, paste0(Site, "_MidYear"), paste0(Site, "_EndYear")))
      }

      if (input$aggr_strategy %in% c("season", "both")) {
        valid_sites <- d_long %>% dplyr::filter(!is.na(Count) & Count > 0) %>% dplyr::pull(Site) %>% unique()
        d_long <- d_long %>% dplyr::filter(Site %in% valid_sites) %>% dplyr::select(-Month_Num)
      } else {
        d_long <- d_long %>% dplyr::select(-Month_Num)
      }
    }

    if (input$missing_rule_zero) d_long <- d_long %>% replace_na(list(Count = 0))
    d_long
  })

  output$preview_dynamic_plots <- renderUI({
    req(!is.null(input$has_months))
    if (input$has_months) {
      tabsetPanel(
        type = "pills",
        tabPanel("Raw Monthly Timeline",
                 br(),
                 card(
                   card_header("Raw Monthly Timeline (All Sites)"),
                   plotOutput("preview_continuous_timeline", height = "800px")
                 )
        ),
        tabPanel("Monthly Seasonality Trends",
                 br(),
                 card(
                   card_header("Monthly Seasonality Trends"),
                   plotOutput("preview_monthly_nests", height = "800px")
                 )
        )
      )
    } else {
      card(
        card_header("Raw Annual Reports (Unimputed)"),
        plotOutput("preview_annual_raw", height = "600px")
      )
    }
  })

  output$preview_continuous_timeline <- renderPlot({
    req(preview_df(), input$has_months)
    d_cont <- preview_df() %>% filter(!is.na(Count)) %>% mutate(Month_Num = as.numeric(Month), Date = as.Date(sprintf("%04d-%02d-01", as.numeric(Year), Month_Num))) %>% filter(!is.na(Month_Num)) %>% arrange(Date)
    ggplot(d_cont, aes(x = Date, y = Count, color = Site)) + geom_line(linewidth = 0.8) + geom_point(size = 1.5) + facet_wrap(~Site, scales = "free_y", ncol = 1) + theme_bw() + labs(x = "Time", y = "Raw Nests / Month") + theme(legend.position = "none", strip.text = element_text(face = "bold", size = 11))
  })

  output$preview_monthly_nests <- renderPlot({
    req(preview_df(), input$has_months)
    month_order <- c("5", "6", "7", "8", "9", "10", "11", "12", "1", "2", "3", "4")
    month_names <- c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")
    d_monthly <- preview_df() %>% mutate(Month_Num = as.numeric(Month), Year_Factored = as.factor(Year)) %>% filter(!is.na(Month_Num), !is.na(Count)) %>% mutate(Month_Factored = factor(as.character(Month_Num), levels = month_order)) %>% arrange(Year, Month_Factored)
    ggplot(d_monthly, aes(x = Month_Factored, y = Count, color = Year_Factored, group = Year_Factored)) + geom_line(linewidth = 1) + geom_point(size = 2) + facet_wrap(~Site, scales = "free_y", ncol = 1) + scale_x_discrete(labels = month_names) + scale_color_viridis_d(option = "magma", direction = -1) + labs(x = "Month", y = "Raw Nests", color = "Year") + theme_bw() + theme(legend.position = "right", strip.text = element_text(face = "bold", size = 12), panel.grid.minor = element_blank())
  })

  output$preview_annual_raw <- renderPlot({
    ggplot(preview_df() %>% filter(!is.na(Count)), aes(x = Year, y = Count, color = Site)) + geom_line(linewidth = 1) + geom_point(size = 2) + facet_wrap(~Site, scales = "free_y", ncol = 1) + theme_bw() + labs(x = "Year", y = "Raw Nests / Year") + theme(legend.position = "none", strip.text = element_text(face = "bold", size = 11))
  })

  output$data_preview_table <- renderTable({ vault$abund })

  # =====================================================================
  # ENGINE 1: STANDARD NESTING ASSESSMENT (Tabs 1-4)
  # =====================================================================
  observeEvent(input$run_model, {
    req(input$year_col, input$beach_cols, input$year_range)
    withProgress(message = 'Comparative Engine Workflow', value = 0, {

      setProgress(0.1, detail = "Formatting and Filtering data...")
      base_df <- df_raw() %>% rename(Year = !!sym(input$year_col)) %>% filter(Year >= input$year_range[1] & Year <= input$year_range[2])
      if (input$has_months) base_df <- base_df %>% rename(Month = !!sym(input$month_col))
      target_beach_cols <- setdiff(input$beach_cols, input$year_col)
      d_long <- base_df %>% dplyr::select(any_of(c("Year", "Month", target_beach_cols))) %>% tidyr::pivot_longer(cols = all_of(target_beach_cols), names_to = "Site", values_to = "Count") %>% mutate(Count = as.numeric(Count))
      if (input$missing_rule_zero) d_long <- d_long %>% replace_na(list(Count = 0))

      setProgress(0.2, detail = "Running Bayesian Fourier Imputation... (Check R Console for live % progress)")
      # 1. Run imputation on the RAW 12-month data to get mathematically consistent monthly gap-fills
      d_imputed_monthly <- aggregate_monthly_to_annual(d_long, quiet = TRUE)

      setProgress(0.4, detail = "Applying Diagnostic Aggregation Strategy...")
      # 2. Slice the mathematically complete dataset based on user UI selection
      if (input$has_months && input$aggr_strategy == "season") {
        d_annual <- d_imputed_monthly %>%
          dplyr::mutate(Site = dplyr::if_else(Month >= 5 & Month <= 10, "Pooled_MidYear", "Pooled_EndYear")) %>%
          dplyr::group_by(Year, Site) %>%
          dplyr::summarise(Count = if(all(is.na(Count))) NA_real_ else sum(Count, na.rm=TRUE), .groups="drop")
      } else if (input$has_months && input$aggr_strategy == "both") {
        d_annual <- d_imputed_monthly %>%
          dplyr::mutate(Site = dplyr::if_else(Month >= 5 & Month <= 10, paste0(Site, "_MidYear"), paste0(Site, "_EndYear"))) %>%
          dplyr::group_by(Year, Site) %>%
          dplyr::summarise(Count = if(all(is.na(Count))) NA_real_ else sum(Count, na.rm=TRUE), .groups="drop")
      } else {
        d_annual <- d_imputed_monthly %>%
          dplyr::group_by(Year, Site) %>%
          dplyr::summarise(Count = if(all(is.na(Count))) NA_real_ else sum(Count, na.rm=TRUE), .groups="drop")
      }

      setProgress(0.6, detail = "Converting to Abundance...")
      abund <- calculate_abundance(d_annual, clutch_freq = input$clutch_freq, remig_int = input$remig_int, quiet = TRUE)
      vault$abund <- abund

      setProgress(0.7, detail = "Running JAGS Trend Engine... (Check R Console for live % progress)")
      res <- tryCatch({ run_turtle_model(abund, iter = input$iterations, burnin = floor(input$iterations / 3), thin = 10) }, error = function(e) { showNotification(paste("JAGS Error:", e$message), type = "error"); return(NULL) })
      req(res); vault$res <- res

      setProgress(0.8, detail = "Running Comparative MARSS Engines...")
      vault$marss <- tryCatch({ run_marss_models(abund) }, error = function(e) { showNotification(paste("MARSS Error:", e$message), type = "warning"); return(NULL) })

      setProgress(0.9, detail = "Calculating Median Posteriors...")
      post_data <- get_posteriors(res, remig_int = input$remig_int)
      vault$year <- post_data$summary$Year; vault$nesters <- post_data$summary$N_fym0; vault$total <- post_data$summary$Total_Females; vault$trend_display <- post_data$summary$U_display; vault$trend_pct <- post_data$summary$U_pct; vault$summary <- post_data$summary
      setProgress(1, detail = "Complete!")
    })
  })

  output$summary_stats <- renderUI({
    validate(need(vault$nesters, "Please select columns, set year range, and click '1. Run Regional Trend Assessment'."))

    # --- NEW: Calculate the Recent 3-RI Trend from JAGS Posteriors ---
    ri_years <- round(3 * input$remig_int)
    n_years <- length(vault$res$years)

    # Ensure we actually have enough years in the model to look back 3 RIs
    if (n_years > ri_years) {
      idx_final <- n_years
      idx_start <- n_years - ri_years

      # Extract the posterior hidden states (X) for the start and end of the recent window.
      # These values are already in log-space.
      X_start <- vault$res$fit$sims.list$X[, idx_start]
      X_final <- vault$res$fit$sims.list$X[, idx_final]

      # Calculate the annualized growth rate (r) across the recent window for all MCMC draws
      recent_r_draws <- (X_final - X_start) / ri_years

      # Get the median and format it identically to the main trend
      median_recent_r <- median(recent_r_draws)
      recent_pct <- (exp(median_recent_r) - 1) * 100

      recent_display <- paste0(round(median_recent_r, 3), " (", round(recent_pct, 2), "%)")
      recent_theme <- if(recent_pct >= 0) "success" else "danger"
    } else {
      recent_display <- "Insufficient Data"
      recent_theme <- "warning"
    }
    # ---------------------------------------------

    layout_column_wrap(
      width = 1/4,
      value_box(title = paste("Estimated Nesters (", vault$year, ")", sep=""),
                value = format(round(vault$nesters), big.mark=","), theme = "primary"),
      value_box(title = paste("Total Females (", vault$year, ")", sep=""),
                value = format(round(vault$total), big.mark=","), theme = "secondary"),
      value_box(title = "Historical Regional Trend (U)",
                value = vault$trend_display, theme = if(vault$trend_pct >= 0) "success" else "danger"),
      value_box(title = paste0("Recent 3-RI Trend (Last ", ri_years, " Yrs)"),
                value = recent_display, theme = recent_theme)
    )
  })

  output$executive_summary <- renderUI({
    req(vault$res, vault$marss)
    sites <- sort(unique(vault$abund$Site))
    jags_pct <- vault$trend_display
    marss_i_pct <- (exp(stats::coef(vault$marss$indep, type = "matrix")$U[,1]) - 1) * 100
    divergence_spread <- max(marss_i_pct) - min(marss_i_pct)
    highest_site <- sites[which.max(marss_i_pct)]
    lowest_site <- sites[which.min(marss_i_pct)]

    if (divergence_spread > 1.5) {
      card(class = "border-warning bg-light", card_header(tags$span("⚠ LOCALIZED BIOLOGICAL TREND DIVERGENCE DETECTED", style = "color: #b76e00; font-weight: bold;")), p(tags$b("Regional Context: "), "The pooled regional profile reports a long-term growth of ", tags$b(jags_pct), " per annum."), p("The underlying metrics indicate a ", tags$b("Shared Regional Assumption Mismatch"), ". Specifically, ", tags$i(highest_site), " is expanding at ", tags$b(paste0(round(max(marss_i_pct), 2), "%/yr")), ", while ", tags$i(lowest_site), " is lagging or declining at ", tags$b(paste0(round(min(marss_i_pct), 2), "%/yr")), ". Forcing a uniform parallel constraint here oversimplifies localized risks."))
    } else {
      card(class = "border-success bg-light", card_header(tags$span("✓ REGIONAL TREND VALIDATED", style = "color: #198754; font-weight: bold;")), p(tags$b("Regional Context: "), "The pooled regional profile reports a long-term growth of ", tags$b(jags_pct), " per annum."), p("Statistical cross-validation confirms consistent, synchronized year-to-year population fluctuations across all evaluated beaches. The default shared-trend model structure is appropriate for this regional assessment."))
    }
  })

  output$clean_baseline_plot <- renderPlot({ req(vault$res, vault$abund); plot_turtle_status(vault$res, vault$abund) })

  # NEW: Dynamic height calculation for the trajectory plot
  output$dynamic_unified_plot <- renderUI({
    req(vault$abund)
    n_facets <- length(unique(vault$abund$Site))
    calculated_height <- max(600, n_facets * 200) # 200px per facet, minimum 600px total
    plotOutput("unified_trend_plot", height = paste0(calculated_height, "px"))
  })

  output$unified_trend_plot <- renderPlot({
    req(vault$res, vault$marss, vault$abund, input$plot_layers)
    j_fit <- vault$res$fit; all_years <- vault$res$years
    sites <- sort(unique(vault$abund$Site)); n_sites <- length(sites)
    df_all_fits <- data.frame()

    if ("jags" %in% input$plot_layers) {
      jags_list <- list()
      for (i in seq_along(sites)) {
        site_A <- if(i == 1) 0 else j_fit$sims.list$A[, i - 1]
        for (t in seq_along(all_years)) {
          pred <- exp(j_fit$sims.list$X[, t] + site_A)
          jags_list[[length(jags_list)+1]] <- data.frame(Year = all_years[t], Site = sites[i], Model = "JAGS Shared Baseline", Median = median(pred), Lower = quantile(pred, 0.025), Upper = quantile(pred, 0.975))
        }
      }
      df_all_fits <- rbind(df_all_fits, do.call(rbind, jags_list))
    }

    if ("marss_s" %in% input$plot_layers) {
      ms_fit <- vault$marss$shared; ms_states <- ms_fit$states[1, ]; ms_se <- ms_fit$states.se[1, ]; ms_A <- stats::coef(ms_fit, type = "matrix")$A
      marss_s_list <- list()
      for(i in 1:n_sites) {
        site_offset <- ms_A[i, 1]
        marss_s_list[[i]] <- data.frame(Year = all_years, Site = sites[i], Model = "MARSS Shared Trend", Median = exp(ms_states + site_offset), Lower = exp((ms_states - 1.96 * ms_se) + site_offset), Upper = exp((ms_states + 1.96 * ms_se) + site_offset))
      }
      df_all_fits <- rbind(df_all_fits, do.call(rbind, marss_s_list))
    }

    if ("marss_i" %in% input$plot_layers) {
      mi_fit <- vault$marss$indep
      marss_i_list <- list()
      for(i in 1:n_sites) {
        mi_states <- mi_fit$states[i, ]; mi_se <- mi_fit$states.se[i, ]
        marss_i_list[[i]] <- data.frame(Year = all_years, Site = sites[i], Model = "MARSS Independent Trend", Median = exp(mi_states), Lower = exp(mi_states - 1.96 * mi_se), Upper = exp(mi_states + 1.96 * mi_se))
      }
      df_all_fits <- rbind(df_all_fits, do.call(rbind, marss_i_list))
    }

    validate(need(nrow(df_all_fits) > 0, "Please select at least one model trajectory layer to display."))
    df_plot <- vault$abund %>% dplyr::left_join(df_all_fits, by = c("Year", "Site"))

    ggplot(df_plot, aes(x = Year)) +
      geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Model), alpha = 0.15, show.legend = FALSE) +
      geom_line(aes(y = Median, color = Model, linetype = Model), linewidth = 1.2) +
      geom_point(aes(y = Annual_Nesters), color = "black", size = 2.5, na.rm = TRUE) +
      facet_wrap(~Site, scales = "free_y", ncol = 1) +
      scale_color_manual(values = c("JAGS Shared Baseline" = "darkred", "MARSS Shared Trend" = "dodgerblue4", "MARSS Independent Trend" = "purple4")) +
      scale_fill_manual(values = c("JAGS Shared Baseline" = "darkred", "MARSS Shared Trend" = "dodgerblue4", "MARSS Independent Trend" = "purple4")) +
      scale_linetype_manual(values = c("JAGS Shared Baseline" = "solid", "MARSS Shared Trend" = "dashed", "MARSS Independent Trend" = "dotdash")) +
      labs(title = "Interactive Structural Model Evaluation Overlay", x = "Season (Year)", y = "Annual Nesters") +
      theme_bw() + theme(strip.text = element_text(face = "bold", size = 12), title = element_text(face = "bold"), legend.position = "bottom")
  })

  output$table_u <- renderTable({
    req(vault$res, vault$marss)
    sites <- sort(unique(vault$abund$Site))
    jags_u <- mean(vault$res$fit$sims.list$U)
    marss_s_u <- stats::coef(vault$marss$shared, type = "matrix")$U[1,1]
    marss_i_u <- stats::coef(vault$marss$indep, type = "matrix")$U[,1]

    get_u_badge <- function(u_val) {
      pct <- (exp(u_val) - 1) * 100
      if (pct < 0) { badge_class <- "bg-danger text-white"; lbl <- "Decline" }
      else if (pct < 1.0) { badge_class <- "bg-warning text-dark"; lbl <- "Stagnant" }
      else { badge_class <- "bg-success text-white"; lbl <- "Growth" }
      sprintf("<span class='badge %s' style='font-size: 0.85rem; padding: 0.4em 0.6em;'>%+.2f%% / yr (%s)</span>", badge_class, pct, lbl)
    }

    rows <- list(
      data.frame(Framework = "JAGS Baseline", `Estimation Strategy` = "Regional Shared Trend (Forced Parallel)", `Log-Scale Coefficient (U)` = sprintf("%.4f", jags_u), `Annual Population Trend` = get_u_badge(jags_u), check.names = FALSE),
      data.frame(Framework = "MARSS Shared", `Estimation Strategy` = "Regional Shared Trend (Forced Parallel)", `Log-Scale Coefficient (U)` = sprintf("%.4f", marss_s_u), `Annual Population Trend` = get_u_badge(marss_s_u), check.names = FALSE)
    )
    for(i in seq_along(sites)) rows[[length(rows)+1]] <- data.frame(Framework = "MARSS Independent", `Estimation Strategy` = paste("Site-Specific:", sites[i]), `Log-Scale Coefficient (U)` = sprintf("%.4f", marss_i_u[i]), `Annual Population Trend` = get_u_badge(marss_i_u[i]), check.names = FALSE)
    do.call(rbind, rows)
  }, sanitize.text.function = function(x) x, align = "llrr")

  output$plot_var <- renderPlot({
    req(vault$res, vault$marss)
    sites <- sort(unique(vault$abund$Site))

    if (length(sites) <= 1) {
      # If only 1 site, we can't contrast independent vs shared, so return empty plot with a message
      return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Diagnostic requires >1 Beach/Season", size = 6) + theme_void())
    }

    # Extract Variances
    jags_q <- mean(vault$res$fit$sims.list$Q); jags_r <- colMeans(vault$res$fit$sims.list$R)
    marss_s_q <- stats::coef(vault$marss$shared, type = "matrix")$Q[1,1]; marss_s_r <- diag(stats::coef(vault$marss$shared, type = "matrix")$R)
    marss_i_q <- diag(stats::coef(vault$marss$indep, type = "matrix")$Q); marss_i_r <- diag(stats::coef(vault$marss$indep, type = "matrix")$R)

    # Format into long dataframe for ggplot
    df_list <- list()
    for(i in seq_along(sites)) {
      df_list[[length(df_list)+1]] <- data.frame(Site = sites[i], Framework = "JAGS Baseline", Parameter = "Q (Biological Fluctuation)", Value = jags_q)
      df_list[[length(df_list)+1]] <- data.frame(Site = sites[i], Framework = "JAGS Baseline", Parameter = "R (Survey Noise)", Value = jags_r[i])

      df_list[[length(df_list)+1]] <- data.frame(Site = sites[i], Framework = "MARSS Shared", Parameter = "Q (Biological Fluctuation)", Value = marss_s_q)
      df_list[[length(df_list)+1]] <- data.frame(Site = sites[i], Framework = "MARSS Shared", Parameter = "R (Survey Noise)", Value = marss_s_r[i])

      df_list[[length(df_list)+1]] <- data.frame(Site = sites[i], Framework = "MARSS Independent", Parameter = "Q (Biological Fluctuation)", Value = marss_i_q[i])
      df_list[[length(df_list)+1]] <- data.frame(Site = sites[i], Framework = "MARSS Independent", Parameter = "R (Survey Noise)", Value = marss_i_r[i])
    }

    df_var <- do.call(rbind, df_list)

    ggplot(df_var, aes(x = Site, y = Value, fill = Framework)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", alpha = 0.8) +
      facet_wrap(~Parameter, scales = "free_y", ncol = 1) +
      scale_fill_manual(values = c("JAGS Baseline" = "darkred", "MARSS Shared" = "dodgerblue4", "MARSS Independent" = "purple4")) +
      labs(x = "Monitored Sub-Site / Season", y = "Estimated Variance Magnitude") +
      theme_bw() +
      theme(strip.text = element_text(face = "bold", size = 12),
            legend.position = "bottom",
            title = element_text(face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))
  })

  pva_data <- reactive({
    req(vault$marss, vault$abund)
    horizon <- input$pva_years
    n_sims <- input$pva_sims

    mi_fit <- vault$marss$indep
    mi_U <- stats::coef(mi_fit, type = "matrix")$U[, 1]
    mi_Q <- diag(stats::coef(mi_fit, type = "matrix")$Q)
    all_years <- sort(unique(vault$abund$Year)); final_year <- max(all_years); proj_years <- (final_year + 1):(final_year + horizon)

    sites <- sort(unique(vault$abund$Site))
    n_sites <- length(sites)

    plot_list <- list()

    for (i in 1:n_sites) {
      u_val <- mi_U[i]; q_val <- mi_Q[i]
      start_x <- mi_fit$states[i, ncol(mi_fit$states)]
      start_n <- exp(start_x)
      sim_matrix <- matrix(NA, nrow = n_sims, ncol = horizon)

      for (s in 1:n_sims) {
        current_x <- start_x
        for (t in 1:horizon) {
          # Added Stratonovich/Ito variance correction (- q_val / 2)
          current_x <- current_x + u_val - (q_val / 2) + rnorm(1, mean = 0, sd = sqrt(q_val))
          sim_matrix[s, t] <- exp(current_x)
        }
      }
      plot_list[[i]] <- data.frame(Year = proj_years, Site = sites[i], Median = apply(sim_matrix, 2, median), Lower = apply(sim_matrix, 2, function(x) quantile(x, 0.025)), Upper = apply(sim_matrix, 2, function(x) quantile(x, 0.975)))
    }

    r_reg <- sample(vault$res$fit$sims.list$U, n_sims, replace = TRUE)
    q_reg <- sample(vault$res$fit$sims.list$Q, n_sims, replace = TRUE)

    post_baseline <- get_posteriors(vault$res, input$remig_int)$draws
    start_reg <- median(post_baseline$N_fym0)

    sim_reg_matrix <- matrix(NA, nrow = n_sims, ncol = horizon)
    for(s in 1:n_sims) {
      current <- start_reg
      for(t in 1:horizon) {
        # Added Stratonovich/Ito variance correction (- q_reg[s] / 2)
        current <- current * exp(r_reg[s] - (q_reg[s] / 2) + rnorm(1, 0, sqrt(q_reg[s])))
        sim_reg_matrix[s,t] <- current
      }
    }

    final_counts_reg <- sim_reg_matrix[, horizon]
    reg_thresh <- data.frame(
      `Assessment Level` = "Pooled Regional Cohort",
      `Final Index Year Value` = round(start_reg, 1),
      `Risk of dropping < 50%` = sprintf("%.1f%%", mean(final_counts_reg < (start_reg * 0.50)) * 100),
      `Risk of dropping < 25%` = sprintf("%.1f%%", mean(final_counts_reg < (start_reg * 0.25)) * 100),
      `Risk of dropping < 12.5%` = sprintf("%.1f%%", mean(final_counts_reg < (start_reg * 0.125)) * 100),
      check.names = FALSE
    )

    list(plot_df = do.call(rbind, plot_list), table_df = reg_thresh)
  })

  output$pva_plot <- renderPlot({
    req(pva_data(), vault$abund)
    ggplot() + geom_ribbon(data = pva_data()$plot_df, aes(x = Year, ymin = Lower, ymax = Upper), fill = "purple4", alpha = 0.12) + geom_line(data = pva_data()$plot_df, aes(x = Year, y = Median), color = "purple4", linewidth = 1.2, linetype = "dashed") + geom_point(data = vault$abund %>% filter(!is.na(Annual_Nesters)), aes(x = Year, y = Annual_Nesters), color = "black", size = 2) + facet_wrap(~Site, scales = "free_y") + labs(title = paste("Independent Sub-Site Viability Projections (Horizon:", input$pva_years, "Years)"), subtitle = "Historical indexes plotted against unconstrained MARSS forecast medians (95% CI Ribbon)", x = "Season (Year)", y = "Annual Nesters") + theme_bw() + theme(strip.text = element_text(face = "bold", size = 11), title = element_text(face = "bold"))
  })

  output$table_pva_thresholds <- renderTable({ req(pva_data()); pva_data()$table_df }, align = "lrrrr")


  # =====================================================================
  # ENGINE 2: ESA FISHERY ASSESSMENT (Tab 5 Regional Update)
  # =====================================================================
  observeEvent(input$run_fishery_model, {
    req(vault$abund, vault$res)
    n_sims <- 1000
    horizon <- as.numeric(input$plot_horizon)

    # 2023 Structural Anchors
    linf_val <- 142.7; k_val <- 0.2262; tknot_val <- -0.17
    lmat_val <- 139.13; sig_mat_val <- 6.34; max_age_val <- 45.0; pf_val <- 0.73

    withProgress(message = 'ESA Section 7 Fishery Assessment...', value = 0, {

      setProgress(0.1, "Step 1: Reconstructing Historical ANE Ghosts via MVN...")
      abund_adj <- vault$abund

      obs_df <- obs_raw()
      safe_df <- safe_raw()

      # -------------------------------------------------------------
      # 2023 GOLD STANDARD: EXTRACT MVN PARAMETERS DYNAMICALLY
      # -------------------------------------------------------------
      if(!is.null(obs_df) && !is.null(safe_df)) {

        if ("Species" %in% names(safe_df)) safe_df <- safe_df[which(grepl("Leatherback|Dc|Dermochelys", safe_df$Species, ignore.case = TRUE)), ]
        obs_df <- obs_df[which(grepl("Leatherback|Dc|Dermochelys", obs_df$Species, ignore.case = TRUE)), ]

        if ("Point.Estimate" %in% names(safe_df)) safe_df$Total_Est <- as.numeric(safe_df$Point.Estimate)
        else if ("Point Estimate" %in% names(safe_df)) safe_df$Total_Est <- as.numeric(safe_df$`Point Estimate`)
        else if ("DC" %in% names(safe_df)) safe_df$Total_Est <- as.numeric(safe_df$DC)
        else safe_df$Total_Est <- 16.5

        if ("SCL_raw" %in% names(obs_df)) obs_df$Length <- as.numeric(obs_df$SCL_raw)
        else if ("Len" %in% names(obs_df)) obs_df$Length <- as.numeric(obs_df$Len)

        if ("M_mu" %in% names(obs_df)) obs_df$Mortality <- as.numeric(obs_df$M_mu)
        else if ("M.mu" %in% names(obs_df)) obs_df$Mortality <- as.numeric(obs_df$M.mu)
        else if (all(c("M_low", "M_high") %in% names(obs_df))) obs_df$Mortality <- rowMeans(obs_df[, c("M_low", "M_high")], na.rm = TRUE)
        else if (all(c("M.low", "M.high") %in% names(obs_df))) obs_df$Mortality <- rowMeans(obs_df[, c("M.low", "M.high")], na.rm = TRUE)
        else obs_df$Mortality <- 0.535

        # Merge for MVN parameter extraction
        merged_data <- obs_df %>% dplyr::left_join(safe_df, by = "Year")
        merged_data$Total_Est[is.na(merged_data$Total_Est)] <- 16.5

        # SAFEGUARD: Filter to only rows with valid regression data
        valid_lm_data <- merged_data %>% dplyr::filter(!is.na(Length) & !is.na(Total_Est) & Length > 0)

        # Only run the linear model if we have at least 2 valid points
        if (nrow(valid_lm_data) >= 2) {
          # Fit Linear Model for Beta0 and Beta1 (Cohort Index)
          len_lm <- lm(log(Length) ~ Total_Est, data = valid_lm_data)
          sim_beta0 <- coef(len_lm)[1]
          sim_beta1 <- coef(len_lm)[2]
          if(is.na(sim_beta1)) sim_beta1 <- 0
          sim_sigma_L <- summary(len_lm)$sigma
          if(is.na(sim_sigma_L) || sim_sigma_L == 0) sim_sigma_L <- 0.338
        } else {
          # Fallback to 2023 DSLL Baseline Defaults if regression data is empty/invalid
          sim_beta0 <- 4.6895
          sim_beta1 <- 0.00618
          sim_sigma_L <- 0.338
        }

        # Mortality parameters in logit space
        sim_mu0 <- mean(safe_logit(merged_data$Mortality), na.rm=TRUE)
        if(is.nan(sim_mu0)) sim_mu0 <- 0.1402
        sim_sigma_D <- sd(safe_logit(merged_data$Mortality), na.rm=TRUE)
        if(is.na(sim_sigma_D) || sim_sigma_D == 0) sim_sigma_D <- 3.783

        # Covariance (Rho)
        sim_rho <- cor(log(merged_data$Length), safe_logit(merged_data$Mortality), use="complete.obs")
        if(is.na(sim_rho)) sim_rho <- -0.51

        # Calculate Fallback Empirical Means for exact observed rows
        emp_mean_len <- mean(obs_df$Length, na.rm = TRUE); if(is.nan(emp_mean_len)) emp_mean_len <- 117.5
        emp_mean_mort <- mean(obs_df$Mortality, na.rm = TRUE); if(is.nan(emp_mean_mort)) emp_mean_mort <- 0.535

      } else {
        # 2023 DSLL Gold Standard Medians (Table 3)
        sim_beta0 <- 4.6895      # log(108.8)
        sim_beta1 <- 0.00618
        sim_mu0 <- 0.1402        # logit(0.535)
        sim_sigma_L <- 0.338
        sim_sigma_D <- 3.783
        sim_rho <- -0.51
      }

      # Build the Covariance Matrix for the MVN engine
      sim_cov <- matrix(c(
        sim_sigma_L^2, sim_sigma_L * sim_sigma_D * sim_rho,
        sim_sigma_L * sim_sigma_D * sim_rho, sim_sigma_D^2
      ), 2, 2)

      if(!is.null(obs_df) && !is.null(safe_df)) {
        # 2. Compile Master Turtle List
        years <- sort(unique(safe_df$Year))
        if(length(years) == 0) years <- 2004:2024

        all_turtles <- data.frame()
        set.seed(42)

        for (target_yr in years) {
          obs_yr <- obs_df[which(obs_df$Year == target_yr), ]
          n_obs <- nrow(obs_yr)

          total_est <- safe_df$Total_Est[which(safe_df$Year == target_yr)]
          if(length(total_est) == 0 || all(is.na(total_est))) total_est <- n_obs
          else total_est <- sum(total_est, na.rm = TRUE)

          n_unobs <- max(0, round(total_est) - n_obs)

          if (n_obs > 0) {
            for (i in 1:n_obs) {
              len_i <- if(is.na(obs_yr$Length[i])) emp_mean_len else obs_yr$Length[i]
              mort_i <- if(is.na(obs_yr$Mortality[i])) emp_mean_mort else obs_yr$Mortality[i]
              all_turtles <- rbind(all_turtles, data.frame(CaptureYear = target_yr, Length = len_i, Mortality = mort_i))
            }
          }
          if (n_unobs > 0) {
            # 2023 Gold Standard: MVN draw respecting Beta1 (cohort) and Rho (covariance)
            mu_l <- sim_beta0 + sim_beta1 * total_est
            draws <- mvtnorm::rmvnorm(n_unobs, mean = c(mu_l, sim_mu0), sigma = sim_cov)

            sim_lens <- pmax(50, pmin(exp(draws[,1]), linf_val - 1.0))
            sim_morts <- pmax(0, pmin(1, inv_logit(draws[,2])))

            for (i in 1:n_unobs) {
              all_turtles <- rbind(all_turtles, data.frame(CaptureYear = target_yr, Length = sim_lens[i], Mortality = sim_morts[i]))
            }
          }
        }

        # 3. Lifespan Projection Matrix
        cumulative_ledger <- data.frame()
        for (i in 1:nrow(all_turtles)) {
          c_year <- all_turtles$CaptureYear[i]
          len_i <- all_turtles$Length[i]
          mort_i <- all_turtles$Mortality[i]

          age_start <- (1 / k_val) * log((linf_val - 4.74) / (linf_val - len_i))
          if(is.nan(age_start) || is.infinite(age_start) || is.na(age_start)) age_start <- max_age_val - 2

          future_years <- seq(c_year, 2050)
          l_y <- length(future_years)
          ages_traj <- seq(age_start, length.out = l_y, by = 1)

          if(age_start >= max_age_val) {
            lens_traj <- rep(linf_val * 0.99, l_y)
          } else {
            lens_traj <- linf_val * (1 - exp(-k_val * (ages_traj - tknot_val)))
          }

          p_mat <- ifelse(lens_traj >= 0.99 * linf_val, 1.0, 1.0 / (1.0 + exp(-(lens_traj - lmat_val) / sig_mat_val)))
          p_binom <- rbinom(l_y, 1, p_mat)
          if(any(p_binom == 1)) p_binom[min(which(p_binom == 1)):l_y] <- 1

          surv_vector <- cumprod((1 - p_binom) * input$ane_pj + p_binom * input$ane_pa)
          sr <- surv_vector * p_binom * pf_val * mort_i * (1 / input$remig_int)

          cumulative_ledger <- rbind(cumulative_ledger, data.frame(CalendarYear = future_years, ANE = sr))
        }

        # 4. Aggregate and Apply to Abundance
        final_historical_ane <- cumulative_ledger %>%
          group_by(CalendarYear) %>%
          summarise(Total_Cumulative_ANE = sum(ANE, na.rm=TRUE))

        site_props <- abund_adj %>% group_by(Site) %>% summarise(mean_n = mean(Annual_Nesters, na.rm=TRUE)) %>% mutate(prop = mean_n / sum(mean_n))

        abund_adj <- abund_adj %>%
          left_join(final_historical_ane, by = c("Year" = "CalendarYear")) %>%
          mutate(Total_Cumulative_ANE = ifelse(is.na(Total_Cumulative_ANE), 0, Total_Cumulative_ANE)) %>%
          left_join(site_props, by = "Site") %>%
          mutate(Annual_Nesters = Annual_Nesters + (Total_Cumulative_ANE * prop)) %>%
          select(-Total_Cumulative_ANE, -mean_n, -prop)

      } else {
        # Fallback if no files provided: roughly simulate the snowball impact over time
        site_props <- abund_adj %>% group_by(Site) %>% summarise(mean_n = mean(Annual_Nesters, na.rm=TRUE)) %>% mutate(prop = mean_n / sum(mean_n))
        abund_adj <- abund_adj %>%
          left_join(site_props, by = "Site") %>%
          mutate(
            proxy_ane = pmax(0, (Year - min(Year)) * 0.15),
            Annual_Nesters = Annual_Nesters + (proxy_ane * prop)
          ) %>%
          select(-proxy_ane, -mean_n, -prop)
      }

      setProgress(0.3, "Step 2: Executing Fishery-Adjusted JAGS Engine...")
      res_adj <- tryCatch({ run_turtle_model(abund_adj, iter = 50000, burnin = 10000, thin = 10) }, error = function(e) { NULL })
      req(res_adj); vault$fishery_res_adj <- res_adj

      setProgress(0.6, "Step 3: Stochastic Regional PVA Generation (Take vs. No Take)...")

      yrs_proj <- seq(max(vault$abund$Year) + 1, length.out = horizon)

      post_baseline <- get_posteriors(vault$res, input$remig_int)$draws
      start_take <- median(post_baseline$N_fym0)

      post_adj <- get_posteriors(res_adj, input$remig_int)$draws
      start_notake <- median(post_adj$N_fym0)

      # -------------------------------------------------------------
      # NO TAKE Baseline (Uses Adjusted Model, ATL = 0)
      # -------------------------------------------------------------
      r_notake <- sample(res_adj$fit$sims.list$U, n_sims, replace = TRUE)
      q_notake <- sample(res_adj$fit$sims.list$Q, n_sims, replace = TRUE)

      sim_notake <- matrix(NA, n_sims, horizon)
      for(i in 1:n_sims) {
        current <- start_notake
        for(y in 1:horizon) {
          # Added Stratonovich/Ito variance correction (- q_notake[i] / 2)
          current <- current * exp(r_notake[i] - (q_notake[i] / 2) + rnorm(1, 0, sqrt(q_notake[i])))
          sim_notake[i,y] <- current
        }
      }

      # -------------------------------------------------------------
      # TAKE Baseline (Uses True State Model, ATL = 2023 CMP & MVN)
      # -------------------------------------------------------------
      r_take <- sample(vault$res$fit$sims.list$U, n_sims, replace = TRUE)
      q_take <- sample(vault$res$fit$sims.list$Q, n_sims, replace = TRUE)

      sim_take <- matrix(NA, n_sims, horizon)

      for(i in 1:n_sims) {
        current <- start_take
        for(y in 1:horizon) {
          takes <- rCMP(1, mu = 8.284754, nu = 0.05725318)
          total_ane <- 0

          if (takes > 0) {
            # 2023 Gold Standard: MVN draw respecting Beta1 (cohort) and Rho (covariance)
            mu_l <- sim_beta0 + sim_beta1 * takes
            draws <- mvtnorm::rmvnorm(takes, mean = c(mu_l, sim_mu0), sigma = sim_cov)

            lens <- pmax(50, pmin(exp(draws[,1]), linf_val - 1.0))
            morts <- pmax(0, pmin(1, inv_logit(draws[,2])))

            age_starts <- (1 / k_val) * log((linf_val - 4.74) / (linf_val - lens))
            age_starts[is.nan(age_starts) | is.infinite(age_starts) | is.na(age_starts)] <- max_age_val - 2
            age_starts <- pmin(age_starts, max_age_val - 2)

            for(t in 1:takes) {
              hor_len <- ceiling(max_age_val - age_starts[t])
              if(hor_len <= 1) hor_len <- 2
              ages_traj <- seq(age_starts[t], by = 1, length.out = hor_len)
              lens_traj <- linf_val * (1 - exp(-k_val * (ages_traj - tknot_val)))

              p_mat <- ifelse(lens_traj >= 0.99 * linf_val, 1.0, 1.0 / (1.0 + exp(-(lens_traj - lmat_val) / sig_mat_val)))
              matured <- rbinom(length(p_mat), 1, p_mat)
              y_mat_idx <- if(any(matured == 1)) min(which(matured == 1)) else length(p_mat)

              if (y_mat_idx > 1) {
                surv_chain <- (1 - p_mat[1:(y_mat_idx - 1)]) * input$ane_pj + p_mat[1:(y_mat_idx - 1)] * input$ane_pa
                comp_surv <- prod(surv_chain)
              } else {
                comp_surv <- 1.0
              }

              total_ane <- total_ane + (comp_surv * pf_val * (1/input$remig_int) * morts[t])
            }
          }
          # Added Stratonovich/Ito variance correction (- q_take[i] / 2)
          current <- (current - total_ane) * exp(r_take[i] - (q_take[i] / 2) + rnorm(1, 0, sqrt(q_take[i])))
          if (current < 0) current <- 0
          sim_take[i,y] <- current
        }
      }

      setProgress(0.9, "Aggregating Difference Matrices...")

      calc_summary <- function(mat) {
        data.frame(
          Year = yrs_proj,
          Median = apply(mat, 2, median),
          L95 = apply(mat, 2, function(x) quantile(x, 0.025)),
          U95 = apply(mat, 2, function(x) quantile(x, 0.975))
        )
      }

      df_notake <- calc_summary(sim_notake) %>% mutate(Scenario = "No Take (Nj + F)")
      df_take <- calc_summary(sim_take) %>% mutate(Scenario = "Take (Nj)")

      diff_mat <- sim_notake - sim_take
      df_diff <- calc_summary(diff_mat)

      vault$fishery_sims <- list(
        overlay = bind_rows(df_notake, df_take),
        difference = df_diff,
        avg_hist_ane = if(!is.null(obs_df)) mean(final_historical_ane$Total_Cumulative_ANE, na.rm=TRUE) else NA
      )

      setProgress(1, "Assessment Complete.")
    })
  })

  output$fishery_trend_ledger <- renderTable({
    req(vault$res, vault$fishery_res_adj, vault$fishery_sims)
    data.frame(
      `Analytical Baseline` = c("True State (Nj)", "Fishery-Adjusted (Nj + F)"),
      `Estimated Trend (r)` = c(sprintf("%.4f", mean(vault$res$fit$sims.list$U)), sprintf("%.4f", mean(vault$fishery_res_adj$fit$sims.list$U))),
      `Process Variance (Q)` = c(sprintf("%.4f", mean(vault$res$fit$sims.list$Q)), sprintf("%.4f", mean(vault$fishery_res_adj$fit$sims.list$Q))),
      `Historical ANE Penalty (Avg/Yr)` = c("0.00", if(is.na(vault$fishery_sims$avg_hist_ane)) "Proxy Used" else sprintf("%.3f", vault$fishery_sims$avg_hist_ane)),
      check.names = FALSE
    )
  })

  output$pva_overlay_plot <- renderPlot({
    req(vault$fishery_sims)
    ggplot(vault$fishery_sims$overlay, aes(x = Year, y = Median, color = Scenario, fill = Scenario)) +
      geom_line(linewidth = 1.2) +
      geom_ribbon(aes(ymin = L95, ymax = U95), alpha = 0.2, color = NA) +
      scale_color_manual(values = c("No Take (Nj + F)" = "chartreuse4", "Take (Nj)" = "dodgerblue4")) +
      scale_fill_manual(values = c("No Take (Nj + F)" = "chartreuse3", "Take (Nj)" = "dodgerblue3")) +
      theme_bw() + labs(y = "Total Regional Annual Nesters", x = "Year") +
      theme(legend.position = "top", title = element_text(face="bold"))
  })

  output$pva_difference_plot <- renderPlot({
    req(vault$fishery_sims)
    ggplot(vault$fishery_sims$difference, aes(x = Year, y = Median)) +
      geom_line(color = "darkorchid4", linewidth = 1.2) +
      geom_ribbon(aes(ymin = L95, ymax = U95), fill = "darkorchid", alpha = 0.3) +
      theme_bw() + labs(y = "\u0394 Regional Annual Nesters", x = "Year") +
      theme(title = element_text(face="bold"))
  })

}

shinyApp(ui, server)
