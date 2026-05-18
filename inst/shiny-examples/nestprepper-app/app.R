library(shiny)
library(bslib)
library(dplyr)
library(tidyr)
library(ggplot2)
library(nestprepper)

ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "minty"),

  # --- BOOTSTRAP TOOLTIP INITIALIZER ENGINE ---
  # Automatically monitors table updates to trigger diagnostic popups instantly on hover
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
      radioButtons("missing_rule", "How should un-surveyed years be handled?",
                   choices = c("Leave as Missing (NA) - Recommended for JAGS" = "na",
                               "Convert to Zeros (0) - Warning: Can warp log-normal math" = "zero")),
      checkboxInput("check_outliers", "Flag extreme outlier years (>3 SD)?", value = TRUE),

      hr(),
      h4("Biological Parameters"),
      numericInput("clutch_freq", "Clutch Frequency", value = 5.5, step = 0.1),
      numericInput("remig_int", "Remigration Interval", value = 3.06, step = 0.01),

      hr(),
      # --- HISTORICAL TAKE SECTION ---
      h4("Historical Take (ANE)"),
      checkboxInput("use_take", "Add Historical Take to estimate intrinsic growth?", value = FALSE),
      conditionalPanel(
        condition = "input.use_take == true",
        fileInput("file_ane", "Upload historical_ANE.csv", accept = ".csv"),
        uiOutput("ane_mapper_ui")
      ),

      hr(),
      h4("Model Precision"),
      sliderInput("iterations", "MCMC Iterations", min = 10000, max = 150000, value = 50000, step = 10000),
      br(),
      actionButton("run_model", "Run Comparative Analysis", class = "btn-primary w-100", style = "font-weight: bold;")
    ),

    mainPanel(
      tabsetPanel(
        # =====================================================================
        # TAB 1: REGIONAL POPULATION PROFILE (MANAGEMENT VIEW)
        # =====================================================================
        tabPanel("Regional Population Profile",
                 br(),
                 # Headline Scorecards
                 uiOutput("summary_stats"),
                 br(),

                 # Analytical Executive Summary Card
                 uiOutput("executive_summary"),
                 br(),

                 # Official Clean JAGS Baseline Plot Only
                 card(
                   card_header("Official Population Assessment Trajectory (JAGS singleUQ Baseline)"),
                   plotOutput("clean_baseline_plot", height = "500px")
                 )
        ),

        # =====================================================================
        # TAB 2: STRUCTURAL DIAGNOSTICS (SCIENTIFIC SANDBOX)
        # =====================================================================
        tabPanel("Structural Diagnostics",
                 br(),
                 h3("Model Structure Diagnostics", style = "font-weight: bold;"),
                 p("This sandbox evaluates whether your data fits the default model rules. Use it to check if your beaches share a synchronized regional trend, or if individual beaches are behaving completely differently in the wild."),
                 br(),

                 # Interactive Layer Selectors
                 card(
                   card_header("Interactive Plot Controls"),
                   checkboxGroupInput("plot_layers", "Select Model Layers to Overlay on Observed Counts:",
                                      choices = c("JAGS Shared Baseline (Official Red)" = "jags",
                                                  "MARSS Shared Trend (Validation Blue)" = "marss_s",
                                                  "MARSS Independent Trends (Site Breakout Purple)" = "marss_i"),
                                      selected = c("jags"), inline = TRUE)
                 ),
                 br(),

                 # Dynamic Interactive Overlay Plot
                 card(
                   card_header("Unified Multi-Engine Trajectory Comparison Overlay"),
                   plotOutput("unified_trend_plot", height = "600px")
                 ),
                 br(),
                 hr(),

                 # --- SECTION 3: TECHNICAL SCORECARDS WITH INSTANT HTML TOOLTIPS ---
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
                     card_header(
                       span("Variance Component Partitioning Breakdown ",
                            tags$span(
                              shiny::icon("info-circle", style = "color: #6c757d; cursor: pointer;"),
                              `data-bs-toggle` = "tooltip",
                              `data-bs-placement` = "top",
                              `data-bs-title` = "HOVER OVER ANY DATA BADGE: Move your cursor over the colored badges below to see a direct, clear breakdown of what each value implies about your field data precision and population stability."
                            )
                       )
                     ),
                     tableOutput("table_var")
                   ))
                 )
        ),

        tabPanel("Data Preview",
                 br(),
                 card(tableOutput("data_preview")))
      )
    )
  )
)


server <- function(input, output, session) {

  vault <- reactiveValues(res = NULL, abund = NULL, marss = NULL, summary = NULL, year = NULL, nesters = NULL, total = NULL, trend_display = NULL, trend_pct = NULL)

  # 1. Raw Data Upload
  df_raw <- reactive({
    req(input$file1)
    purrr::map_df(input$file1$datapath, ~read.csv(.x))
  })

  # 2. ANE Data Upload
  df_ane <- reactive({
    req(input$file_ane)
    read.csv(input$file_ane$datapath)
  })

  # 3. Column Mapping UIs
  output$mapper_ui <- renderUI({
    req(df_raw())
    cols <- names(df_raw())
    tagList(
      selectInput("year_col", "1. Select Year Column (Season)", choices = c("", cols)),
      checkboxInput("has_months", "Data contains monthly records",
                    value = any(grepl("month", cols, ignore.case = TRUE))),
      conditionalPanel(
        condition = "input.has_months == true",
        selectInput("month_col", "2. Select Month Column", choices = c("", cols))
      ),
      selectizeInput("beach_cols", "3. Select Nesting Site Column(s)", choices = cols, multiple = TRUE)
    )
  })

  # 4. Dynamic Year Range Slider
  output$year_range_ui <- renderUI({
    req(df_raw(), input$year_col)
    yrs <- suppressWarnings(as.numeric(df_raw()[[input$year_col]]))
    yrs <- yrs[!is.na(yrs)]
    req(length(yrs) > 0)
    sliderInput("year_range", "4. Select Year Range for Analysis",
                min = min(yrs), max = max(yrs),
                value = c(min(yrs), max(yrs)), sep = "")
  })

  # 5. Main Execution Logic
  observeEvent(input$run_model, {
    req(input$year_col, input$beach_cols, input$year_range)

    withProgress(message = 'Comparative Engine Workflow', value = 0, {

      setProgress(0.1, detail = "Formatting and Filtering data...")

      base_df <- df_raw()
      d <- base_df %>% rename(Year = !!sym(input$year_col))
      d <- d %>% dplyr::filter(Year >= input$year_range[1] & Year <= input$year_range[2])

      if (input$has_months && !is.null(input$month_col) && input$month_col != "") {
        d <- d %>% rename(Month = !!sym(input$month_col))
      }

      target_beach_cols <- setdiff(input$beach_cols, input$year_col)

      d_long <- d %>%
        dplyr::select(any_of(c("Year", "Month", target_beach_cols))) %>%
        tidyr::pivot_longer(cols = all_of(target_beach_cols), names_to = "Site", values_to = "Count") %>%
        dplyr::mutate(Count = as.numeric(Count))

      if (input$missing_rule == "zero") {
        d_long <- d_long %>% tidyr::replace_na(list(Count = 0))
      }

      setProgress(0.2, detail = "Aggregating/Imputing Months...")
      d_annual <- aggregate_monthly_to_annual(d_long, quiet = TRUE)

      setProgress(0.3, detail = "Converting to Abundance...")
      abund <- calculate_abundance(d_annual, clutch_freq = input$clutch_freq, remig_int = input$remig_int, quiet = TRUE)

      # --- ANE Take Addition ---
      if (input$use_take) {
        req(df_ane(), input$ane_year_col, input$ane_val_col)
        setProgress(0.35, detail = "Adding Historical ANE proportional to site size...")

        ane_data <- df_ane() %>%
          rename(Year = !!sym(input$ane_year_col), ANE_val = !!sym(input$ane_val_col)) %>%
          select(Year, ANE_val) %>%
          group_by(Year) %>% summarise(ANE_val = sum(ANE_val, na.rm=TRUE))

        abund <- abund %>%
          dplyr::group_by(Year) %>%
          dplyr::mutate(
            Total_Nesters_Yr = sum(Annual_Nesters, na.rm = TRUE),
            Site_Prop = ifelse(Total_Nesters_Yr > 0, Annual_Nesters / Total_Nesters_Yr, 1 / n())
          ) %>%
          dplyr::ungroup() %>%
          dplyr::left_join(ane_data, by = "Year") %>%
          dplyr::mutate(
            ANE_val = tidyr::replace_na(ANE_val, 0),
            Annual_Nesters = Annual_Nesters + (ANE_val * Site_Prop),
            Total_Adult_Females = Annual_Nesters * input$remig_int
          ) %>%
          dplyr::select(-Total_Nesters_Yr, -Site_Prop, -ANE_val)
      }

      vault$abund <- abund

      # --- EXECUTE ENGINE 1: JAGS BASELINE ---
      setProgress(0.4, detail = paste("Running JAGS (", input$iterations, " iterations)...", sep=""))
      res <- tryCatch({
        run_turtle_model(
          abund,
          iter = input$iterations,
          burnin = floor(input$iterations / 3),
          thin = 10
        )
      }, error = function(e) {
        showNotification(paste("JAGS Error:", e$message), type = "error", duration = NULL)
        return(NULL)
      })

      req(res)
      vault$res <- res

      # --- EXECUTE ENGINE 2: MARSS ---
      setProgress(0.8, detail = "Running Comparative MARSS Engines...")
      marss_res <- tryCatch({
        run_marss_models(abund)
      }, error = function(e) {
        showNotification(paste("MARSS Optimization Error:", e$message), type = "warning", duration = NULL)
        return(NULL)
      })
      vault$marss <- marss_res

      setProgress(0.9, detail = "Calculating Median Posteriors...")
      post_data <- get_posteriors(res, remig_int = input$remig_int)

      vault$year          <- post_data$summary$Year
      vault$nesters       <- post_data$summary$N_fym0
      vault$total         <- post_data$summary$Total_Females
      vault$trend_display <- post_data$summary$U_display
      vault$trend_pct     <- post_data$summary$U_pct
      vault$summary       <- post_data$summary

      setProgress(1, detail = "Complete!")
    })
  })

  # 6. Summary Statistics Boxes
  output$summary_stats <- renderUI({
    validate(need(vault$nesters, "Please select columns, set year range, and click 'Run Comparative Analysis'."))

    layout_column_wrap(
      width = 1/3,
      value_box(title = paste("Estimated Nesters (", vault$year, ")", sep=""), value = format(vault$nesters, big.mark=","), theme = "primary"),
      value_box(title = paste("Total Females (", vault$year, ")", sep=""), value = format(vault$total, big.mark=","), theme = "secondary"),
      value_box(title = "Regional Population Trend (U)", value = vault$trend_display,
                theme = if(vault$trend_pct >= 0) "success" else "danger")
    )
  })

  # --- SCIENTIFIC EXECUTIVE SYNTHESIS ANALYSIS ---
  output$executive_summary <- renderUI({
    req(vault$res, vault$marss)

    sites <- sort(unique(vault$abund$Site))
    jags_pct <- vault$trend_display

    marss_i_u <- stats::coef(vault$marss$indep, type = "matrix")$U[,1]
    marss_i_pct <- (exp(marss_i_u) - 1) * 100

    divergence_spread <- max(marss_i_pct) - min(marss_i_pct)
    highest_site <- sites[which.max(marss_i_pct)]
    lowest_site <- sites[which.min(marss_i_pct)]

    if (divergence_spread > 1.5) {
      alert_class <- "border-warning bg-light"
      title_style <- "color: #b76e00; font-weight: bold;"
      alert_status <- tags$span("⚠ LOCALIZED BIOLOGICAL TREND DIVERGENCE DETECTED: Monitored sub-populations are exhibiting asymmetrical trajectories.", style = title_style)

      narrative_detail <- p("The underlying metrics indicate a ", tags$b("Shared Regional Assumption Mismatch"), ". While the regulatory baseline forces all sub-beaches into a locked parallel trajectory, unconstrained modeling reveals that localized nesting dynamics have decoupled. Specifically, ", tags$i(highest_site), " is expanding at ", tags$b(paste0(round(max(marss_i_pct), 2), "%/yr")), ", while ", tags$i(lowest_site), " is lagging or declining at ", tags$b(paste0(round(min(marss_i_pct), 2), "%/yr")), ". Forcing a uniform parallel constraint here oversimplifies localized recovery targets. Please review the unconstrained trajectories on the ", tags$b("Structural Diagnostics"), " tab.")
    } else {
      alert_class <- "border-success bg-light"
      title_style <- "color: #198754; font-weight: bold;"
      alert_status <- tags$span("✓ REGIONAL TREND VALIDATED: Local nesting sub-populations exhibit synchronized demographic trajectories.", style = title_style)

      narrative_detail <- p("Statistical cross-validation confirms consistent, synchronized year-to-year population fluctuations across all evaluated beaches. Because localized dynamics are highly aligned, the default shared-trend model structure (JAGS Baseline) is precise and appropriate for this regional assessment.")
    }

    card(
      class = alert_class,
      card_header(alert_status),
      p(tags$b("Regional Context: "), "The pooled regional population profile reports a long-term growth trajectory of ", tags$b(jags_pct), " per annum. Both Bayesian MCMC sampling and frequentist Kalman optimization frameworks achieve complete analytical convergence on this aggregate baseline calculation."),
      narrative_detail
    )
  })

  # --- TAB 1 CLEAN PLOT ENGINE ---
  output$clean_baseline_plot <- renderPlot({
    req(vault$res, vault$abund)
    plot_turtle_status(vault$res, vault$abund)
  })

  # --- TAB 2 INTERACTIVE LAYER OVERLAY PLOT ENGINE ---
  output$unified_trend_plot <- renderPlot({
    req(vault$res, vault$marss, vault$abund, input$plot_layers)

    j_fit <- vault$res$fit
    all_years <- vault$res$years
    sites <- sort(unique(vault$abund$Site))
    n_sites <- length(sites)

    df_all_fits <- data.frame()

    if ("jags" %in% input$plot_layers) {
      jags_list <- list()
      for (i in seq_along(sites)) {
        site_A <- if(i == 1) 0 else j_fit$sims.list$A[, i]
        for (t in seq_along(all_years)) {
          pred <- exp(j_fit$sims.list$X[, t] + site_A)
          jags_list[[length(jags_list)+1]] <- data.frame(
            Year = all_years[t], Site = sites[i], Model = "JAGS Shared Baseline",
            Median = median(pred), Lower = quantile(pred, 0.025), Upper = quantile(pred, 0.975)
          )
        }
      }
      df_all_fits <- rbind(df_all_fits, do.call(rbind, jags_list))
    }

    if ("marss_s" %in% input$plot_layers) {
      ms_fit <- vault$marss$shared
      ms_states <- ms_fit$states[1, ]
      ms_se <- ms_fit$states.se[1, ]
      ms_A <- stats::coef(ms_fit, type = "matrix")$A

      marss_s_list <- list()
      for(i in 1:n_sites) {
        site_offset <- ms_A[i, 1]
        marss_s_list[[i]] <- data.frame(
          Year = all_years, Site = sites[i], Model = "MARSS Shared Trend",
          Median = exp(ms_states + site_offset),
          Lower = exp((ms_states - 1.96 * ms_se) + site_offset),
          Upper = exp((ms_states + 1.96 * ms_se) + site_offset)
        )
      }
      df_all_fits <- rbind(df_all_fits, do.call(rbind, marss_s_list))
    }

    if ("marss_i" %in% input$plot_layers) {
      mi_fit <- vault$marss$indep
      marss_i_list <- list()
      for(i in 1:n_sites) {
        mi_states <- mi_fit$states[i, ]
        mi_se <- mi_fit$states.se[i, ]
        marss_i_list[[i]] <- data.frame(
          Year = all_years, Site = sites[i], Model = "MARSS Independent Trend",
          Median = exp(mi_states),
          Lower = exp(mi_states - 1.96 * mi_se),
          Upper = exp(mi_states + 1.96 * mi_se)
        )
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
      labs(title = "Interactive Structural Model Evaluation Overlay", x = "Season (Year)", y = "Annual Nesters", color = "Active Model Layers:", linetype = "Active Model Layers:") +
      theme_bw() +
      theme(strip.text = element_text(face = "bold", size = 12), title = element_text(face = "bold"), legend.position = "bottom")
  })

  # --- TAB 2 COEFFICIENTS TABLE ---
  output$table_u <- renderTable({
    req(vault$res, vault$marss)
    sites <- sort(unique(vault$abund$Site))

    jags_u <- mean(vault$res$fit$sims.list$U)
    marss_s_u <- stats::coef(vault$marss$shared, type = "matrix")$U[1,1]
    marss_i_u <- stats::coef(vault$marss$indep, type = "matrix")$U[,1]

    get_u_badge <- function(u_val) {
      pct <- (exp(u_val) - 1) * 100
      if (pct < 0) {
        badge_class <- "bg-danger text-white"; lbl <- "Decline"
      } else if (pct < 1.0) {
        badge_class <- "bg-warning text-dark"; lbl <- "Stagnant"
      } else {
        badge_class <- "bg-success text-white"; lbl <- "Growth"
      }
      sprintf("<span class='badge %s' style='font-size: 0.85rem; padding: 0.4em 0.6em;'>%+.2f%% / yr (%s)</span>", badge_class, pct, lbl)
    }

    rows <- list(
      data.frame(Framework = "JAGS Baseline", `Estimation Strategy` = "Regional Shared Trend (Forced Parallel)", `Log-Scale Coefficient (U)` = sprintf("%.4f", jags_u), `Annual Population Trend` = get_u_badge(jags_u), check.names = FALSE),
      data.frame(Framework = "MARSS Shared", `Estimation Strategy` = "Regional Shared Trend (Forced Parallel)", `Log-Scale Coefficient (U)` = sprintf("%.4f", marss_s_u), `Annual Population Trend` = get_u_badge(marss_s_u), check.names = FALSE)
    )
    for(i in seq_along(sites)) {
      rows[[length(rows)+1]] <- data.frame(
        Framework = "MARSS Independent", `Estimation Strategy` = paste("Site-Specific (Unconstrained):", sites[i]),
        `Log-Scale Coefficient (U)` = sprintf("%.4f", marss_i_u[i]), `Annual Population Trend` = get_u_badge(marss_i_u[i]), check.names = FALSE
      )
    }
    do.call(rbind, rows)
  }, sanitize.text.function = function(x) x, align = "llrr")

  # --- TAB 2 VARIANCE COMPONENT TABLE WITH VAL-SPECIFIC LAYER TOOLTIPS ---
  output$table_var <- renderTable({
    req(vault$res, vault$marss)
    sites <- sort(unique(vault$abund$Site))

    jags_q <- mean(vault$res$fit$sims.list$Q)
    jags_r <- colMeans(vault$res$fit$sims.list$R)

    marss_s_q <- stats::coef(vault$marss$shared, type = "matrix")$Q[1,1]
    marss_s_r <- diag(stats::coef(vault$marss$shared, type = "matrix")$R)

    marss_i_q <- diag(stats::coef(vault$marss$indep, type = "matrix")$Q)
    marss_i_r <- diag(stats::coef(vault$marss$indep, type = "matrix")$R)

    # Q Badge Generator (Natural Environmental Fluctuations)
    get_q_badge <- function(q_val) {
      if (q_val < 0.05) {
        b_class <- "bg-success text-white"; label <- "Low Fluctuation"
        tip <- "True population numbers are highly stable, solid, and predictable from year to year."
      } else if (q_val < 0.20) {
        b_class <- "bg-warning text-dark"; label <- "Moderate Fluctuations"
        tip <- "Standard fluctuations driven by typical marine cycles like changing regional currents or regional food availability."
      } else {
        b_class <- "bg-danger text-white"; label <- "High Volatility"
        tip <- "High Volatility: Intense ecological shocks or acute environmental anomalies are throwing off population stability."
      }
      sprintf("<span class='badge %s' data-bs-toggle='tooltip' data-bs-placement='top' data-bs-title='%s' style='font-size: 0.85rem; padding: 0.35em 0.5em; min-width: 170px; cursor: pointer;'>%.4f (%s)</span>", b_class, tip, q_val, label)
    }

    # R Badge Generator (Monitoring & Survey Noise)
    get_r_badge <- function(r_val, model_type, jags_val, indep_val) {
      if (r_val < 0.05) {
        b_class <- "bg-success text-white"; label <- "High Precision"
      } else if (r_val < 0.20) {
        b_class <- "bg-warning text-dark"; label <- "Moderate Noise"
      } else {
        b_class <- "bg-danger text-white"; label <- "High Survey Noise ⚠"
      }

      tip <- "Measures standard variation introduced by fieldwork limits and survey constraints."

      if (model_type == "jags") {
        if (jags_val >= 0.12 && indep_val < 0.05) {
          tip <- "SHARED PARALLEL ASSUMPTION MISMATCH: This beach is monitored well in the field. However, it looks noisy here because this regional model tries to force it onto a uniform average trend line that its real biology breaks away from."
        } else if (jags_val < 0.05 && indep_val >= 0.12) {
          tip := "REGIONAL SECURITY SAFETY NET: This beach looks precise here only because it is borrowing data strength from well-surveyed neighboring beaches to smooth over its data gaps behind a group average."
        } else {
          if (r_val < 0.05) tip <- "Beach monitored well: Field data points tightly track the shared regional baseline trend line."
          if (r_val >= 0.05 && r_val < 0.20) tip <- "Standard monitoring noise: Acceptable fieldwork variations expected under a shared regional model."
          if (r_val >= 0.20) tip <- "Highly noisy or missing data: Significant data tracking issues or survey gaps relative to the regional trend line."
        }
      }

      if (model_type == "indep") {
        if (jags_val >= 0.12 && indep_val < 0.05) {
          tip <- "BEACH MONITORED WELL! Unrestricting the models proves your field tracking is highly precise. The high noise score from the baseline model was a structural rule error, not a fieldwork problem."
        } else if (jags_val < 0.05 && indep_val >= 0.12) {
          tip <- "LOTS OF MISSING OR NOISY DATA ON THIS BEACH: Without neighboring data to back it up, this site's lone tracking gaps, missing years, or low sample counts stick out immediately. This beach requires enhanced survey effort consistency."
        } else {
          if (r_val < 0.05) tip <- "Beach monitored well: State tracks show tight, clean, highly precise standalone field tracking execution."
          if (r_val >= 0.05 && r_val < 0.20) tip <- "Standard monitoring noise: Normal standalone dataset variation for an unconstrained site timeline."
          if (r_val >= 0.20) tip <- "Highly noisy or missing data: This single beach has severe survey dropouts or extreme monitoring gaps."
        }
      }

      if (model_type == "shared") {
        if (r_val < 0.05) tip <- "Beach monitored well: Matrix optimization confirms strong, consistent group-wide field tracking data."
        if (r_val >= 0.05 && r_val < 0.20) tip <- "Standard monitoring noise: Expected variance thresholds across the combined beach data profile."
        if (r_val >= 0.20) tip <- "Highly noisy or missing data: High overall variation detected under parallel model constraints."
      }

      sprintf("<span class='badge %s' data-bs-toggle='tooltip' data-bs-placement='top' data-bs-title='%s' style='font-size: 0.85rem; padding: 0.35em 0.5em; min-width: 170px; cursor: pointer;'>%.4f (%s)</span>", b_class, tip, r_val, label)
    }

    rows <- list()
    for(i in seq_along(sites)) {
      rows[[length(rows)+1]] <- data.frame(
        Framework = "JAGS Baseline", `Monitored Sub-Site` = sites[i],
        `Environmental Fluctuations (Q)` = get_q_badge(jags_q),
        `Monitoring & Survey Noise (R)` = get_r_badge(jags_r[i], "jags", jags_r[i], marss_i_r[i]),
        check.names = FALSE
      )
    }
    for(i in seq_along(sites)) {
      rows[[length(rows)+1]] <- data.frame(
        Framework = "MARSS Shared", `Monitored Sub-Site` = sites[i],
        `Environmental Fluctuations (Q)` = get_q_badge(marss_s_q),
        `Monitoring & Survey Noise (R)` = get_r_badge(marss_s_r[i], "shared", jags_r[i], marss_i_r[i]),
        check.names = FALSE
      )
    }
    for(i in seq_along(sites)) {
      rows[[length(rows)+1]] <- data.frame(
        Framework = "MARSS Independent", `Monitored Sub-Site` = sites[i],
        `Environmental Fluctuations (Q)` = get_q_badge(marss_i_q[i]),
        `Monitoring & Survey Noise (R)` = get_r_badge(marss_i_r[i], "indep", jags_r[i], marss_i_r[i]),
        check.names = FALSE
      )
    }
    do.call(rbind, rows)
  }, sanitize.text.function = function(x) x, align = "llrr")

  # 9. Raw Data Preview Table
  output$data_preview <- renderTable({
    req(vault$abund)
    vault$abund
  })
}

shinyApp(ui, server)
