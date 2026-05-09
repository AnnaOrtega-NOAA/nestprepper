#' Run Bayesian Trend Analysis (Stabilized for Imputation)
#' @export
run_turtle_model <- function(df, iter = 100000, parallel = FALSE) {

  # 1. Ensure the grid is full for imputation
  all_years <- min(df$Year):max(df$Year)
  all_sites <- unique(df$Site)
  full_grid <- expand.grid(Year = all_years, Site = all_sites)

  prepared_df <- full_grid %>%
    dplyr::left_join(df, by = c("Year", "Site")) %>%
    dplyr::group_by(Year, Site) %>%
    dplyr::summarise(Annual_Nesters = sum(Annual_Nesters, na.rm = FALSE), .groups = "drop")

  wide_data <- prepared_df %>%
    tidyr::pivot_wider(names_from = Site, values_from = Annual_Nesters)

  # 2. Matrix Prep with Log-Offset (Prevents Invalid Parent Values)
  Y <- t(log(as.matrix(wide_data[,-1]) + 0.01))
  Y[is.na(Y)] <- NA # Restore NAs where data is actually missing

  n_yrs <- ncol(Y)
  n_timeseries <- nrow(Y)

  jags_data <- list(Y = Y, n.yrs = n_yrs, n.timeseries = n_timeseries)

  # 3. Stabilized Model String (Mirroring 2026 methodology)
  model_string <- "
  model {
    # Slightly broader prior to help convergence during burn-in
    U ~ dnorm(0, 10)

    # Process Variance (Q) - Changed to a more stable Gamma prior
    tauQ ~ dgamma(0.1, 0.1)
    Q <- 1/tauQ

    # State Process
    X[1] ~ dnorm(0, 0.001)
    for(t in 2:n.yrs) {
      X[t] ~ dnorm(X[t-1] + U, tauQ)
    }

    # Observation Process (A offsets)
    for(i in 1:n.timeseries) {
      A[i] ~ dnorm(0, 0.001)
      for(t in 1:n.yrs) {
        # Fixed precision 100 matches your singleUQ.txt standard
        Y[i,t] ~ dnorm(X[t] + A[i], 100)
      }
    }
  }"

  mod_file <- tempfile(fileext = ".txt")
  writeLines(model_string, mod_file)

  fit <- jagsUI::jags(
    data = jags_data,
    parameters.to.save = c("U", "Q", "X", "A"),
    model.file = mod_file,
    n.chains = 3, n.iter = iter, n.burnin = floor(iter/2), n.thin = 50, parallel = parallel
  )

  unlink(mod_file)
  return(list(fit = fit, years = wide_data$Year))
}
