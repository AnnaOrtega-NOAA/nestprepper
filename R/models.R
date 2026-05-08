#' Run Bayesian Trend Analysis
#' @param df Dataframe from calculate_abundance
#' @param iter Number of MCMC iterations
#' @param parallel Logical. Run on multiple cores?
#' @export
run_turtle_model <- function(df, iter = 1000, parallel = FALSE) {

  if (!requireNamespace("jagsUI", quietly = TRUE)) {
    stop("Package 'jagsUI' is needed. Please install it.")
  }

  cat("\n--- Initializing Bayesian Trend Model ---\n")

  # Prepare Data
  wide_data <- df %>%
    dplyr::select(Year, Site, Annual_Nesters) %>%
    tidyr::pivot_wider(names_from = Site, values_from = Annual_Nesters)

  log_mat <- log(as.matrix(wide_data[,-1]))
  log_mat[is.infinite(log_mat) | is.nan(log_mat)] <- NA

  jags_data <- list(
    y = t(log_mat),
    n_yrs = nrow(wide_data),
    n_sites = ncol(log_mat),
    u_mean = 0, u_sd = 0.5,
    q_alpha = 0.01, q_beta = 0.01
  )

  model_string <- "
  model {
    U ~ dnorm(u_mean, 1/(u_sd^2))
    Q_inv ~ dgamma(q_alpha, q_beta)
    Q <- 1/Q_inv
    X[1] ~ dnorm(0, 0.01)
    for(t in 2:n_yrs) {
      X[t] ~ dnorm(X[t-1] + U, Q_inv)
    }
    for(i in 1:n_sites) {
      A[i] ~ dnorm(0, 0.01)
      for(t in 1:n_yrs) {
        y[i,t] ~ dnorm(X[t] + A[i], 100)
      }
    }
  }"

  mod_file <- tempfile(fileext = ".txt")
  writeLines(model_string, mod_file)

  # Explicitly using parallel = FALSE for the first test
  fit <- jagsUI::jags(
    data = jags_data,
    parameters.to.save = c("U", "Q", "X", "A"),
    model.file = mod_file,
    n.chains = 3,
    n.iter = iter,
    n.burnin = floor(iter/2),
    parallel = parallel
  )

  unlink(mod_file)
  return(list(fit = fit, years = wide_data$Year))
}
