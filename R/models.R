#' Run Bayesian Trend Analysis (Exact singleUQ Baseline)
#' @import dplyr tidyr jagsUI
#' @export
run_turtle_model <- function(df, iter = 150000, burnin = 50000, thin = 10) {

  all_years <- min(df$Year, na.rm = TRUE):max(df$Year, na.rm = TRUE)
  sites <- sort(unique(df$Site))
  n.yrs <- length(all_years)
  n.timeseries <- length(sites)

  # 1. Build Y Matrix safely
  prep_df <- df %>%
    dplyr::group_by(Year, Site) %>%
    dplyr::summarise(Annual_Nesters = if(all(is.na(Annual_Nesters))) NA_real_ else sum(Annual_Nesters, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(Year = all_years, Site = sites) %>%
    dplyr::arrange(Year) %>%
    tidyr::pivot_wider(names_from = Site, values_from = Annual_Nesters) %>%
    dplyr::select(all_of(sites))

  mat_data <- as.matrix(prep_df)
  mat_data[mat_data == 0 | is.nan(mat_data)] <- NA
  Y_matrix <- t(log(mat_data))

  x0_mean <- mean(Y_matrix, na.rm = TRUE)
  if(is.na(x0_mean)) x0_mean <- 5

  # 2. HARDCODED PIFSC-95 PRIORS
  q_alpha <- 0.01; q_beta <- 0.01
  r_alpha <- 0.01; r_beta <- 0.01

  jags_data <- list(
    Y = Y_matrix,
    n.yrs = n.yrs,
    n.timeseries = n.timeseries,
    a_mean = 0, a_tau = 1 / (4^2),
    u_mean = 0, u_tau = 1 / (0.5^2),
    q_alpha = q_alpha, q_beta = q_beta,
    r_alpha = r_alpha, r_beta = r_beta,
    x0_mean = x0_mean, x0_tau = 1 / (10^2)
  )

  # 3. Exact Boyd 2016 State-Space Model
  model_string <- "
  model {
    A[1] <- 0
    for(i in 2:n.timeseries) {
      A[i] ~ dnorm(a_mean, a_tau)
    }
    for(i in 1:n.timeseries) {
      tauR[i] ~ dgamma(r_alpha, r_beta)
      R[i] <- 1/tauR[i]
      for(t in 1:n.yrs) {
        Y[i,t] ~ dnorm(X[t] + A[i], tauR[i])
      }
    }
    tauQ ~ dgamma(q_alpha, q_beta)
    Q <- 1/tauQ
    U ~ dnorm(u_mean, u_tau)
    X0 ~ dnorm(x0_mean, x0_tau)

    X[1] ~ dnorm(X0 + U, tauQ)
    for(t in 2:n.yrs) {
      X[t] ~ dnorm(X[t-1] + U, tauQ)
    }
  }"

  fit <- jagsUI::jags(
    data = jags_data,
    parameters.to.save = c("U", "Q", "R", "X0", "X", "A"),
    model.file = textConnection(model_string),
    n.chains = 2, n.iter = iter, n.burnin = burnin, n.thin = thin, parallel = FALSE
  )

  return(list(fit = fit, years = all_years))
}
