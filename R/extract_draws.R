#' Extract Posteriors and Calculate Final Abundance (Technical Memo Standard)
#' @export
get_posteriors <- function(res, remig_int = 3.0) {

  fit <- res$fit
  years <- res$years
  fy <- length(years)

  # 1. Handle site logic safely
  A_sims <- fit$sims.list$A
  if (is.null(dim(A_sims))) {
    A_sims <- matrix(A_sims, ncol = 1)
  }
  n_sites <- ncol(A_sims)

  # 2. Extract final year state (X)
  X_final <- fit$sims.list$X[, fy]

  # 3. Calculate N_fym0 (Total Annual Nesters in Final Year)
  # Replicates: sum(exp(X + A))
  n_draws <- length(X_final)
  X_fym0 <- rep(0, n_draws)
  for (i in 1:n_sites) {
    X_fym0 <- X_fym0 + exp(X_final + A_sims[, i])
  }

  # 4. Extract parameters for the Tech Memo standard
  # Instantaneous growth rate (r) [cite: 336, 1499]
  median_r <- median(fit$sims.list$U)

  # Annual percentage change interpretation [cite: 356, 804, 1015]
  # Formula: (exp(r) - 1) * 100
  pct_change <- round((exp(median_r) - 1) * 100, 2)

  # Format the display string: e.g., "-0.200 (-18.13%)"
  u_display <- paste0(round(median_r, 3), " (", pct_change, "%)")

  # 5. Median Abundance
  median_N_fym0 <- median(X_fym0)
  total_females <- round(median_N_fym0 * remig_int) # [cite: 341]

  # 6. Return standard dataframe for future projections [cite: 731, 2459]
  draws_df <- data.frame(
    U = fit$sims.list$U,
    Q = fit$sims.list$Q,
    N_fym0 = X_fym0
  )

  return(list(
    draws = draws_df,
    summary = list(
      Year = max(years),
      U_display = u_display, # The fixed variable name
      U_pct = pct_change,    # Keep raw % for color-coding logic in UI
      Q_med = round(median(fit$sims.list$Q), 4),
      N_fym0 = round(median_N_fym0),
      Total_Females = total_females
    )
  ))
}
