#' Run Comparative MARSS Models (Shared and Independent)
#' @param df Cleaned abundance data frame from calculate_abundance()
#' @return A list containing both the shared-trend and independent-trend MARSS fit objects
#' @import MARSS dplyr tidyr
#' @export
run_marss_models <- function(df) {
  
  # 1. Reshape the data into the exact wide matrix format required by MARSS
  all_years <- min(df$Year, na.rm = TRUE):max(df$Year, na.rm = TRUE)
  sites <- sort(unique(df$Site))
  n_sites <- length(sites)
  
  prep_df <- df %>%
    dplyr::group_by(Year, Site) %>%
    dplyr::summarise(Annual_Nesters = if(all(is.na(Annual_Nesters))) NA_real_ else sum(Annual_Nesters, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(Year = all_years, Site = sites) %>%
    dplyr::arrange(Year) %>%
    tidyr::pivot_wider(names_from = Site, values_from = Annual_Nesters) %>%
    dplyr::select(all_of(sites))
  
  # MARSS expects time across columns, so we transpose the matrix
  # and work in log-space to handle the exponential growth/decline
  Y_matrix <- t(log(as.matrix(prep_df)))
  Y_matrix[Y_matrix == 0 | is.nan(Y_matrix)] <- NA
  
  #---------------------------------------------------------
  # MODEL 1: MARSS Shared Trend (Mirrors JAGS singleUQ)
  #---------------------------------------------------------
  # Hidden States (X): All sites share 1 underlying true population state trajectory
  # Site Offsets (A): Fixed constants relative to the baseline anchor site ("scaling")
  # Observation Errors (R): Unique to each site
  # Process Variance (Q): A single shared scalar variance parameter
  
  model_shared <- list(
    Z = matrix(1, nrow = n_sites, ncol = 1), # Map all sites to 1 state
    A = "scaling",                          # ALLOWS UNIQUE SITE SIZES RELATIVE TO ANCHOR
    R = "diagonal and unequal",              # Each beach has its own counting error
    Q = matrix("q"),                         # 1 shared process variance parameter
    U = matrix("u")                          # 1 shared regional trend parameter
  )
  
  fit_shared <- MARSS::MARSS(Y_matrix, model = model_shared, silent = TRUE)
  
  #---------------------------------------------------------
  # MODEL 2: MARSS Independent Trends (The Diagnostic Breakdown)
  #---------------------------------------------------------
  # Hidden States (X): Each site gets its own completely unlinked state line
  # Site Offsets (A): Dropped (set to 0) because baseline scaling happens via unique states
  # Observation Errors (R): Unique to each site
  # Process Variance (Q): Completely separate, unlinked variances for each site
  
  model_indep <- list(
    Z = diag(1, n_sites),                    # Each site maps to its own state line
    A = "zero",                              # Absorbed directly into the unique states
    R = "diagonal and unequal",
    Q = "diagonal and unequal",              # Independent environmental fluctuations
    U = "unequal"                            # Separate growth rates (U) calculated for each site
  )
  
  fit_indep <- MARSS::MARSS(Y_matrix, model = model_indep, silent = TRUE)
  
  return(list(
    shared = fit_shared,
    indep = fit_indep,
    years = all_years,
    sites = sites
  ))
}