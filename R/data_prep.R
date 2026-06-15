#' Aggregate Monthly Data using Exact Tomo Eguchi AR1 Bayesian Fourier Imputation
#' @import dplyr tidyr jagsUI
#' @export
aggregate_monthly_to_annual <- function(df, quiet = FALSE) {

  print(">>> SUCCESS: EXACT TOMO JAGS ENGINE INITIATED <<<")

  # 1. Sanitize Input
  df <- df %>%
    dplyr::mutate(Year = as.numeric(Year), Count = as.numeric(Count)) %>%
    dplyr::filter(!is.na(Year), !is.na(Site))

  # If no months are provided, skip imputation
  if (!"Month" %in% names(df)) {
    return(df %>%
             dplyr::group_by(Year, Site) %>%
             dplyr::summarise(Count = if(all(is.na(Count))) NA_real_ else sum(Count, na.rm = TRUE), .groups = "drop"))
  }

  if (!quiet) message("--- Running Tomo's Exact AR1 Bayesian Fourier Imputation ---")

  df <- df %>% dplyr::mutate(Month = as.numeric(Month)) %>% dplyr::filter(!is.na(Month))

  # CRITICAL GUARD: Track years that have ZERO real data
  valid_years <- df %>%
    dplyr::group_by(Year, Site) %>%
    dplyr::summarise(valid_months = sum(!is.na(Count) & Count >= 0), .groups = "drop")

  df_clean <- df %>%
    dplyr::group_by(Year, Month, Site) %>%
    dplyr::summarise(Count = if(all(is.na(Count))) NA_real_ else sum(Count, na.rm = TRUE), .groups = "drop")

  all_years <- min(df_clean$Year, na.rm = TRUE):max(df_clean$Year, na.rm = TRUE)
  sites <- sort(unique(df_clean$Site))
  n_years <- length(all_years)
  n_timeseries <- length(sites)

  # 2. Recreate Calendar Grid
  full_grid <- expand.grid(Month = 1:12, Year = all_years, Site = sites, stringsAsFactors = FALSE)

  prep_df <- full_grid %>%
    dplyr::left_join(df_clean, by = c("Year", "Month", "Site")) %>%
    dplyr::arrange(Year, Month) %>%
    tidyr::pivot_wider(names_from = Site, values_from = Count) %>%
    dplyr::select(all_of(sites))

  y_matrix <- as.matrix(prep_df)
  y_matrix[y_matrix == 0] <- NA
  y_matrix <- log(y_matrix)

  # 3. Period Logic
  periods <- rep(12, n_timeseries)
  for(i in 1:n_timeseries) {
    # Added your three new beaches to trigger the 6-month engine automatically
    if(grepl("W_|Wermon|Wamlana|Waspait|Waenibe", sites[i], ignore.case = TRUE)) {
      periods[i] <- 6
    }
  }

  jags_data <- list(
    y = y_matrix,
    m = rep(1:12, times = n_years),
    n.steps = nrow(y_matrix),
    n.months = 12,
    pi = pi,
    period = periods,
    n.timeseries = n_timeseries,
    n.years = n_years
  )

  # 4. EXACT TOMO EGUCHI MODEL STRING
  model_string <- "
  model{
    for(j in 1:n.timeseries) {
       predX0[j] ~ dnorm(5, 0.1)
       predX[1,j] <- c[j, m[1]] + predX0[j]
       X[1,j] ~ dnorm(predX[1,j], tau.X[j])
       y[1,j] ~  dnorm(X[1,j], tau.y[j])

       for (t in 2:n.steps){
           predX[t,j] <-  c[j,m[t]] + X[t-1, j]
           X[t,j] ~ dnorm(predX[t,j], tau.X[j])
           y[t,j] ~  dnorm(X[t,j], tau.y[j])
        }

        for (y in 1:n.years){
           for (mm in 1:12){
              tmp2[y, mm, j] <- exp(X[(y*12 - mm + 1), j])
           }
           N[y, j] <- log(sum(tmp2[y,,j]))
        }
    }

    for (j in 1:n.timeseries){
        for (k in 1:n.months){
            c.const[j, k] <-  2 * pi * k / period[j]
            c[j, k] <- beta.cos[j] * cos(c.const[j,k]) + beta.sin[j] * sin(c.const[j,k])
        }

        sigma.y[j] ~ dgamma(2, 0.5)
        tau.y[j] <- 1/(sigma.y[j] * sigma.y[j])
        beta.cos[j] ~ dnorm(0, 1)
        beta.sin[j] ~ dnorm(0, 1)

        sigma.X[j] ~ dgamma(2, 0.5)
        tau.X[j] <- 1/(sigma.X[j] * sigma.X[j])
    }
  }"

  # 5. Execute - NOW SAVING MONTHLY IMPUTATIONS (tmp2)
  jm <- jagsUI::jags(
    data = jags_data,
    parameters.to.save = c("N", "tmp2"),
    model.file = textConnection(model_string),
    n.chains = 5, n.iter = 100000, n.burnin = 50000, n.thin = 5, parallel = FALSE
  )

  # Extract the 3D array of imputed monthly counts [Year, Month, Site]
  tmp2_sims <- jm$q50$tmp2

  # Format into a clean dataframe
  dimnames(tmp2_sims) <- list(all_years, 1:12, sites)
  res_monthly <- as.data.frame(as.table(tmp2_sims))
  colnames(res_monthly) <- c("Year", "Month", "Site", "Count")
  res_monthly$Year <- as.numeric(as.character(res_monthly$Year))
  res_monthly$Month <- as.numeric(as.character(res_monthly$Month))
  res_monthly$Site <- as.character(res_monthly$Site)

  # 6. CRITICAL QA/QC: Trim Pre-Monitoring Hallucinations
  first_valid <- df %>%
    dplyr::filter(!is.na(Count) & Count > 0) %>%
    dplyr::group_by(Site) %>%
    dplyr::summarise(first_year = min(Year, na.rm=TRUE), .groups = "drop")

  res_monthly <- res_monthly %>%
    dplyr::left_join(first_valid, by = "Site") %>%
    dplyr::mutate(
      Count = ifelse(Year < first_year, NA_real_, Count)
    ) %>%
    dplyr::select(-first_year)

  if (!quiet) message("--- Imputation Complete ---")
  return(res_monthly) # Now returning mathematically complete monthly data
}
