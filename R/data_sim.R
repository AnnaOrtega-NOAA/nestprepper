#' Simulate Sea Turtle Nesting Data
#' @param n_years Number of years to simulate
#' @param n_sites Number of nesting beaches
#' @param true_trend The annual growth rate (e.g., -0.03 is 3% decline)
#' @export
simulate_turtle_data <- function(n_years = 25, 
                                 n_sites = 3, 
                                 true_trend = -0.03, 
                                 process_error = 0.1, 
                                 nest_to_female = 5.5) {
  
  years <- 2001:(2000 + n_years)
  log_abundance <- numeric(n_years)
  log_abundance[1] <- log(500) 
  
  for(t in 2:n_years) {
    log_abundance[t] <- log_abundance[t-1] + true_trend + rnorm(1, 0, process_error)
  }
  
  data_list <- list()
  for(s in 1:n_sites) {
    site_log_obs <- log_abundance + rnorm(n_years, 0, 0.2)
    site_df <- data.frame(
      Year = years,
      Site = paste("Beach", LETTERS[s]),
      Count = round(exp(site_log_obs) * nest_to_female)
    )
    # Simulate random missing monitoring years
    site_df$Count[sample(1:n_years, size = floor(n_years * 0.1))] <- NA
    data_list[[s]] <- site_df
  }
  
  return(do.call(rbind, data_list))
}