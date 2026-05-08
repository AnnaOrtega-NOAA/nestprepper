#' Plot Population Trends and Status
#' @param model_output Result from run_turtle_model
#' @param abundance_df Result from calculate_abundance
#' @param species_name String for plot title
#' @export
plot_turtle_status <- function(model_output, abundance_df, species_name = "Sea Turtle") {
  
  # 1. Extract the Median and Credible Intervals from JAGS
  # X is the hidden state (log abundance)
  fit <- model_output$fit
  years <- model_output$years
  
  # Convert log-state back to natural scale for plotting
  # We use the median of the posterior distribution
  est_log <- fit$summary[grep("X", rownames(fit$summary)), ]
  est_nat <- exp(est_log[, c("2.5%", "50%", "97.5%")])
  
  plot_data <- data.frame(
    Year = years,
    Lower = est_nat[, 1],
    Median = est_nat[, 2],
    Upper = est_nat[, 3]
  )
  
  # 2. Create the "Ghost Plot" using ggplot2
  library(ggplot2)
  
  p <- ggplot() +
    # The Shaded Uncertainty (The Ghost)
    geom_ribbon(data = plot_data, aes(x = Year, ymin = Lower, ymax = Upper), 
                fill = "grey80", alpha = 0.5) +
    # The Median Trend Line
    geom_line(data = plot_data, aes(x = Year, y = Median), 
              color = "darkblue", size = 1.2) +
    # The Raw Observed Data Points (Annual Nesters)
    geom_point(data = abundance_df, aes(x = Year, y = Annual_Nesters, color = Site), 
               alpha = 0.6, size = 2) +
    labs(title = paste(species_name, "DPS Status & Trend"),
         subtitle = paste0("Annual Trend (U): ", round(fit$mean$U, 3), 
                           " [95% CI: ", round(fit$q2.5$U, 3), " to ", round(fit$q97.5$U, 3), "]"),
         y = "Annual Nesters (Estimated)",
         x = "Season") +
    theme_minimal() +
    scale_color_viridis_d()
  
  print(p)
  
  # 3. Decision Support Question
  cat("\n--- Status Assessment ---")
  if(fit$mean$U < 0) {
    message("\nTrend is NEGATIVE. Would you like to calculate the probability of a 50% decline?")
  } else {
    message("\nTrend is POSITIVE or stable.")
  }
  
  return(p)
}