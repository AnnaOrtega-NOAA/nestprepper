#' Plot Turtle Population Status (Visualizing Imputed Trend)
#' @export
plot_turtle_status <- function(res, abund, title = "Population Trend") {
  fit <- res$fit
  years <- res$years
  nsims <- length(fit$sims.list$U)

  # Calculate Total Population from Posterior (Imputation happens here)
  # We sum across all sites using Site Offsets (A) and State (X)
  n_sites <- if(is.matrix(fit$sims.list$A)) ncol(fit$sims.list$A) else 1
  total_abundance_samples <- matrix(0, nrow = nsims, ncol = length(years))

  for(t in 1:length(years)) {
    yearly_sum <- rep(0, nsims)
    for(i in 1:n_sites) {
      current_A <- if(n_sites == 1) fit$sims.list$A else fit$sims.list$A[, i]
      yearly_sum <- yearly_sum + exp(fit$sims.list$X[, t] + current_A)
    }
    total_abundance_samples[, t] <- yearly_sum
  }

  est_nat <- apply(total_abundance_samples, 2, quantile, probs = c(0.025, 0.5, 0.975))

  plot_df <- data.frame(
    Year = years,
    Estimate = est_nat["50%", ],
    Low = est_nat["2.5%", ],
    High = est_nat["97.5%", ]
  )

  # Final Plot: Continuous Ribbon + Line (The Imputation) vs Points (Observations)
  library(ggplot2)
  ggplot() +
    geom_ribbon(data = plot_df, aes(x = Year, ymin = Low, ymax = High), fill = "grey80", alpha = 0.5) +
    geom_line(data = plot_df, aes(x = Year, y = Estimate), color = "darkblue", linewidth = 1.2) +
    geom_point(data = abund, aes(x = Year, y = Annual_Nesters, color = Site), alpha = 0.6) +
    labs(title = title, subtitle = "The blue line and gray band include imputed values for missing data years.",
         y = "Annual Nesters", x = "Year") +
    theme_minimal()
}
