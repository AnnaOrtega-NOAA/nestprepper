#' Plot Turtle Population Status (Exact singleUQ Aesthetics)
#' @export
plot_turtle_status <- function(res, abund) {
  fit <- res$fit
  years <- res$years
  n_yrs <- length(years)

  # 1. Handle site logic
  A_sims <- fit$sims.list$A
  if (is.null(dim(A_sims))) {
    A_sims <- matrix(A_sims, ncol = 1)
  }
  n_sites <- ncol(A_sims)

  # 2. Calculate Totals (Sum exp(X + A))
  X0_sims <- fit$sims.list$X0
  X_sims  <- fit$sims.list$X

  X0_total <- rep(0, length(X0_sims))
  X_total  <- matrix(0, nrow = nrow(X_sims), ncol = ncol(X_sims))

  for (i in 1:n_sites) {
    X0_total <- X0_total + exp(X0_sims + A_sims[, i])
    X_total  <- X_total + exp(sweep(X_sims, 1, A_sims[, i], "+"))
  }

  # 3. Get Log Quantiles
  X_q <- apply(log(X_total), 2, quantile, probs = c(0.025, 0.5, 0.975))

  med_line  <- c(median(log(X0_total)), X_q[2, ])
  low_line  <- c(quantile(log(X0_total), 0.025), X_q[1, ])
  high_line <- c(quantile(log(X0_total), 0.975), X_q[3, ])

  # 4. Process Observed Data safely
  # Force matching of observed years to the model timeline
  obs_summary <- abund %>%
    dplyr::group_by(Year) %>%
    dplyr::summarise(Annual_Total = sum(Annual_Nesters, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Annual_Total = ifelse(Annual_Total == 0, NA, Annual_Total))

  log_obs <- log(obs_summary$Annual_Total)
  obs_years <- obs_summary$Year

  # 5. Set plot limits
  yrange_vals <- c(log_obs, low_line, high_line)
  yrange <- range(yrange_vals[is.finite(yrange_vals)], na.rm = TRUE)
  plot_years <- c(min(years) - 1, years)

  # 6. Density Data
  den_X0 <- density(log(X0_total), adj = 2)
  den_X0$y2 <- den_X0$y / max(den_X0$y)
  den_N0 <- density(log(X_total[, n_yrs]), adj = 2)
  den_N0$y2 <- den_N0$y / max(den_N0$y)

  # ---------------------------------------------------------
  # 7. BUILD THE BASE R PLOT
  # ---------------------------------------------------------
  par(mar = c(4, 4, 1, 1))

  plot(plot_years, med_line, type = "n", las = 1,
       ylab = "log(Annual Nesters)", xlab = "Season",
       ylim = yrange,
       xlim = c(min(plot_years) - 0.1, max(plot_years) + 0.5),
       xaxs = 'i')

  polygon(c(plot_years, rev(plot_years)),
          c(low_line, rev(high_line)),
          col = 'grey85', border = FALSE)

  lines(plot_years, med_line, lwd = 3, col = 'gray50')

  # X0 Ghost
  q_X0 <- quantile(log(X0_total), probs = c(0.025, 0.975))
  xid <- sapply(q_X0, function(x) which.min(abs(x - den_X0$x)))
  polygon(x = c(rep(plot_years[1], length(xid[1]:xid[2])),
                rev(den_X0$y2[xid[1]:xid[2]] + plot_years[1])),
          y = c(den_X0$x[xid[1]:xid[2]], rev(den_X0$x[xid[1]:xid[2]])),
          border = FALSE, col = rgb(30/255, 144/255, 255/255, 0.3))

  # N_final Ghost
  q_N0 <- quantile(log(X_total[, n_yrs]), probs = c(0.025, 0.975))
  xid_n0 <- sapply(q_N0, function(x) which.min(abs(x - den_N0$x)))
  polygon(x = c(rep(max(plot_years), length(xid_n0[1]:xid_n0[2])),
                rev(-den_N0$y2[xid_n0[1]:xid_n0[2]] + max(plot_years))),
          y = c(den_N0$x[xid_n0[1]:xid_n0[2]], rev(den_N0$x[xid_n0[1]:xid_n0[2]])),
          border = FALSE, col = rgb(153/255, 50/255, 204/255, 0.3))

  # FIX: Observed Data Points plotted specifically against obs_years
  points(obs_years, log_obs, pch = 16, col = "black")

  # Median Points
  points(plot_years, med_line, pch = 16,
         col = c('dodgerblue3', rep('red', length(years) - 1), "darkorchid3"))

  legend("topright",
         legend = c("Observed Data", "Median r", "95% CI Ribbon", "T0 Density", "N_final Density"),
         pch = c(16, NA, 15, 15, 15), lwd = c(NA, 3, NA, NA, NA),
         pt.cex = c(1, NA, 2, 2, 2),
         col = c("black", "gray50", "grey85", "dodgerblue3", "darkorchid3"),
         bty = "n")
}
