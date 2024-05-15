total_bias_plot <- function(object) {
  print("Generating Total Bias Plot ...")
  
  # Extract the objects from the output
  bias <- object$bias
  data_agg <- object$agg
  nb_simul <- object$nb_simul
  
  # Retrieve useful params for simulation
  sim_params <- object$sim_params
  
  m1 <- matrix(data = sim_params$model_4_coef, nrow = 1)
  v1 <- matrix(data = sim_params$model_4_cov, ncol = 2)
  
  # Initiate variables
  sim_max_d <- vector(mode = "list", length = nb_simul)
  
  for (j in 1:nb_simul) {
    blup_x_j <- rnorm(dim(data_agg)[1], mean = data_agg$fitted_y2,
                      sd = data_agg$sd_blup)
    biases_j <- rockchalk::mvrnorm(dim(data_agg)[1], mu = m1, Sigma = v1)
    
    v_bias_j <- (biases_j[, 2] - 1)^2 * data_agg$v_blup + v1[1, 1] +
      v1[2, 2] * (data_agg$v_blup + blup_x_j^2) + 2 * v1[1, 2] * blup_x_j
    bias_j <- biases_j[, 1] + blup_x_j * (biases_j[, 2] - 1)
    
    d_j <- abs(bias_j - data_agg$bias_y1) / sqrt(v_bias_j)
    max_d_j <- max(d_j)
    
    sim_max_d[[j]] <- max_d_j
  }
  
  crit_value2 <- quantile(unlist(sim_max_d), c(0.95))
  
  data_agg$bias_y1_lo <- data_agg$bias_y1 -
    crit_value2 * data_agg$se_bias_y1
  data_agg$bias_y1_up <- data_agg$bias_y1 +
    crit_value2 * data_agg$se_bias_y1
  
  fp <- function(...) mfp::fp(...)
  
  frac_poly_bias_y1_lo <- mfp::mfp(bias_y1_lo ~ fp(fitted_y2, df = 4),
                              data = data_agg)
  frac_poly_bias_y1_up <- mfp::mfp(bias_y1_up ~ fp(fitted_y2, df = 4),
                              data = data_agg)
  
  data_agg$bias_y1_lo_fit <- predict(frac_poly_bias_y1_lo)
  data_agg$bias_y1_up_fit <- predict(frac_poly_bias_y1_up)
  
  min_bias <- min(data_agg$bias_y1_lo_fit, data_agg$bias_y1_up_fit,
                  data_agg$bias_y1, na.rm = TRUE)
  max_bias <- max(data_agg$bias_y1_lo_fit, data_agg$bias_y1_up_fit,
                  data_agg$bias_y1, na.rm = TRUE)
  
  min_bias <- floor(min_bias)
  range <- max_bias - min_bias
  max_bias <- ceiling(max_bias + range * 0.2)
  
  # Order data for plot
  data_agg <- data_agg[order(data_agg$y2_hat), ]
  
  par(mar = c(3.5, 3.5, 3, 4) + 0.1)
  # Plot the bias
  plot(data_agg$y2_hat, data_agg$bias, xlab = "", ylab = "", axes = FALSE,
       col = "red", type = "l", lwd = 2, ylim = c(min_bias, max_bias))
  title(main = "Total bias plot", cex.main = 0.9)
  
  # Add the subtitle
  subtitle <- paste("Differential bias: ", round(bias[1, 1], 3), "; ",
                    "Proportional bias: ", round(bias[2, 1], 3), sep = "")
  mtext(subtitle, side = 3, cex = 0.8)
  
  # Confidence bands
  points(data_agg$y2_hat, data_agg$bias_y1_lo_fit, col = "red",
         type = "l", lty = 2)
  points(data_agg$y2_hat, data_agg$bias_y1_up_fit, col = "red",
         type = "l", lty = 2)
  
  # Zero horizontal line
  abline(h = 0, col = "wheat2")
  
  # y-axis
  axis(2, col = "black", las = 1)
  mtext("Bias", side = 2, line = 2)
  box(col = "black")
  
  # x-axis
  axis(1)
  mtext("BLUP of x", side = 1, col = "black", line = 2)
  
  # Legend
  legend("top", legend = c("Bias", "95%CB"), pch = c(1, 19),
         col = c("red", "red"), pt.cex = c(0, 0), y.intersp = 0.7, yjust = 0.2,
         lty = c(1, 2), bty = "n", cex = 0.8)
}