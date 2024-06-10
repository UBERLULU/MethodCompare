#' Plot the precision of the methods
#' 
#' This plot allows the visual comparison of the precision (i.e. standard 
#' deviation) of the new measurement method with the reference standard by 
#' creating a scatter plot of the estimated standard deviations, along with 
#' their 95% confidence bands, against the best 
#' linear prediction (BLUP) of the true latent trait x.
#'
#' @inheritParams total_bias_plot
#' 
#' @export
#'
#' @examples
#' ### Load the data
#' data(data1)
#' ### Analysis
#' measure_model <- measure_compare(data1)
#' ### Plot the precision of the two methods
#' precision_plot(measure_model)
precision_plot <- function(object, object2 = NULL) {
  print("Generating Precision Plot ...")
  
  two_objects <- !is.null(object2)
  
  # Simulate for the first dataset
  object_sim <- precision_simulation(object)
  min_y <- object_sim$min_y
  max_y <- object_sim$max_y
  
  # Simulate for the second dataset (if there is one)
  if (two_objects) {
    object2_sim <- precision_simulation(object2)
    min_y <- min(min_y, object2_sim$min_y)
    max_y <- max(max_y, object2_sim$max_y)
  }
  
  graphics::par(mar = c(3.5, 3.5, 3, 4) + 0.1)
  # Plot the precision
  plot(object_sim$data_agg$y2_hat, exp(object_sim$data_agg$log_sig_res_y1_corr),
       xlab = "", ylab = "", axes = FALSE, col = "red", type = "l", lwd = 2,
       ylim = c(min_y, max_y))
  graphics::title(main = "Precision plot", cex.main = 0.9)
  
  if (two_objects) {
    points(object2_sim$data_agg$y2_hat,
           exp(object2_sim$data_agg$log_sig_res_y1_corr), xlab = "", ylab = "",
           col = "blue", type = "l", lwd = 2)
  }
  
  # Add the subtitle
  subtitle <- "with 95% confidence bands (after recalibration)"
  graphics::mtext(subtitle, side = 3, cex = 0.8)
  
  # Confidence bands
  graphics::points(object_sim$data_agg$y2_hat,
                   object_sim$data_agg$sig_e1_corr_lo_fit, col = "red",
                   type = "l", lty = 2)
  graphics::points(object_sim$data_agg$y2_hat,
                   object_sim$data_agg$sig_e1_corr_up_fit, col = "red",
                   type = "l", lty = 2)
  
  if (two_objects) {
    graphics::points(object2_sim$data_agg$y2_hat,
                     object2_sim$data_agg$sig_e1_corr_lo_fit, col = "blue",
                     type = "l", lty = 2)
    graphics::points(object2_sim$data_agg$y2_hat,
                     object2_sim$data_agg$sig_e1_corr_up_fit, col = "blue",
                     type = "l", lty = 2)
  }
  
  # y2
  graphics::points(object_sim$data_agg$y2_hat,
                   exp(object_sim$data_agg$log_sig_res_y2), col = "black",
                   type = "l", lwd = 2)
  
  # Confidence bands
  graphics::points(object_sim$data_agg$y2_hat,
                   object_sim$data_agg$sig_e2_lo_fit, col = "black",
                   type = "l", lty = 2)
  graphics::points(object_sim$data_agg$y2_hat,
                   object_sim$data_agg$sig_e2_up_fit, col = "black",
                   type = "l", lty = 2)
  
  # y-graphics::axis
  graphics::axis(2, col = "black", las = 1)
  graphics::mtext("Standard deviation of the measurementr errors", side = 2, line = 2)
  graphics::box(col = "black")
  
  # x-graphics::axis
  graphics::axis(1)
  graphics::mtext("BLUP of x", side = 1, col = "black", line = 2)
  
  # Legend
  if (!two_objects) {
    graphics::legend("top",
                     legend = c("Reference standard: y2", 
                                "Recalibrated new method: y1_corr"),
                     pch = c(1, 19), col = c("black", "red"), pt.cex = c(0, 0),
                     y.intersp = 0.7, yjust = 0.2, lty = c(1, 1), bty = "n",
                     cex = 0.8)
  } else {
    graphics::legend("top",
                     legend = c("Reference standard: y2", 
                                "Recalibrated new method of model 1: y1_corr_1",
                                "Recalibrated new method of model 2: y1_corr_2"),
                     pch = c(1, 19), col = c("black", "red", "blue"), pt.cex = c(0, 0),
                     y.intersp = 0.7, yjust = 0.2, lty = c(1, 1, 1), bty = "n",
                     cex = 0.8)
  }
}