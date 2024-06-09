#' Plot total bias
#' 
#' Plot total bias along with 95% confidence bands
#'
#' @param object list returned by \link{measure_compare} function.
#' @param object (optional) returned by \link{measure_compare} function. 
#' If provided, will plot a second total bias estimate.
#' 
#' @export
#'
#' @examples
#' ### Load the data
#' data(data1)
#' ### Analysis
#' measure_model <- measure_compare(data1)
#' ### Plot the total bias
#' total_bias_plot(measure_model)
total_bias_plot <- function(object, object2 = NULL) {
  print("Generating Total Bias Plot ...")
  
  two_objects <- !is.null(object2)
  
  # Simulate for the first dataset
  object_sim <- total_bias_simulation(object)
  min_bias <- object_sim$min_bias
  max_bias <- object_sim$max_bias
  
  # Simulate for the second dataset (if there is one)
  if (two_objects) {
    object2_sim <- total_bias_simulation(object2)
    min_bias <- min(min_bias, object2_sim$min_bias)
    max_bias <- max(max_bias, object2_sim$max_bias)
  }
  
  graphics::par(mar = c(3.5, 3.5, 3, 4) + 0.1)
  # Plot the bias
  plot(object_sim$data_agg$y2_hat, object_sim$data_agg$bias, xlab = "", 
       ylab = "", axes = FALSE, col = "red", type = "l", lwd = 2, 
       ylim = c(min_bias, max_bias))
  graphics::title(main = "Total bias plot", cex.main = 0.9)
  
  if (two_objects) {
    points(object2_sim$data_agg$y2_hat, object2_sim$data_agg$bias, xlab = "", 
           ylab = "", col = "blue", type = "l", lwd = 2)
  }
  
  if (!two_objects) {
    # Add the subtitle
    subtitle <- paste("Differential bias: ",
                      round(object_sim$bias[1, 1], 3), "; ",
                      "Proportional bias: ", round(object_sim$bias[2, 1], 3),
                      sep = "")
    graphics::mtext(subtitle, side = 3, cex = 0.8)
  }
  
  # Confidence bands
  graphics::points(object_sim$data_agg$y2_hat,
                   object_sim$data_agg$bias_y1_lo_fit, col = "red",
                   type = "l", lty = 2)
  graphics::points(object_sim$data_agg$y2_hat,
                   object_sim$data_agg$bias_y1_up_fit, col = "red",
                   type = "l", lty = 2)
  
  if (two_objects) {
    graphics::points(object2_sim$data_agg$y2_hat,
                     object2_sim$data_agg$bias_y1_lo_fit, col = "blue",
                     type = "l", lty = 2)
    graphics::points(object2_sim$data_agg$y2_hat,
                     object2_sim$data_agg$bias_y1_up_fit, col = "blue",
                     type = "l", lty = 2)
  }
  
  # Zero horizontal line
  graphics::abline(h = 0, col = "wheat2")
  
  # y-graphics::axis
  graphics::axis(2, col = "black", las = 1)
  graphics::mtext("Bias", side = 2, line = 2)
  graphics::box(col = "black")
  
  # x-graphics::axis
  graphics::axis(1)
  graphics::mtext("BLUP of x", side = 1, col = "black", line = 2)
  
  # Legend
  if (!two_objects) {
    graphics::legend("top", legend = c("Bias model 1", "95%CB"),
                     pch = c(1, 19), col = c("red", "red"),
                     pt.cex = c(0, 0), y.intersp = 0.7, yjust = 0.2,
                     lty = c(1, 2), bty = "n", cex = 0.8)
  } else {
    graphics::legend("top", legend = c("Bias model 1", "95%CB", "Bias model 2",
                                       "95%CB"),
                     pch = c(1, 19), col = c("red", "red", "blue", "blue"),
                     pt.cex = c(0, 0), y.intersp = 0.7, yjust = 0.2,
                     lty = c(1, 2, 1, 2), bty = "n", cex = 0.8)
  }
}
