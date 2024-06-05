#' Plot the bias and measurements
#' 
#' Plot the bias and measurements from the two methods along with their 
#' regression lines.
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
#' ### Plot the bias
#' bias_plot(measure_model)
bias_plot <- function(object) {
  print("Generating Bias Plot ...")
  
  # Extract the objects from the output
  bias <- object$bias
  data_sub <- object$sub
  data_new <- object$new
  models <- object$models
  subtitle <- paste("Differential bias:", round(bias[1, 1], 3), "; ",
                    "Proportional bias:", round(bias[2, 1], 3), sep = "")
  
  data_new$fit_y1 <- predict(models[[4]], data_new)
  data_new$fit_y2 <- predict(models[[2]], data_new)
  
  # Compute min and max values for y-axis
  min_y <- min(data_sub$y2, data_sub$y1, na.rm = TRUE)
  max_y <- max(data_sub$y2, data_sub$y1, na.rm = TRUE)
  min_y <- floor(min_y)
  range <- max_y - min_y
  max_y <- ceiling(max_y + range * 0.2)
  
  # Order data for plot
  data_sub <- data_sub[order(data_sub$y2_hat), ]
  data_new <- data_new[order(data_new$y2_hat), ]
  
  par(mar = c(3.5, 3.5, 3, 4) + 0.1)
  ### plot scatter for y1 and y2 with respect to BLUP of x
  plot(data_sub$y2_hat, data_sub$y2, axes = FALSE, xlab = "", ylab = "",
       col = "grey", ylim = c(min_y, max_y), cex = 0.5)
  title(main = "Bias", cex.main = 0.9)
  points(data_sub$y2_hat, data_sub$y1, col = "blue", pch = 19, cex = 0.5)
  
  ### Add the fitted line from model 2 and model 4,
  abline(a = 0, b = models[[2]]$coef, lwd = 2, lty = 1, col = "black")
  abline(models[[4]]$coef, lwd = 2, lty = 2, col = "blue")
  
  ### Add the subtitle
  mtext(subtitle, side = 3, cex = 0.8)
  ### Add the left y axis
  axis(2, col = "black", las = 1)
  mtext("y1 and y2", side = 2, line = 2)
  box(col = "black")
  
  ### Add second plot: bias axis
  par(new = TRUE)
  ## Plot the bias plot and put axis scale on right
  plot(data_new$y2_hat, data_new$bias, xlab = "", ylab = "", axes = FALSE,
       col = "red", type = "l", lty = 1, lwd = 2)
  abline(h = 0, col = "black", lwd = 1)
  ## Add the right y axis and label
  mtext("Bias", side = 4, col = "red", line = 2.5)
  axis(4, col = "red", col.axis = "red", las = 1)
  ### Draw the x axis and add the label
  axis(1)
  mtext("BLUP of x", side = 1, col = "black", line = 2)
  ### Add the legend
  lm_bias <- lm(data_new$bias ~ data_new$y2_hat)
  if (coef(lm_bias)[2] < 0) {
    legend("top", legend = c("Reference method (y2)", "New method (y1)",
                             "Bias"),
           pch = c(1, 19), col = c("black", "blue", "red"),
           pt.cex = c(1, 1, 0.0001), y.intersp = 0.7, yjust = 0.2,
           lty = c(1, 2, 1), bty = "n", cex = 0.8)
  } else {
    legend("topleft", legend = c("Reference method (y2)", "New method (y1)",
                                 "Bias"),
           pch = c(1, 19), col = c("black", "blue", "red"),
           pt.cex = c(1, 1, 0.0001), y.intersp = 0.7, yjust = 0.2,
           lty = c(1, 2, 1), bty = "n", cex = 0.8)
  }
}