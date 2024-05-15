compare_plot <- function(object) {
  data_sub <- object$sub
  models <- object$models
  
  max <- max(data_sub$y2, data_sub$y1, data_sub$y1_corr, na.rm = TRUE)
  min <- min(data_sub$y2, data_sub$y1, data_sub$y1_corr, na.rm = TRUE)
  range <- max - min
  
  # Order data for plot
  data_sub <- data_sub[order(data_sub$y2_hat), ]
  
  par(mar = c(3.5, 3.5, 2, 2) + 0.1)
  plot(data_sub$y2_hat, data_sub$y2, pch = 1, cex = 0.5, col = "grey",
       axes = FALSE, xlab = "", ylab = "", ylim = c(min - range * 0.1,
                                                    max + range * 0.2))
  title(main = "Comparison of the methods", cex.main = 0.9)
  ### Add the y axis
  axis(2, col = "black", las = 1)
  mtext("Measurement method", side = 2, line = 2)
  box(col = "black")
  ### Add the x axis
  axis(1)
  mtext("BLUP of x", side = 1, col = "black", line = 2)
  points(data_sub$y2_hat, data_sub$y1, pch = 19, col = "blue", cex = 0.5)
  points(data_sub$y2_hat, data_sub$y1_corr, pch = 19, col = "red",
         cex = 0.5)
  abline(c(0, 1), lwd = 2)
  abline(models[[4]]$coefficients, lwd = 2, lty = 1, col = "blue")
  abline(models[[6]]$coefficients, lwd = 2, lty = 2, col = "red")
  legend("topleft", legend = c("Reference method (y2)", "New method (y1)",
                               "New method(corrected)"),
         pch = c(1, 19, 19), lty = c(1, 1, 2), col = c("black", "blue", "red"),
         y.intersp = 0.7, yjust = 0.2, bty = "n", cex = 0.8)
}