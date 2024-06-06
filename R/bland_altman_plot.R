#' Bland and Atlman's limits of agreement plot
#'
#' This function produces Bland and Altman's Limits of Agreement plot
#' (LoA) when there are repeated measurements with possibly heteroscedastic
#' measurement errors.
#'
#'
#' @param data a dataframe contains the object identification number (id),
#' the measurement values from new measurement method (y1) and those from
#' the reference standard (y2)
#' @param new specify the variable name or location for the new measurement
#' method
#' @param ref specify the variable name or location for the reference
#' measurement method
#' @param id specify the variable name for location for the subject
#' identification number (id)
#' @param fill logical. if \code{TRUE} use the average value for new methods to
#' fill out the missing value (only useful for drawing a plot with all the
#' measurements by the reference standard)
#'
#' @author  Mingkai Peng
#' @details  This functions computes the limits of agreement (LoA) when there
#' are repeated measurements and possibly the measurement error are
#' heteroscedastic
#' @export
#' @importFrom stats lm
#' @examples
#' ### Load the data
#' data(data1)
#' ### Bland and Altman's plot
#' bland_altman_plot(data1)

bland_altman_plot <- function(data, new = "y1", ref = "y2", id = "id",
                              fill = TRUE) {
  print("Generating Bland and Altman extended LoA Plot")
  
  data_sub <- (data[, c(id, new, ref)])
  colnames(data_sub) <- c("id", "y1", "y2")
  #### calculate the difference and average
  data_sub$diff_m <- data_sub$y1 - data_sub$y2
  data_sub$avg_m <- (data_sub$y1 + data_sub$y2) /  2
  #### Calculate the harmonic means
  y1 <- as.numeric(data.frame(table(stats::na.omit(data_sub)$id))$Freq)
  y2 <- as.numeric(data.frame(table(data_sub$id))$Freq)
  coef_1 <- 1 / mean(1 / y1)
  coef_1 <- 1 - 1 / coef_1
  coef_2 <- 1 / mean(1 / y2)
  coef_2 <- 1 - 1 / coef_2
  #### With-object variance
  model_y2 <- lme4::lmer(y2 ~ 1 + (1 | id), data = data_sub,
                         na.action = stats::na.exclude)
  vif <- sigma(model_y2)^2 * coef_2
  if (max(y1) > 1) {
    model_y1 <- lme4::lmer(y1 ~ 1 + (1 | id), data = data_sub,
                           na.action = stats::na.exclude)
    vif <- sigma(model_y2)^2 * coef_2 + sigma(model_y1)^2 * coef_1
  }
  ####    model on difference based on average
  data_sub_1 <- stats::na.omit(data_sub)
  model_1 <- lm(diff_m ~ avg_m, data = data_sub_1)
  data_sub_1$fitted <- fitted(model_1)
  data_sub_1$resid_fitted <- residuals(model_1)
  data_sub_1$resid_fitted_abs <- abs(data_sub_1$resid_fitted)
  #####   Regression on absolute residuals based on average
  model_2 <- lm(resid_fitted_abs ~ avg_m, data = data_sub_1)
  data_sub_1$sig2_abs_res <- fitted(model_2) * sqrt(pi / 2)
  data_sub_1$sig2_abs_res <- data_sub_1$sig2_abs_res^2
  data_sub_1$upper <- data_sub_1$fitted +
    1.96 * sqrt(data_sub_1$sig2_abs_res + vif)
  data_sub_1$lower <- data_sub_1$fitted -
    1.96 * sqrt(data_sub_1$sig2_abs_res + vif)
  ####y
  if (fill) {
    mean_y1 <- aggregate(y1 ~ id, data = data_sub, mean)
    colnames(mean_y1)[2] <- "y1.mean"
    data_sub <- merge(data_sub, mean_y1, by = "id")
    data_sub$y1 <- ifelse(is.na(data_sub$y1), data_sub$y1.mean, data_sub$y1)
    data_sub$diff_m <- data_sub$y1 - data_sub$y2
    data_sub$avg_m <- (data_sub$y1 + data_sub$y2) / 2
  }
  ### model for plot
  model_3 <- lm(upper ~ avg_m, data = data_sub_1)
  model_4 <- lm(lower ~ avg_m, data = data_sub_1)
  max <- max(abs(data_sub$diff_m), na.rm = TRUE)
  max <- max + max * 0.2
  ##### final plot
  graphics::par(mar = c(3.5, 3.5, 2, 2) + 0.1)
  plot(data_sub$avg_m, data_sub$diff_m, xlab = "", ylab = "", axes = FALSE,
       col = "grey", ylim = c(-max, max))
  graphics::title(main = "Extended Bland & Altman LoA plot",
        cex.main = 0.9)
  ### Add the y graphics::axis
  graphics::axis(2, col = "black", las = 1)
  graphics::mtext("Difference of y1 and y2 (y1-y2)", side = 2, line = 2)
  graphics::box(col = "black")
  ### Add the x graphics::axis
  graphics::axis(1)
  graphics::mtext("Mean of y1 and y2 [(y1+y2)/2]", side = 1, col = "black", line = 2)
  graphics::abline(model_1$coefficients, col = "blue", lwd = 2)
  graphics::abline(h = 0, col = "black", lwd = 2)
  graphics::abline(model_3$coefficients, col = "blue", lty = 2, lwd = 2)
  graphics::abline(model_4$coefficients, col = "blue", lty = 2, lwd = 2)
  graphics::legend("topleft", legend = c("Regression line", "95% LoA"),
         col = c("blue", "blue"), y.intersp = 0.7, yjust = 0.2, lty = c(1, 2),
         bty = "n", cex = 0.8)
}
