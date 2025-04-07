#' Simple Slope Plot Function
#'
#' This function generates a simple slope plot using ggplot2 for a moderated regression model.
#'
#' @param b0 Intercept of the outcome variable. Note: the default intercept is 0, which corresponds to the default identification in Mplus for latent outcomes. Adjust b0 according to your own model results. In manifest models, the intercept is usually non-zero.
#' @param b1 Slope coefficient for the independent variable (x).
#' @param b2 Slope coefficient for the moderator (w).
#' @param b3 Interaction coefficient between x and w.
#' @param mean_x Mean value of the independent variable (x).
#' @param var_x Variance of the independent variable (x).
#' @param mean_w Mean value of the moderator (w).
#' @param var_w Variance of the moderator (w).
#' @param legend_position Position of the legend. Options are: "UL" (inside upper left), "UR" (inside upper right), "BL" (inside bottom left), "BR" (inside bottom right), or "right" (default: outside center right).
#' @param fontsize Base font size for the plot.
#'
#' @details
#' The function creates six points (two for each of three moderator levels: M - 1SD, Mean, and M + 1SD) and connects them with lines.
#' You can change the title and axis names by modifying the \code{labs()} function within the code.
#' Keep in mind that the default intercept (b0) is set to 0 in Mplus for latent outcome identification.
#' In practice, especially in manifest models, the intercept should often be non-zero. Adjust b0 as needed.
#'
#' @return A ggplot2 object representing the simple slope plot.
#'
#' @examples
#' \dontrun{
#'   # Example usage:
#'   set.seed(123)
#'   p <- simple_slope(b0 = 1, b1 = 0.5, b2 = 0.3, b3 = 0.2,
#'                     mean_x = 0, var_x = 1,
#'                     mean_w = 0, var_w = 1,
#'                     legend_position = "UL", fontsize = 16)
#'   print(p)
#' }
#' @import ggplot2
#' @export
simple_slope <- function(b0 = 0, b1, b2, b3,
                         mean_x, var_x,
                         mean_w, var_w,
                         legend_position = "right",
                         fontsize = 12) {
  # Compute standard deviations
  sd_x <- sqrt(var_x)
  sd_w <- sqrt(var_w)

  # Define low and high values for the independent variable (x)
  x_low <- mean_x - sd_x
  x_high <- mean_x + sd_x

  # Define moderator levels: low, mean, and high (w)
  w_low  <- mean_w - sd_w
  w_mean <- mean_w
  w_high <- mean_w + sd_w

  # Create a data frame with six combinations:
  # Two x-values for each of the three w-levels.
  df <- data.frame(
    x = rep(c(x_low, x_high), times = 3),
    w = rep(c(w_low, w_mean, w_high), each = 2)
  )

  # Create a factor for the moderator levels with clear labels
  df$w_level <- factor(rep(c("M - 1SD", "Mean", "M + 1SD"), each = 2),
                       levels = c("M - 1SD", "Mean", "M + 1SD"))

  # Compute the predicted outcome using the moderated regression model:
  # y = b0 + b1*x + b2*w + b3*x*w
  df$y <- with(df, b0 + b1 * x + b2 * w + b3 * x * w)

  # Determine legend position and justification based on the new argument.
  if (legend_position == "UL") {
    leg_pos <- c(0.05, 0.95)
    leg_just <- c(0, 1)
  } else if (legend_position == "UR") {
    leg_pos <- c(0.95, 0.95)
    leg_just <- c(1, 1)
  } else if (legend_position == "BL") {
    leg_pos <- c(0.05, 0.05)
    leg_just <- c(0, 0)
  } else if (legend_position == "BR") {
    leg_pos <- c(0.95, 0.05)
    leg_just <- c(1, 0)
  } else {
    # Default outside at center right
    leg_pos <- "right"
    leg_just <- NULL
  }

  # Create the ggplot: each moderator level gets its own line connecting the two x-values.
  p <- ggplot(df, aes(x = x, y = y, group = w_level,
                      linetype = w_level, shape = w_level)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    labs(x = "Independent Variable (x)",
         y = "Predicted Outcome (y)",
         title = "Simple Slope Plot at Different Moderator Levels",
         linetype = "Moderator Level",
         shape = "Moderator Level") +
    scale_linetype_manual(values = c("M - 1SD" = "solid",
                                     "Mean"    = "dotted",
                                     "M + 1SD" = "solid")) +
    scale_shape_manual(values = c("M - 1SD" = 15,  # square
                                  "Mean"    = 16,  # triangle
                                  "M + 1SD" = 17)) +  # circle
    theme_classic(base_size = fontsize) +
    theme(
      legend.position = leg_pos,
      legend.justification = leg_just
    )

  return(p)
}
