#' Simple Slope Plot Based on MplusAutomation Output
#'
#' This function extracts model parameters from an MplusAutomation object and creates
#' a simple slope plot for a moderated regression model:
#'
#' \deqn{y = b0 + b1*x + b2*w + b3*(x*w) + \varepsilon}
#'
#' It automatically retrieves regression coefficients, variances, and simple slopes
#' at different moderator levels (mean, mean - 1SD, mean + 1SD), and generates a plot
#' with annotated legend entries showing the slope estimates, standard errors, and p-values.
#'
#' @param model An MplusAutomation model object containing parameters (typically output from `readModels()`).
#' @param x Name of the independent variable in Mplus output.
#' @param w Name of the moderator variable.
#' @param int Name of the interaction term (x*w).
#' @param y Name of the dependent variable.
#' @param slp_l Name of the parameter for the simple slope at low moderator value under model constraint (M - 1SD).
#' @param slp_h Name of the parameter for the simple slope at high moderator value  under model constraint (M + 1SD).
#' @param mean_x Mean value of the independent variable (default = 0).
#' @param mean_w Mean value of the moderator (default = 0).
#' @param legend_position Position of the legend. Options are:
#' \itemize{
#'   \item `"UL"` = Upper Left inside plot,
#'   \item `"UR"` = Upper Right inside plot,
#'   \item `"BL"` = Bottom Left inside plot,
#'   \item `"BR"` = Bottom Right inside plot,
#'   \item `"right"` = Default, outside plot on right side.
#' }
#' @param fontsize Font size for the plot text (default = 14).
#'
#' @details
#' The function extracts:
#' \itemize{
#'   \item The intercept from the \code{Intercepts} section.
#'   \item The coefficients for the independent variable (x), moderator (w), and interaction (x*w) from the \code{y.ON} section.
#'   \item The variances of x and w from the \code{Variances} section.
#'   \item The simple slopes at M - 1SD and M + 1SD from the \code{New.Additional.Parameters} section.
#' }
#'
#' The legend labels dynamically include the simple slope estimates, standard errors, and p-values, separated by a new line (`\\n`).
#'
#' @return A ggplot2 object displaying the simple slope plot.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' simple_slope(model,
#'              x = "x",
#'              w = "w",
#'              int = "int",
#'              y = "y",
#'              slp_h = "slp_hi",
#'              slp_l = "slp_lo",
#'              legend_position = "UL")
#' }
#'
#' @import ggplot2
#' @export

simple_slope <- function(model, x, w, int, y, slp_l, slp_h,
                         mean_x = 0, mean_w = 0,
                         legend_position = "right",
                         legend_name = "Moderator Level",
                         fontsize = 14) {
  # Retrieve the unstandardized parameters data.frame from the model
  params <- model[["parameters"]][["unstandardized"]]

  # --- Extract regression coefficients ---
  # Intercept: if found, use the "INT" row; otherwise default to 0
  intercept_row <- subset(params, (paramHeader == "intercept" & param == toupper(y)))
  b0 <- if (nrow(intercept_row) > 0) intercept_row$est[1] else 0

  # Independent variable (iv) coefficient: search for "x" (default to 0)
  iv_row <- subset(params, (paramHeader == paste0(toupper(y), ".ON") & param == toupper(x)))
  b1 <- if (nrow(iv_row) > 0) iv_row$est[1] else 0


  # Moderator coefficient: search for "w" (default to 0)
  mod_row <- subset(params, (paramHeader == paste0(toupper(y), ".ON") & param == toupper(w)))
  b2 <- if (nrow(mod_row) > 0) mod_row$est[1] else 0

  # Interaction coefficient: search for "int" (default to 0)
  intx_row <- subset(params, (paramHeader == paste0(toupper(y), ".ON") & param == toupper(int)))
  b3 <- if (nrow(intx_row) > 0) intx_row$est[1] else 0

  # --- Extract variances for x and w ---
  varx_row <- subset(params, (paramHeader == "Variances" & param == toupper(x)))
  varw_row <- subset(params, (paramHeader == "Variances" & param == toupper(w)))
  var_x <- if(nrow(varx_row) > 0) varx_row$est[1] else 1
  var_w <- if(nrow(varw_row) > 0) varw_row$est[1] else 1

  # --- Retrieve simple slope labels ---
  slp_l_row <- subset(params, (paramHeader == "New.Additional.Parameters" & param == toupper(slp_l)))
  slp_h_row <- subset(params, (paramHeader == "New.Additional.Parameters" & param == toupper(slp_h)))

  slp_lo <- if (nrow(slp_l_row) > 0) paste0("b = ", sprintf("%.3f", slp_l_row$est[1]), ", SE = ", sprintf("%.3f", slp_l_row$se[1]), ", p = ", sprintf("%.3f", slp_l_row$pval[1])) else ""
  slp_mean <- if (nrow(iv_row) > 0) paste0("b = ", sprintf("%.3f", iv_row$est[1]), ", SE = ", sprintf("%.3f", iv_row$se[1]), ", p = ", sprintf("%.3f", iv_row$pval[1])) else ""
  slp_hi <- if (nrow(slp_h_row) > 0) paste0("b = ", sprintf("%.3f", slp_h_row$est[1]), ", SE = ", sprintf("%.3f", slp_h_row$se[1]), ", p = ", sprintf("%.3f", slp_h_row$pval[1])) else ""

  # Compute standard deviations
  sd_x <- sqrt(var_x)
  sd_w <- sqrt(var_w)

  # --- Create data for plotting ---
  # Define two x values: mean ± SD (using provided mean_x)
  x_low <- mean_x - sd_x
  x_high <- mean_x + sd_x

  # Define three moderator levels: one SD below, at, and one SD above (using provided mean_w)
  w_low  <- mean_w - sd_w
  w_mean <- mean_w
  w_high <- mean_w + sd_w

  # Create a data frame: two x-values for each of the three moderator levels
  df <- data.frame(
    x = rep(c(x_low, x_high), times = 3),
    w = rep(c(w_low, w_mean, w_high), each = 2)
  )

  w_low_lable = paste0("M - 1SD", "\n", slp_lo)
  w_mean_lable = paste0("Mean", "\n", slp_mean)
  w_high_lable = paste0("M + 1SD", "\n", slp_hi)

  df$w_level <- factor(rep(c(w_low_lable, w_mean_lable, w_high_lable), each = 2),
                       levels = c(w_low_lable, w_mean_lable, w_high_lable))

  # Compute predicted outcome using: y = b0 + b1*x + b2*w + b3*x*w
  df$y <- with(df, b0 + b1 * x + b2 * w + b3 * x * w)

  # --- Legend positioning ---
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
    leg_pos <- "right"  # default: outside center right
    leg_just <- NULL
  }

  # --- Create ggplot ---
  p <- ggplot(df, aes(x = x, y = y, group = w_level,
                      linetype = w_level, shape = w_level)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    labs(x = "Independent Variable (x)",
         y = "Predicted Outcome (y)",
         title = "Simple Slope Plot at Different Moderator Levels",
         linetype = legend_name,
         shape = legend_name) +
    theme_classic(base_size = fontsize) +
    theme(text = element_text(family = "serif"),
          legend.position = leg_pos,
          legend.key.height = unit(2, "lines"),
          legend.justification = leg_just,
          legend.background = element_rect(fill = "white", color = "black"))

  return(p)
}


