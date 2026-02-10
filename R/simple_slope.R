#' Simple Slope Plot Based on MplusAutomation Output
#'
#' This function extracts model parameters from an MplusAutomation object and creates
#' a simple slope plot for a moderated regression model:
#'
#' \deqn{y = b_0 + b_1 x + b_2 w + b_3 (x \times w) + \varepsilon}
#'
#' It automatically retrieves regression coefficients, means, variances, and simple slopes
#' at different moderator levels (mean, mean - 1SD, mean + 1SD), and generates a plot
#' with annotated legend entries showing the slope estimates, standard errors, and p-values.
#'
#' The function supports both observed variable models (path analysis) and latent variable
#' models (SEM with latent interactions). It determines whether the predictor and moderator
#' are observed or latent by checking the sample statistics. For observed variables, means
#' and variances are extracted from sample statistics. For latent variables, means default
#' to 0 (as is standard in Mplus), and variances are extracted from the model parameter
#' estimates.
#'
#' @param model An MplusAutomation model object (typically output from
#'   \code{MplusAutomation::readModels()}).
#' @param x Character string. Name of the independent variable (predictor) as it appears
#'   in the Mplus output. Case-insensitive.
#' @param w Character string. Name of the moderator variable as it appears in the Mplus
#'   output. Case-insensitive.
#' @param int Character string. Name of the interaction term (x*w) as it appears in the
#'   Mplus output. For observed variable models, this is typically created using
#'   \code{DEFINE: int = x * w;} in the Mplus input. For latent variable models, this is
#'   defined using \code{int | x XWITH w;} in the MODEL command. Case-insensitive.
#' @param y Character string. Name of the dependent variable as it appears in the Mplus
#'   output. Can be observed or latent. Case-insensitive.
#' @param slp_l Character string. Name of the parameter for the simple slope at the low
#'   moderator level (M - 1SD), defined under \code{MODEL CONSTRAINT} in the Mplus input.
#'   Case-insensitive.
#' @param slp_h Character string. Name of the parameter for the simple slope at the high
#'   moderator level (M + 1SD), defined under \code{MODEL CONSTRAINT} in the Mplus input.
#'   Case-insensitive.
#' @param mean_x Numeric or \code{NULL}. Mean of the independent variable. If \code{NULL}
#'   (default), the function attempts to extract it from the sample statistics for observed
#'   variables, or defaults to 0 for latent variables.
#' @param mean_w Numeric or \code{NULL}. Mean of the moderator variable. Same extraction
#'   logic as \code{mean_x}.
#' @param var_x Numeric or \code{NULL}. Variance of the independent variable. If \code{NULL}
#'   (default), the function attempts to extract it from the sample statistics for observed
#'   variables, or from the model parameter estimates (Variances section) for latent variables.
#' @param var_w Numeric or \code{NULL}. Variance of the moderator variable. Same extraction
#'   logic as \code{var_x}.
#' @param legend_position Character string. Position of the legend. Options are:
#' \itemize{
#'   \item \code{"UL"}: Upper Left inside plot
#'   \item \code{"UR"}: Upper Right inside plot
#'   \item \code{"BL"}: Bottom Left inside plot
#'   \item \code{"BR"}: Bottom Right inside plot
#'   \item \code{"right"}: Default, outside plot on the right side
#' }
#' @param legend_name Character string. Title displayed above the legend entries.
#'   Default is \code{"Moderator Level"}.
#' @param fontsize Numeric. Base font size for the plot text. Default is 14.
#' @param return_data Logical. If \code{TRUE}, return the data frame used for plotting
#'   instead of the ggplot object. The returned data frame contains columns \code{x},
#'   \code{w}, and \code{y}. Default is \code{FALSE}.
#'
#' @details
#' The function performs the following steps:
#'
#' \strong{1. Coefficient extraction}
#'
#' The intercept of \code{y} is extracted from the \code{Intercepts} section. The
#' regression coefficients for \code{x}, \code{w}, and \code{int} are extracted from
#' the \code{y.ON} section of the unstandardized parameter estimates.
#'
#' \strong{2. Mean and variance extraction}
#'
#' When \code{mean_x}, \code{mean_w}, \code{var_x}, or \code{var_w} are left as
#' \code{NULL}, the function first checks whether the variable appears in
#' \code{model[["sampstat"]][["univariate.sample.statistics"]]}. If found, the variable
#' is treated as observed and its sample mean and variance are used. If not found, the
#' variable is treated as latent: its mean is set to 0 and its variance is read from
#' the \code{Variances} section of the model parameter estimates. Users can override
#' any of these values by passing explicit numeric arguments.
#'
#' \strong{3. Simple slope annotation}
#'
#' The simple slopes at M - 1SD and M + 1SD are retrieved from the
#' \code{New.Additional.Parameters} section. The simple slope at the mean is taken
#' from the main effect of \code{x}. These estimates, along with their standard errors
#' and p-values, are displayed in the legend. P-values are formatted without a leading
#' zero (e.g., \code{p = .057}).
#'
#' \strong{4. Plotting}
#'
#' Predicted values of \code{y} are computed at two levels of \code{x} (mean - 1SD and
#' mean + 1SD) crossed with three levels of \code{w} (mean - 1SD, mean, mean + 1SD),
#' producing six data points and three regression lines.
#'
#' \strong{Mplus input requirements}
#'
#' The Mplus input file should include \code{OUTPUT: SAMPSTAT;} so that sample statistics
#' are available for observed variable models. Simple slopes must be defined under
#' \code{MODEL CONSTRAINT} using the \code{NEW} statement. For example:
#'
#' \preformatted{
#' MODEL CONSTRAINT:
#'   NEW(slp_lo slp_hi);
#'   slp_lo = b1 + b3 * sqrt(var_w) * (-1);
#'   slp_hi = b1 + b3 * sqrt(var_w);
#' }
#'
#' @return If \code{return_data = FALSE} (default), a \code{ggplot2} object displaying the
#'   simple slope plot with three lines representing the moderator at M - 1SD, Mean, and
#'   M + 1SD. Each legend entry is annotated with the corresponding slope estimate,
#'   standard error, and p-value.
#'
#'   If \code{return_data = TRUE}, a data frame with three columns:
#'   \describe{
#'     \item{x}{Values of the independent variable (mean - 1SD and mean + 1SD).}
#'     \item{w}{Values of the moderator (mean + 1SD, mean, and mean - 1SD).}
#'     \item{y}{Predicted values of the dependent variable.}
#'   }
#'
#' @examples
#' \dontrun{
#' library(MplusAutomation)
#'
#' # --- Observed variable model (path analysis) ---
#' res <- readModels("path/to/output", filefilter = "manifest_mod")
#'
#' # Plot with legend in upper left corner
#' simple_slope(res,
#'              x = "parent", w = "res", int = "int", y = "se",
#'              slp_l = "slp_lo", slp_h = "slp_hi",
#'              legend_position = "UL")
#'
#' # --- Latent variable model (SEM with latent interaction) ---
#' res_sem <- readModels("path/to/output", filefilter = "sem_mod")
#'
#' # Plot with default legend position
#' simple_slope(res_sem,
#'              x = "parent", w = "res", int = "int", y = "atl",
#'              slp_l = "slp_lo", slp_h = "slp_hi")
#'
#' # Return plotting data instead of the figure
#' plot_df <- simple_slope(res_sem,
#'                         x = "parent", w = "res", int = "int", y = "atl",
#'                         slp_l = "slp_lo", slp_h = "slp_hi",
#'                         return_data = TRUE)
#'
#' # Override automatically extracted values
#' simple_slope(res_sem,
#'              x = "parent", w = "res", int = "int", y = "atl",
#'              slp_l = "slp_lo", slp_h = "slp_hi",
#'              mean_x = 0, mean_w = 0,
#'              var_x = 0.357, var_w = 0.386)
#' }
#'
#' @import ggplot2
#' @export
simple_slope <- function(model, x, w, int, y, slp_l, slp_h,
                         mean_x = NULL, mean_w = NULL,
                         var_x = NULL, var_w = NULL,
                         legend_position = "right",
                         legend_name = "Moderator Level",
                         fontsize = 14,
                         return_data = FALSE) {
  # Retrieve the unstandardized parameters data.frame from the model
  params <- model[["parameters"]][["unstandardized"]]

  # Retrieve sample statistics (means and variances)
  sampstat <- model[["sampstat"]][["univariate.sample.statistics"]]

  # Trim whitespace from rownames to avoid matching failures
  if (!is.null(sampstat)) {
    rownames(sampstat) <- trimws(rownames(sampstat))
  }

  # --- Helper: extract mean and variance for a variable ---
  # If found in sampstat, treat as observed variable and read from there.
  # If not found in sampstat, treat as latent variable:
  #   mean defaults to 0, variance is read from parameters$unstandardized.
  extract_mv <- function(varname, user_mean, user_var) {
    vname_upper <- toupper(varname)
    in_sampstat <- !is.null(sampstat) && vname_upper %in% rownames(sampstat)

    # Mean
    if (!is.null(user_mean)) {
      m <- user_mean
    } else if (in_sampstat) {
      m <- sampstat[vname_upper, "Mean"]
    } else {
      m <- 0
    }

    # Variance
    if (!is.null(user_var)) {
      v <- user_var
    } else if (in_sampstat) {
      v <- sampstat[vname_upper, "Variance"]
    } else {
      var_row <- subset(params, (paramHeader == "Variances" & param == vname_upper))
      if (nrow(var_row) > 0) {
        v <- var_row$est[1]
      } else {
        v <- 1
        warning("Variance of '", varname, "' not found in model output. Defaulting to 1.")
      }
    }

    list(mean = m, var = v)
  }

  mv_x <- extract_mv(x, mean_x, var_x)
  mv_w <- extract_mv(w, mean_w, var_w)
  mean_x <- mv_x$mean
  mean_w <- mv_w$mean
  var_x  <- mv_x$var
  var_w  <- mv_w$var

  # --- Extract regression coefficients ---
  intercept_row <- subset(params, (tolower(paramHeader) == "intercepts" & param == toupper(y)))
  b0 <- if (nrow(intercept_row) > 0) intercept_row$est[1] else 0

  iv_row <- subset(params, (paramHeader == paste0(toupper(y), ".ON") & param == toupper(x)))
  b1 <- if (nrow(iv_row) > 0) iv_row$est[1] else 0

  mod_row <- subset(params, (paramHeader == paste0(toupper(y), ".ON") & param == toupper(w)))
  b2 <- if (nrow(mod_row) > 0) mod_row$est[1] else 0

  intx_row <- subset(params, (paramHeader == paste0(toupper(y), ".ON") & param == toupper(int)))
  b3 <- if (nrow(intx_row) > 0) intx_row$est[1] else 0

  # --- Helper functions for formatting ---
  num_fmt <- function(x, d) sprintf(paste0("%0.", d, "f"), x)

  p_fmt <- function(p) {
    if (p < .001) {
      "p < .001"
    } else {
      paste0("p = ", sub("^0\\.", ".", num_fmt(p, 3)))
    }
  }

  # --- Retrieve simple slope labels ---
  slp_l_row <- subset(params, (paramHeader == "New.Additional.Parameters" & param == toupper(slp_l)))
  slp_h_row <- subset(params, (paramHeader == "New.Additional.Parameters" & param == toupper(slp_h)))

  slp_lo   <- if (nrow(slp_l_row) > 0)
    paste0("b = ", num_fmt(slp_l_row$est[1], 3),
           ", SE = ", num_fmt(slp_l_row$se[1], 3),
           ", ", p_fmt(slp_l_row$pval[1])) else ""

  slp_mean <- if (nrow(iv_row) > 0)
    paste0("b = ", num_fmt(iv_row$est[1], 3),
           ", SE = ", num_fmt(iv_row$se[1], 3),
           ", ", p_fmt(iv_row$pval[1])) else ""

  slp_hi   <- if (nrow(slp_h_row) > 0)
    paste0("b = ", num_fmt(slp_h_row$est[1], 3),
           ", SE = ", num_fmt(slp_h_row$se[1], 3),
           ", ", p_fmt(slp_h_row$pval[1])) else ""

  # Compute standard deviations
  sd_x <- sqrt(var_x)
  sd_w <- sqrt(var_w)

  # --- Create data for plotting ---
  x_low <- mean_x - sd_x
  x_high <- mean_x + sd_x

  w_low  <- mean_w - sd_w
  w_mean <- mean_w
  w_high <- mean_w + sd_w

  df <- data.frame(
    x = rep(c(x_low, x_high), times = 3),
    w = rep(c(w_high, w_mean, w_low), each = 2)
  )

  # Compute predicted outcome using: y = b0 + b1*x + b2*w + b3*x*w
  df$y <- with(df, b0 + b1 * x + b2 * w + b3 * x * w)

  # --- Return data if requested ---
  if (return_data) {
    return(df)
  }

  # --- Build w_level factor for legend ---
  w_low_lable  <- paste0("M - 1SD", "\n", slp_lo)
  w_mean_lable <- paste0("Mean", "\n", slp_mean)
  w_high_lable <- paste0("M + 1SD", "\n", slp_hi)

  df$w_level <- factor(rep(c(w_high_lable, w_mean_lable, w_low_lable), each = 2),
                       levels = c(w_high_lable, w_mean_lable, w_low_lable))

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
    leg_pos <- "right"
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
          legend.key.size = unit(1, "cm"),
          legend.key.height = unit(2, "lines"),
          legend.justification = leg_just,
          legend.background = element_rect(fill = "white", color = "black"))

  return(p)
}
