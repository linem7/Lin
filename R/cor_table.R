#' Correlation Table with Descriptive Statistics
#'
#' @description
#' Computes pairwise Pearson correlations for a selected set of numeric variables,
#' appends significance stars, and combines these with means (M) and standard
#' deviations (SD) in an APA-style table.  Descriptive statistics can appear
#' either in the first two columns ("left") or as the last two rows ("bottom").
#'
#' @param data A data frame (or tibble) containing the variables to analyze.
#' @param vars Unquoted variables to include; passed via tidyselect (e.g., `c(mpg, hp, wt)`).
#' @param digits_stat Integer; number of decimal places for M and SD (default: 2).
#' @param digits_corr Integer; number of decimal places for correlation coefficients (default: 3).
#' @param stats_position Character; where to place M and SD.
#'        - `"left"`: in columns after Variable (default)
#'        - `"bottom"`: as extra rows below the correlation grid
#'
#' @return A data frame with columns:
#'   * `Variable`: numbered variable names
#'   * `M`, `SD`: descriptive statistics (when `stats_position = "left"`)
#'   * `1`, `2`, …: formatted correlations with significance stars
#'   * _or_ additional rows (`"M"`, `"SD"`) when `stats_position = "bottom"`.
#'
#' @seealso \code{\link[stats]{cor.test}}
#' @export
#' @importFrom dplyr select
#'
#' @examples
#' library(datasets)
#' data(mtcars)
#'
#' # Default: stats on the left
#' cor_table(mtcars, c(mpg, disp, hp, wt))
#'
#' # Stats on the bottom
#' cor_table(
#'   mtcars,
#'   c(mpg, disp, hp, wt),
#'   stats_position = "bottom"
#' )
#'
#' # Using iris numeric columns
#' data(iris)
#' cor_table(
#'   iris,
#'   c(Sepal.Length, Sepal.Width, Petal.Length),
#'   digits_stat = 1,
#'   digits_corr = 2
#' )
cor_table <- function(data, vars,
                      digits_stat    = 2,
                      digits_corr    = 3,
                      stats_position = c("left", "bottom")) {

  stats_position <- match.arg(stats_position)

  ## correlation formatter (APA 7) ----------------------------------------
  fmt_r <- function(x, digits) {
    sub("^(-?)0\\.", "\\1.", sprintf(paste0("%.", digits, "f"), x))
  }
  ## ----------------------------------------------------------------------

  # select and name variables
  data_sel  <- select(data, {{ vars }})
  var_names <- names(data_sel)
  n         <- length(var_names)

  # helper for significance stars
  get_stars <- function(p_val) {
    if (is.na(p_val))    return("")
    if (p_val < 0.001)   return("***")
    if (p_val < 0.01)    return("**")
    if (p_val < 0.05)    return("*")
    ""
  }

  # descriptive stats
  means <- sapply(data_sel, mean, na.rm = TRUE)
  sds   <- sapply(data_sel, sd,   na.rm = TRUE)

  # compute pairwise cor and p-value (lower triangle)
  corr_matrix <- matrix(NA, n, n)
  p_matrix    <- matrix(NA, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i > j) {
        tst <- cor.test(data_sel[[i]], data_sel[[j]])
        corr_matrix[i, j] <- tst$estimate
        p_matrix[i, j]    <- tst$p.value
      }
    }
  }

  # build character matrix
  formatted_corr <- matrix("", n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i > j) {
        r <- fmt_r(corr_matrix[i, j], digits_corr)  # CHANGED
        s <- get_stars(p_matrix[i, j])
        formatted_corr[i, j] <- paste0(r, s)
      } else if (i == j) {
        formatted_corr[i, j] <- "\u2013"             # diagonal
      }
    }
  }

  # assemble final table
  if (stats_position == "left") {
    corr_df <- as.data.frame(formatted_corr, stringsAsFactors = FALSE)
    colnames(corr_df) <- as.character(seq_len(n))

    final_df <- data.frame(
      Variable = paste0(seq_len(n), ". ", var_names),
      M  = sprintf(paste0("%.", digits_stat, "f"), means),  # zero kept
      SD = sprintf(paste0("%.", digits_stat, "f"), sds),    # zero kept
      corr_df,
      check.names = FALSE, stringsAsFactors = FALSE
    )

  } else {  # stats on bottom
    tbl <- matrix("", n + 2, n + 1)
    colnames(tbl) <- c("Variable", as.character(seq_len(n)))

    for (i in seq_len(n)) {
      tbl[i, 1] <- paste0(i, ". ", var_names[i])
      tbl[i, 2:(n + 1)] <- formatted_corr[i, ]
    }
    tbl[n + 1, 1] <- "M"
    tbl[n + 2, 1] <- "SD"
    for (j in seq_len(n)) {
      tbl[n + 1, j + 1] <- sprintf(paste0("%.", digits_stat, "f"), means[j])
      tbl[n + 2, j + 1] <- sprintf(paste0("%.", digits_stat, "f"), sds[j])
    }

    final_df <- as.data.frame(tbl, stringsAsFactors = FALSE, check.names = FALSE)
  }

  return(final_df)
}

