#' Batch-convert multiple Z-score variables to Likert scales
#'
#' Applies \code{\link{likertize_zscore}} across groups of variables in a data frame,
#' using per-group settings for the number of Likert levels and optional skew.
#' The function does not re-estimate densities from the data; it only uses the chosen
#' reference density to define proportionally balanced cutpoints and then bins.
#'
#' @param data A \code{data.frame} (or tibble) containing Z-score variables to be discretized.
#'   Non-target columns are copied through unchanged.
#' @param var_config A list of configuration lists. Each configuration must contain:
#'   \itemize{
#'     \item \code{vars}: character vector of column names in \code{data} to convert.
#'     \item \code{n_levels}: integer (>= 2), number of Likert levels for these \code{vars}.
#'     \item \code{skew}: numeric, passed to \code{\link{likertize_zscore}}. Here \code{skew}
#'       is the skew-normal shape parameter used to define the reference density (\code{0} for symmetric).
#'   }
#'   Example:
#'   \preformatted{
#'   var_config <- list(
#'     list(vars = c("z_a", "z_b"), n_levels = 5, skew = 0),
#'     list(vars = "z_c",           n_levels = 7, skew = 2)
#'   )
#'   }
#' @param as_ordered Logical. If \code{TRUE}, converts the resulting integer categories to
#'   \code{ordered} factors with fixed levels \code{1:n_levels} for each variable, based on its
#'   configuration. Defaults to \code{FALSE}.
#'
#' @return A data frame of the same class as \code{data}, with the specified variables
#'   replaced by their Likert-coded versions (integer or ordered factor depending on \code{as_ordered}).
#'
#' @details
#' - Variables listed in \code{var_config} that are not found in \code{data} trigger a warning and are skipped.
#' - If a target variable contains \code{NA}, they are preserved.
#' - If different groups in \code{var_config} reference the same variable name, the last occurrence wins.
#' - When \code{as_ordered = TRUE}, factor levels are fixed to \code{1:n_levels} from the final config for that variable.
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   z_a = rnorm(1000),
#'   z_b = rnorm(1000, mean = 0.3),
#'   z_c = rnorm(1000)
#' )
#'
#' cfg <- list(
#'   list(vars = c("z_a", "z_b"), n_levels = 5, skew = 0),
#'   list(vars = "z_c",           n_levels = 7, skew = 2)
#' )
#'
#' out <- convert_to_likert_batch(df, cfg)
#' head(out)
#'
#' out_ord <- convert_to_likert_batch(df, cfg, as_ordered = TRUE)
#' head(out_ord)
#'
#' @seealso \code{\link{likertize_zscore}}, \code{\link{discretize_density}}
#' @export
convert_to_likert_batch <- function(data, var_config, as_ordered = FALSE) {

  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame (or tibble).")
  }
  if (!is.list(var_config) || length(var_config) == 0L) {
    stop("`var_config` must be a non-empty list of configuration lists.")
  }
  if (!is.logical(as_ordered) || length(as_ordered) != 1L) {
    stop("`as_ordered` must be a single logical value.")
  }

  data_likert <- data
  all_vars <- character(0)

  # record the final n_levels per variable (last occurrence wins)
  var_levels <- integer(0)

  for (i in seq_along(var_config)) {
    config <- var_config[[i]]

    if (!is.list(config) ||
        !"vars" %in% names(config) ||
        !"n_levels" %in% names(config) ||
        !"skew" %in% names(config)) {
      stop(sprintf(
        "Each config must be a list with `vars`, `n_levels`, and `skew`. Problem at index %d.",
        i
      ))
    }

    vars     <- config$vars
    n_levels <- config$n_levels
    skew     <- config$skew

    if (!is.character(vars) || length(vars) == 0L) {
      stop(sprintf("`vars` must be a non-empty character vector (config index %d).", i))
    }
    if (length(n_levels) != 1L || !is.finite(n_levels) || n_levels < 2) {
      stop(sprintf("`n_levels` must be a single finite integer >= 2 (config index %d).", i))
    }
    n_levels <- as.integer(n_levels)

    if (length(skew) != 1L || !is.finite(skew)) {
      stop(sprintf("`skew` must be a single finite numeric value (config index %d).", i))
    }

    for (var in vars) {
      if (var %in% names(data)) {
        data_likert[[var]] <- likertize_zscore(
          z_score  = data[[var]],
          n_levels = n_levels,
          skew     = skew
        )
        var_levels[[var]] <- n_levels
      } else {
        warning(sprintf("Variable '%s' not found in `data`; skipping.", var))
      }
    }

    all_vars <- c(all_vars, vars)
  }

  if (isTRUE(as_ordered)) {
    unique_vars <- unique(all_vars)
    for (var in unique_vars) {
      if (var %in% names(data_likert)) {
        nl <- var_levels[[var]]
        if (is.null(nl) || !is.finite(nl) || nl < 2) next
        data_likert[[var]] <- ordered(data_likert[[var]], levels = seq_len(nl))
      }
    }
  }

  data_likert
}
