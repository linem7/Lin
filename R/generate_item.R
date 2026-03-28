#' Generate Likert-Scale Item Scores from Mean Scores
#'
#' Given a numeric column of mean scores (typically from an existing dataset),
#' this function reverse-engineers a set of Likert-scale item-level responses
#' that match the implied total scores while preserving reasonable psychometric
#' properties (unidimensional factor structure, adequate reliability, and
#' approximately normal marginal distributions).
#'
#' The generation process involves three stages. First, a population-level
#' covariance structure is defined via a single-factor model with specified
#' factor loadings, ordinal thresholds, and optional residual covariances.
#' Second, a large candidate pool is simulated from this model using
#' \code{\link[lavaan]{simulateData}}. Third, for each input mean score, a
#' row whose total score matches (or most closely approximates) the target
#' is sampled from the pool and returned.
#'
#' @param data A data frame containing the column specified by \code{var}.
#' @param var Unquoted column name in \code{data} that contains the mean
#'   scores to be matched. Must be numeric.
#' @param n_items Integer. Number of items to generate. Default is \code{8}.
#' @param item_name Character string used as the prefix for generated item
#'   column names. If \code{NULL} (default), the name of \code{var} is used.
#'   Item columns are named \code{<item_name>1}, \code{<item_name>2}, etc.
#' @param likert_points Integer. Number of response categories on the Likert
#'   scale (e.g., \code{5} for a 1--5 scale). Default is \code{5}.
#' @param loading_range Numeric vector of length 2 giving the lower and upper
#'   bounds for uniformly sampled factor loadings. Higher loadings produce
#'   stronger inter-item correlations and higher reliability. Default is
#'   \code{c(0.65, 0.80)}.
#' @param concentration Numeric. Controls the spacing of ordinal thresholds.
#'   Larger values compress the middle categories and spread the tails,
#'   yielding a more peaked (leptokurtic) marginal distribution. Default is
#'   \code{1.5}.
#' @param threshold_shift Numeric. A constant added to all thresholds after
#'   scaling, allowing the entire response distribution to be shifted upward
#'   or downward. Default is \code{0}.
#' @param n_resid_cov Integer. Number of randomly selected item pairs for
#'   which a residual covariance is introduced. This produces realistic,
#'   slightly imperfect model fit. Set to \code{0} to disable. Default is
#'   \code{5}.
#' @param resid_cov_range Numeric vector of length 2 giving the lower and
#'   upper bounds for uniformly sampled residual covariance values. Default
#'   is \code{c(0.08, 0.15)}.
#' @param pool_size Integer. Number of candidate rows to simulate. Larger
#'   pools improve the chance of exact total-score matches but consume more
#'   memory. Default is \code{10000}.
#' @param check Logical. If \code{TRUE}, descriptive statistics, CFA fit
#'   indices, and reliability estimates are printed for the final matched
#'   subset (not the full candidate pool). Requires the \pkg{psych} package
#'   and \code{Lin::fit_table_lavaan} / \code{Lin::reliability_table}.
#'   Default is \code{FALSE}.
#'
#' @return A tibble consisting of all columns in \code{data} with
#'   \code{n_items} new item columns appended on the right. Each item column
#'   contains integer values from \code{1} to \code{likert_points}.
#'
#' @details
#' When the candidate pool does not contain any row with the exact target
#' total score, the row with the closest total score is selected instead.
#' Increasing \code{pool_size} reduces the frequency of such approximations.
#'
#' Because the matched subset is drawn with replacement from a finite pool,
#' duplicate response patterns may occur, especially when \code{pool_size}
#' is small relative to the number of rows in \code{data}.
#'
#' @examples
#' # Example 1: Basic usage with default settings
#' set.seed(42)
#' data <- tibble::tibble(
#'   id = 1:200,
#'   x = round(runif(200, 1.5, 4.5), 2)
#' )
#' result <- generate_item(data, x)
#' head(result)
#'
#' # Example 2: Custom item naming, 7-point scale, and quality check
#' set.seed(123)
#' data2 <- tibble::tibble(
#'   id = 1:300,
#'   satisfaction = round(runif(300, 2.0, 6.5), 2)
#' )
#' result2 <- generate_item(
#'   data2,
#'   satisfaction,
#'   n_items = 6,
#'   item_name = "sat",
#'   likert_points = 7,
#'   loading_range = c(0.70, 0.85),
#'   check = TRUE
#' )
#' head(result2)
#'
#' @seealso \code{\link[lavaan]{simulateData}} for the underlying simulation
#'   engine.
#'
#' @export
generate_item <- function(
    data,
    var,
    n_items = 8,
    item_name = NULL,
    likert_points = 5,
    loading_range = c(0.65, 0.80),
    concentration = 1.5,
    threshold_shift = 0,
    n_resid_cov = 5,
    resid_cov_range = c(0.08, 0.15),
    pool_size = 10000,
    check = FALSE
) {

  # --- Parse column name ---
  var_name <- deparse(substitute(var))
  mean_scores <- data[[var_name]]

  if (is.null(mean_scores)) stop(sprintf("Column '%s' not found in data.", var_name))
  if (!is.numeric(mean_scores)) stop(sprintf("Column '%s' must be numeric.", var_name))

  # --- Item naming: default to the original column name ---
  if (is.null(item_name)) item_name <- var_name
  item_names <- paste0(item_name, seq_len(n_items))

  midpoint <- (1 + likert_points) / 2

  # --- Factor loadings ---
  lambdas <- runif(n_items, loading_range[1], loading_range[2])
  loading_line <- paste0(
    "f =~ ",
    paste(sprintf("%.3f*%s", lambdas, item_names), collapse = " + ")
  )

  # --- Thresholds ---
  base_thresholds <- qnorm(seq_len(likert_points - 1) / likert_points)
  anchor <- median(mean_scores, na.rm = TRUE)
  shift <- -(anchor - midpoint) * concentration * 0.5 + threshold_shift

  threshold_lines <- vapply(item_names, function(item) {
    thresholds <- base_thresholds * concentration + shift
    terms <- paste(sprintf("%.3f*t%d", thresholds, seq_along(thresholds)), collapse = " + ")
    sprintf("%s | %s", item, terms)
  }, character(1))

  # --- Residual covariances ---
  resid_lines <- character(0)
  if (n_items >= 2 && n_resid_cov > 0) {
    all_pairs <- t(combn(item_names, 2))
    n_resid_cov <- min(n_resid_cov, nrow(all_pairs))
    sel <- sample(seq_len(nrow(all_pairs)), size = n_resid_cov, replace = FALSE)
    sel_pairs <- all_pairs[sel, , drop = FALSE]
    resid_vals <- runif(n_resid_cov, resid_cov_range[1], resid_cov_range[2])
    resid_lines <- vapply(seq_len(n_resid_cov), function(i) {
      sprintf("%s ~~ %.3f*%s", sel_pairs[i, 1], resid_vals[i], sel_pairs[i, 2])
    }, character(1))
  }

  # --- Assemble model syntax ---
  model_syn <- paste(
    loading_line,
    paste(threshold_lines, collapse = "\n"),
    if (length(resid_lines) > 0) paste(resid_lines, collapse = "\n") else NULL,
    sep = "\n\n"
  )

  # --- Simulate candidate pool ---
  sim_data <- lavaan::simulateData(model_syn, sample.nobs = pool_size)
  sim_mat <- as.matrix(sim_data) + 1L
  colnames(sim_mat) <- item_names
  row_sums <- rowSums(sim_mat)

  # --- Index by total score for fast lookup ---
  sum_index <- split(seq_len(pool_size), row_sums)

  # --- Match each target total score ---
  target_sums <- round(mean_scores * n_items)

  chosen_rows <- vapply(target_sums, function(ts) {
    ts_char <- as.character(ts)
    candidates <- sum_index[[ts_char]]
    if (!is.null(candidates) && length(candidates) > 0) {
      idx <- sample(candidates, 1)
    } else {
      idx <- which.min(abs(row_sums - ts))
    }
    idx
  }, integer(1))

  items_df <- tibble::as_tibble(sim_mat[chosen_rows, , drop = FALSE])

  # --- Optional quality check (on matched subset) ---
  if (check) {
    check_data <- as.data.frame(items_df)
    cfa_syntax <- paste0("f =~ ", paste(item_names, collapse = " + "))

    cat("=== Descriptives (matched subset) ===\n")
    desc <- psych::describe(check_data) %>%
      dplyr::select(mean, sd, min, max, skew, kurtosis)
    print(desc)

    cat("\n=== CFA Fit (matched subset) ===\n")
    print(Lin::fit_table_lavaan(cfa_syntax, data = check_data))

    cat("\n=== Reliability (matched subset) ===\n")
    print(Lin::reliability_table(cfa_syntax, data = check_data))
  }

  # --- Bind generated items to original data ---
  result <- dplyr::bind_cols(data, items_df)

  return(result)
}
