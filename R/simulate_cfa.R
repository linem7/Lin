#' Simulate Likert-type data for confirmatory factor analysis
#'
#' Generates ordinal (Likert-scale) data from a specified multi-factor CFA
#' population model using \code{lavaan::simulateData()}. The function provides
#' independent control over factor loadings, factor correlations, threshold
#' spacing (distribution shape), and threshold shift (distribution location),
#' along with optional within-factor residual covariances.
#'
#' @param n_obs Integer. Number of observations to simulate. Default is 300.
#' @param loading_range Numeric vector of length 2. The lower and upper bounds
#'   from which standardised factor loadings are drawn uniformly for each item.
#'   Default is \code{c(0.68, 0.80)}.
#' @param factor_cor Factor correlations. Accepts four input forms:
#'   \describe{
#'     \item{Scalar}{A single numeric value in (-1, 1), broadcast to all factor
#'       pairs.}
#'     \item{Unnamed vector}{A numeric vector of length equal to the number of
#'       unique factor pairs, assigned in row-major order of the upper triangle
#'       (i.e., pair 1-2, 1-3, ..., 2-3, 2-4, ...).}
#'     \item{Named vector}{Names should follow the format
#'       \code{"factorA_factorB"}, matching factor names in \code{factor_items}.}
#'     \item{Matrix}{A symmetric numeric matrix whose dimensions match the
#'       number of factors. Upper-triangle values are used.}
#'   }
#'   Default is 0.25.
#' @param threshold_shift Numeric scalar or named numeric vector. Additive
#'   shift applied to all thresholds of items within each factor, after the
#'   \code{concentration} scaling. Negative values shift the observed
#'   distribution toward higher response categories (higher mean); positive
#'   values shift it toward lower categories (lower mean). When a single
#'   unnamed value is supplied, it is recycled to all factors. Default is
#'   \code{c(a = -0.10, b = -0.10, c = -0.10, d = -0.10, e = -0.10)}.
#' @param concentration Positive numeric scalar. Multiplicative scaling factor
#'   applied to the baseline equal-probability thresholds
#'   (\code{qnorm(1:(\emph{k}-1) / \emph{k})}). Larger values push thresholds
#'   further from zero, widening the central response category and producing a
#'   more bell-shaped (peaked) distribution. Smaller values compress thresholds
#'   toward zero, yielding a flatter, more uniform distribution. Default is 1.5.
#' @param likert_points Integer >= 2. Number of response categories on the
#'   Likert scale. Default is 5.
#' @param n_resid_cov Non-negative integer. Number of within-factor residual
#'   covariances to add to the population model. Pairs are sampled randomly
#'   from items belonging to the same factor. If the requested number exceeds
#'   the available within-factor pairs, it is silently truncated. Default is 5.
#' @param factor_items Named list of character vectors. Each element defines the
#'   item names belonging to one factor. List names serve as factor names.
#'   Default provides five factors (a through e) with 8, 4, 5, 6, and 6 items.
#' @param resid_cov_range Numeric vector of length 2. The lower and upper
#'   bounds from which residual covariance magnitudes are drawn uniformly.
#'   Default is \code{c(0.08, 0.15)}.
#'
#' @return A \code{\link[tibble]{tibble}} of simulated responses. All item
#'   columns are numeric, with values ranging from 0 to
#'   \code{likert_points - 1}.
#'
#' @details
#' The data-generating process follows the latent response variable framework.
#' Continuous latent responses are drawn from a multivariate normal distribution
#' whose covariance structure is implied by the factor loadings and factor
#' correlations. These continuous values are then discretised into ordered
#' categories according to the threshold parameters.
#'
#' The threshold for item \eqn{j} at cut-point \eqn{c} is computed as:
#' \deqn{\tau_{jc} = \texttt{concentration} \times \Phi^{-1}(c / K) +
#' \texttt{threshold\_shift}_f}
#' where \eqn{K} is \code{likert_points}, \eqn{\Phi^{-1}} is the standard
#' normal quantile function, and \eqn{f} indexes the factor to which item
#' \eqn{j} belongs.
#'
#' @examples
#' library(lavaan)
#' library(dplyr)
#'
#' # --- Example 1: Scalar factor_cor (two factors) ---
#' set.seed(2024)
#' dat1 <- simulate_cfa(
#'   n_obs           = 500,
#'   loading_range   = c(0.65, 0.80),
#'   factor_cor      = 0.30,
#'   threshold_shift = c(F1 = -0.25, F2 = 0.00),
#'   concentration   = 1.5,
#'   likert_points   = 5,
#'   n_resid_cov     = 2,
#'   factor_items    = list(
#'     F1 = paste0("x", 1:4),
#'     F2 = paste0("y", 1:3)
#'   ),
#'   resid_cov_range = c(0.05, 0.12)
#' )
#' str(dat1)
#'
#' # --- Example 2: Named vector factor_cor (three factors) ---
#' set.seed(100)
#' dat2 <- simulate_cfa(
#'   n_obs           = 400,
#'   loading_range   = c(0.60, 0.78),
#'   factor_cor      = c(A_B = 0.35, A_C = 0.20, B_C = 0.45),
#'   threshold_shift = c(A = -0.15, B = 0.00, C = 0.10),
#'   concentration   = 1.5,
#'   likert_points   = 7,
#'   n_resid_cov     = 3,
#'   factor_items    = list(
#'     A = paste0("a", 1:5),
#'     B = paste0("b", 1:4),
#'     C = paste0("c", 1:4)
#'   ),
#'   resid_cov_range = c(0.06, 0.14)
#' )
#' str(dat2)
#'
#' # --- Example 3: Matrix factor_cor (four factors) ---
#' set.seed(999)
#' cor_mat <- matrix(c(
#'   1.00, 0.35, 0.20, 0.15,
#'   0.35, 1.00, 0.40, 0.25,
#'   0.20, 0.40, 1.00, 0.50,
#'   0.15, 0.25, 0.50, 1.00
#' ), nrow = 4, byrow = TRUE)
#'
#' dat3 <- simulate_cfa(
#'   n_obs           = 500,
#'   loading_range   = c(0.63, 0.79),
#'   factor_cor      = cor_mat,
#'   threshold_shift = c(F1 = -0.20, F2 = 0.00, F3 = 0.10, F4 = -0.10),
#'   concentration   = 1.5,
#'   likert_points   = 5,
#'   n_resid_cov     = 4,
#'   factor_items    = list(
#'     F1 = paste0("f1_", 1:5),
#'     F2 = paste0("f2_", 1:4),
#'     F3 = paste0("f3_", 1:4),
#'     F4 = paste0("f4_", 1:3)
#'   ),
#'   resid_cov_range = c(0.05, 0.12)
#' )
#' str(dat3)
#'
#' @seealso \code{\link[lavaan]{simulateData}}
#' @export
simulate_cfa <- function(
    n_obs = 300,
    loading_range = c(0.68, 0.80),
    factor_cor = 0.25,
    threshold_shift = c(a = -0.10, b = -0.10, c = -0.10, d = -0.10, e = -0.10),
    concentration = 1.5,
    likert_points = 5,
    n_resid_cov = 5,
    factor_items = list(
      a = paste0("a", 1:8),
      b = paste0("b", 1:4),
      c = paste0("c", 1:5),
      d = paste0("d", 1:6),
      e = paste0("e", 1:6)
    ),
    resid_cov_range = c(0.08, 0.15)
) {

  # --- Validate arguments ---
  if (!is.numeric(likert_points) || length(likert_points) != 1 ||
      likert_points < 2 || likert_points != as.integer(likert_points)) {
    stop("`likert_points` must be an integer >= 2.")
  }
  if (length(loading_range) != 2 || loading_range[1] >= loading_range[2]) {
    stop("`loading_range` must be length-2 with lower < upper.")
  }
  if (length(resid_cov_range) != 2 || resid_cov_range[1] >= resid_cov_range[2]) {
    stop("`resid_cov_range` must be length-2 with lower < upper.")
  }
  if (!is.numeric(concentration) || length(concentration) != 1 || concentration <= 0) {
    stop("`concentration` must be a single positive number.")
  }

  factor_names <- names(factor_items)
  if (is.null(factor_names) || any(factor_names == "")) {
    stop("`factor_items` must be a named list, e.g., list(a = ..., b = ...).")
  }
  n_factors <- length(factor_names)

  # --- Parse factor_cor ---
  if (is.matrix(factor_cor)) {
    if (nrow(factor_cor) != n_factors || ncol(factor_cor) != n_factors) {
      stop("`factor_cor` matrix dimensions must match the number of factors.")
    }
    cor_values <- list()
    for (i in seq_len(n_factors - 1)) {
      for (j in (i + 1):n_factors) {
        pair_key <- paste(factor_names[i], factor_names[j], sep = "_")
        cor_values[[pair_key]] <- factor_cor[i, j]
      }
    }
  } else if (is.numeric(factor_cor) && is.null(dim(factor_cor))) {
    if (length(factor_cor) == 1) {
      if (factor_cor <= -1 || factor_cor >= 1) {
        stop("Scalar `factor_cor` must be in (-1, 1).")
      }
      cor_values <- list()
      for (i in seq_len(n_factors - 1)) {
        for (j in (i + 1):n_factors) {
          pair_key <- paste(factor_names[i], factor_names[j], sep = "_")
          cor_values[[pair_key]] <- factor_cor
        }
      }
    } else if (!is.null(names(factor_cor))) {
      cor_values <- as.list(factor_cor)
    } else {
      expected_n <- n_factors * (n_factors - 1) / 2
      if (length(factor_cor) != expected_n) {
        stop(
          sprintf(
            "Unnamed `factor_cor` vector must have length 1 or %d (number of factor pairs).",
            expected_n
          )
        )
      }
      cor_values <- list()
      idx <- 1
      for (i in seq_len(n_factors - 1)) {
        for (j in (i + 1):n_factors) {
          pair_key <- paste(factor_names[i], factor_names[j], sep = "_")
          cor_values[[pair_key]] <- factor_cor[idx]
          idx <- idx + 1
        }
      }
    }
  } else {
    stop("`factor_cor` must be a numeric scalar, named vector, unnamed vector, or matrix.")
  }

  if (any(unlist(cor_values) <= -1) || any(unlist(cor_values) >= 1)) {
    stop("All factor correlations must be in (-1, 1).")
  }

  # --- Align threshold_shift ---
  if (is.null(names(threshold_shift))) {
    if (length(threshold_shift) == 1) {
      threshold_shift <- rep(threshold_shift, n_factors)
    }
    names(threshold_shift) <- factor_names
  }
  if (!all(factor_names %in% names(threshold_shift))) {
    stop("`threshold_shift` must include all factors in `factor_items`.")
  }
  threshold_shift <- threshold_shift[factor_names]

  # --- Loading syntax ---
  lambda_list <- lapply(factor_items, function(items) {
    runif(length(items), loading_range[1], loading_range[2])
  })

  loading_lines <- vapply(
    seq_along(factor_names),
    function(k) {
      f <- factor_names[k]
      items <- factor_items[[k]]
      lam <- lambda_list[[k]]
      paste0(f, " =~ ", paste(sprintf("%.3f*%s", lam, items), collapse = " + "))
    },
    character(1)
  )

  # --- Factor correlation syntax ---
  cor_lines <- character(0)
  if (n_factors >= 2) {
    for (i in seq_len(n_factors - 1)) {
      for (j in (i + 1):n_factors) {
        pair_key <- paste(factor_names[i], factor_names[j], sep = "_")
        val <- cor_values[[pair_key]]
        cor_lines <- c(
          cor_lines,
          sprintf("%s ~~ %.3f*%s", factor_names[i], val, factor_names[j])
        )
      }
    }
  }

  # --- Threshold syntax ---
  base_thresholds <- qnorm(seq_len(likert_points - 1) / likert_points)

  make_threshold_line <- function(item, shift) {
    thresholds <- base_thresholds * concentration + shift
    terms <- paste(
      sprintf("%.2f*t%d", thresholds, seq_along(thresholds)),
      collapse = " + "
    )
    sprintf("%s | %s", item, terms)
  }

  threshold_lines <- unlist(
    lapply(factor_names, function(f) {
      vapply(
        factor_items[[f]], make_threshold_line, character(1),
        shift = threshold_shift[[f]]
      )
    }),
    use.names = FALSE
  )

  # --- Within-factor residual covariance syntax ---
  all_items <- unlist(factor_items, use.names = FALSE)

  within_pairs <- do.call(rbind, lapply(factor_items, function(items) {
    if (length(items) < 2) return(NULL)
    t(combn(items, 2))
  }))

  n_resid_cov <- min(max(0L, as.integer(n_resid_cov)), nrow(within_pairs))

  resid_lines <- character(0)
  if (n_resid_cov > 0 && nrow(within_pairs) > 0) {
    sel <- sample(seq_len(nrow(within_pairs)), size = n_resid_cov, replace = FALSE)
    sel_pairs <- within_pairs[sel, , drop = FALSE]
    resid_vals <- runif(n_resid_cov, resid_cov_range[1], resid_cov_range[2])
    resid_lines <- vapply(
      seq_len(n_resid_cov),
      function(i) sprintf("%s ~~ %.3f*%s", sel_pairs[i, 1], resid_vals[i], sel_pairs[i, 2]),
      character(1)
    )
  }

  # --- Assemble model syntax ---
  model_syn <- paste(
    paste(loading_lines, collapse = "\n"),
    paste(cor_lines, collapse = "\n"),
    paste(threshold_lines, collapse = "\n"),
    paste(resid_lines, collapse = "\n"),
    sep = "\n\n"
  )

  # --- Simulate and return ---
  data_sim <- lavaan::simulateData(model_syn, sample.nobs = n_obs) |>
    tibble::as_tibble() |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(all_items),
      ~ as.numeric(as.character(.x))
    ))

  data_sim
}

