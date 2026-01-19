#' Convert Z-scores to a Likert scale (skew-normal reference)
#'
#' Discretizes continuous Z-scores into \eqn{1..n\_levels} Likert categories so that
#' the **marginal category proportions** follow a chosen reference density:
#' \itemize{
#'   \item standard normal when \code{skew = 0}
#'   \item univariate skew-normal when \code{skew != 0}, where \code{skew} is the skew-normal
#'         shape parameter \code{alpha}
#' }
#'
#' Internally, this function calls \code{discretize_density()} to compute endpoints that
#' partition the real line into \code{n_levels} bins whose areas match the target density,
#' then bins \code{z_score} using those endpoints.
#'
#' @param z_score Numeric vector of (approximately) standard normal values to be discretized.
#'   Missing values are preserved.
#' @param n_levels Integer (>= 2). Number of Likert categories to create. Defaults to 5.
#' @param skew Numeric. The skew-normal shape parameter \code{alpha} used to define the
#'   reference density when computing cutpoints. Use \code{0} for a symmetric normal
#'   reference. Positive values shift probability mass toward higher categories; negative
#'   values shift probability mass toward lower categories.
#'
#'   Practical guidance: \code{abs(skew) <= 3} often yields mild to moderate skew in the
#'   induced category proportions. Larger magnitudes may cause highly imbalanced categories,
#'   including empty or near-empty levels, which can reduce information for downstream analyses.
#'
#' @return An integer vector with values in \code{1:n_levels}. If ordered-factor output is
#'   desired, convert the result using \code{ordered()} outside this function.
#'
#' @details
#' - This function does not estimate any density from the data. It only uses the chosen
#'   reference density (normal or skew-normal) to define proportional cutpoints, then bins
#'   the provided \code{z_score}.
#' - The helper \code{discretize_density(density_fn, n_levels, eps)} must exist in your
#'   environment. It is expected to return a list with numeric vector \code{endp} of length
#'   \code{n_levels + 1}, the endpoints on the Z-scale that define the category boundaries.
#'
#' @examples
#' set.seed(1)
#' z <- rnorm(1000)
#'
#' # 5-point, symmetric reference
#' x5 <- likertize_zscore(z, n_levels = 5, skew = 0)
#' table(x5) / length(x5)
#'
#' # 5-point, skew-normal reference (shape = 2)
#' x5_s2 <- likertize_zscore(z, n_levels = 5, skew = 2)
#' table(x5_s2) / length(x5_s2)
#'
#' # stronger skew (shape = 5), may yield imbalanced categories
#' x5_s5 <- likertize_zscore(z, n_levels = 5, skew = 5)
#' table(x5_s5) / length(x5_s5)
#'
#' @seealso discretize_density
#' @export
#' @importFrom stats dnorm
#' @importFrom sn dsn
likertize_zscore <- function(z_score, n_levels = 5, skew = 0) {

  if (!is.numeric(z_score)) {
    stop("`z_score` must be a numeric vector.")
  }
  if (length(n_levels) != 1L || !is.finite(n_levels) || n_levels < 2) {
    stop("`n_levels` must be a single finite integer >= 2.")
  }
  if (length(skew) != 1L || !is.finite(skew)) {
    stop("`skew` must be a single finite numeric value.")
  }

  # Informational warnings for extreme shape values
  if (abs(skew) > 5) {
    warning("Large `skew` may yield highly imbalanced or empty categories; check the resulting frequency table.")
  } else if (abs(skew) > 3) {
    warning("Moderately large `skew` may yield imbalanced categories; consider abs(skew) <= 3 for typical use.")
  }

  if (isTRUE(all.equal(skew, 0))) {
    result <- discretize_density(
      density_fn = stats::dnorm,
      n_levels   = n_levels,
      eps        = 1e-06
    )
  } else {
    skew_density <- function(x) sn::dsn(x, xi = 0, omega = 1, alpha = skew)

    result <- discretize_density(
      density_fn = skew_density,
      n_levels   = n_levels,
      eps        = 1e-06
    )
  }

  endpoints <- result$endp
  if (!is.numeric(endpoints) || length(endpoints) < 2) {
    stop("`discretize_density()` returned invalid `endp`. Check its implementation.")
  }

  breaks <- c(-Inf, endpoints[2:n_levels], Inf)

  likert <- cut(
    z_score,
    breaks         = breaks,
    labels         = 1:n_levels,
    include.lowest = TRUE
  )

  as.integer(as.character(likert))
}
