#' Convert Z-scores to a Likert scale (with optional skew)
#'
#' Discretizes continuous Z-scores into \eqn{1..n\_levels} Likert categories so that
#' the **marginal category proportions** follow either:
#' (a) the standard normal density (when `skew = 0`), or
#' (b) a univariate skew-normal density (when `skew != 0`).
#'
#' Internally, this function uses a helper \code{discretize_density()} to compute
#' \emph{endpoints} that partition the real line into \code{n_levels} bins whose
#' areas (probabilities) match the target density. It then cuts the input \code{z_score}
#' at those endpoints and relabels the bins as integer Likert levels.
#'
#' @param z_score Numeric vector of (approximately) standard normal values to be discretized.
#'   Missing values are preserved.
#' @param n_levels Integer (>= 2). Number of Likert categories to create. Defaults to 5.
#' @param skew Numeric, controls asymmetry of the **target** density used to set the cuts.
#'   Use \code{0} for symmetric normal (default). Nonzero values induce skewness via the
#'   skew-normal family. This argument is transformed to a skew-normal \code{alpha}
#'   (shape) by \eqn{ \alpha = \frac{\text{skew}}{\sqrt{1-\text{skew}^2}} \sqrt{2/\pi} }.
#'   Reasonable range is \code{-0.95} to \code{0.95}. Extreme values may be unstable.
#'
#' @return An integer vector with values in \code{1:n_levels}. If \code{as.ordered} behavior
#'   is desired, convert the result to \code{ordered()} outside this function.
#'
#' @details
#' - This function **does not** re-estimate any density from the data. It only uses the
#'   chosen reference density (normal or skew-normal) to define **proportionally balanced
#'   cutpoints** and then bins the provided \code{z_score}.
#' - The helper \code{discretize_density(density_fn, n_levels, eps)} must exist in your
#'   environment (usually elsewhere in the same package). It is expected to return a list
#'   with a numeric vector \code{endp} of length \code{n_levels + 1}, i.e., the cumulative
#'   probability endpoints (on the Z-scale) that define the category boundaries.
#'
#' @examples
#' # Basic: 5-point Likert with symmetric category proportions
#' set.seed(1)
#' z <- rnorm(1000)
#' x5 <- likertize_zscore(z, n_levels = 5)
#' table(x5) / length(x5)  # ~balanced under normal reference
#'
#' # Skewed category proportions (right-skew in the reference density)
#' x5_skew <- likertize_zscore(z, n_levels = 5, skew = 0.5)
#' table(x5_skew) / length(x5_skew)
#'
#' # 7-point scale
#' x7 <- likertize_zscore(z, n_levels = 7)
#' head(x7)
#'
#' @seealso discretize_density
#' @export
#' @importFrom stats dnorm
#' @importFrom sn dsn
likertize_zscore <- function(z_score, n_levels = 5, skew = 0) {

  # ---- Basic input validation ------------------------------------------------
  if (!is.numeric(z_score)) {
    stop("`z_score` must be a numeric vector.")
  }
  if (length(n_levels) != 1L || !is.finite(n_levels) || n_levels < 2) {
    stop("`n_levels` must be a single finite integer >= 2.")
  }
  if (length(skew) != 1L || !is.finite(skew)) {
    stop("`skew` must be a single finite numeric value.")
  }
  # Mild safeguard against pathological alpha transformation (division by ~0)
  if (abs(skew) >= 0.999) {
    stop("`skew` is too extreme; use a value in (-0.999, 0.999).")
  }

  # ---- Choose the reference density used to compute balanced cutpoints -------
  if (skew == 0) {
    # Symmetric standard normal density
    result <- discretize_density(
      density_fn = stats::dnorm,  # <- target density for proportional binning
      n_levels  = n_levels,
      eps       = 1e-06
    )
  } else {
    # Skew-normal density: convert user-friendly 'skew' to shape 'alpha'
    # NOTE: 'skew' here is mapped to the skew-normal shape parameter alpha.
    #       The transformation below mirrors the author's original intent.
    alpha <- skew / sqrt(1 - skew^2) * sqrt(2 / pi)

    # Density function handle for sn::dsn with mean 0 and scale 1
    skew_density <- function(x) sn::dsn(x, xi = 0, omega = 1, alpha = alpha)

    result <- discretize_density(
      density_fn = skew_density,
      n_levels   = n_levels,
      eps        = 1e-06
    )
  }

  # ---- Translate endpoints into cut breaks and bin the data ------------------
  endpoints <- result$endp
  # Expected: length(endpoints) == n_levels + 1
  if (!is.numeric(endpoints) || length(endpoints) < 2) {
    stop("`discretize_density()` returned invalid `endp`. Check its implementation.")
  }

  # Build open-ended breaks so extreme values are captured
  breaks <- c(-Inf, endpoints[2:n_levels], Inf)

  # Cut into 1..n_levels and preserve NA from input
  likert <- cut(
    z_score,
    breaks         = breaks,
    labels         = 1:n_levels,
    include.lowest = TRUE
  )

  # Return as plain integers (1..n_levels)
  as.integer(as.character(likert))
}
