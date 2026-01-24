#' Convert Z-scores to a Likert scale (Metric/Equidistant approach)
#'
#' Discretizes continuous Z-scores into \eqn{1..n\_levels} Likert categories by defining
#' **equidistant intervals** over the effective range of a reference density.
#'
#' Unlike quantile-based discretization (which forces a uniform distribution), this
#' "metric" approach partitions the underlying Z-scale into equal-width bins. Consequently:
#' \itemize{
#'   \item If the reference is standard normal (\code{skew = 0}), the resulting Likert
#'     distribution will be bell-shaped (unimodal, symmetric).
#'   \item If the reference is skew-normal (\code{skew != 0}), the intervals shift,
#'     creating a skewed Likert distribution that mimics the reference shape.
#' }
#'
#' @param z_score Numeric vector of (approximately) standard normal values to be discretized.
#'   Missing values are preserved.
#' @param n_levels Integer (>= 2). Number of Likert categories to create. Defaults to 5.
#' @param skew Numeric. The skew-normal shape parameter \code{alpha} used to define the
#'   reference effective range. Use \code{0} for a symmetric normal reference.
#'   \itemize{
#'     \item \strong{Positive skew (\code{skew > 0}):} shifts the range limits rightward.
#'       Given standard normal data, this results in mass concentrating in **lower** categories.
#'     \item \strong{Negative skew (\code{skew < 0}):} shifts the range limits leftward.
#'       Given standard normal data, this results in mass concentrating in **higher** categories.
#'   }
#'   Practical guidance: \code{abs(skew) <= 3} often yields mild to moderate skew.
#'
#' @return An integer vector with values in \code{1:n_levels}.
#'
#' @details
#' \strong{Strategy: Cut the Metric, not the Density}
#'
#' Instead of finding cutpoints that integrate to specific probabilities (which yields
#' uniform counts), this function:
#' \enumerate{
#'   \item Determines the "effective range" of the reference distribution (covering
#'     the 0.1\% to 99.9\% quantiles).
#'   \item Divides this numeric range into \code{n_levels} bins of equal width.
#'   \item Maps the provided \code{z_score} into these bins.
#' }
#' This preserves the relative distances between data points, ensuring that central
#' values remain in the central categories (bell-shaped) rather than being spread
#' out to fill the tails.
#'
#' @examples
#' set.seed(1)
#' z <- rnorm(1000)
#'
#' # 1. Symmetric reference (skew = 0)
#' # Result: Bell-shaped distribution (few 1s/5s, many 3s)
#' x5 <- likertize_zscore(z, n_levels = 5, skew = 0)
#' barplot(table(x5), main = "Skew = 0 (Bell-shaped)")
#'
#' # 2. Positive skew reference (skew = 2)
#' # Result: Right-skewed (Long tail right -> Mass on left/lower categories)
#' x5_s2 <- likertize_zscore(z, n_levels = 5, skew = 2)
#' barplot(table(x5_s2), main = "Skew = 2 (Right-skewed)")
#'
#' @export
#' @importFrom stats qnorm
#' @importFrom sn qsn
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

  # 1. Define the coverage probability to determine the "Effective Range"
  # We use 0.001 to 0.999 to cover 99.8% of the theoretical mass,
  # treating values outside this as the "infinite" tails.
  probs <- c(0.001, 0.999)

  # 2. Calculate range limits based on reference density
  if (isTRUE(all.equal(skew, 0))) {
    # Reference: Standard Normal
    range_lims <- stats::qnorm(probs)
  } else {
    # Reference: Skew-Normal
    # Note: xi=0, omega=1 keeps it on the standard scale, only alpha varies
    range_lims <- sn::qsn(probs, xi = 0, omega = 1, alpha = skew)
  }

  # 3. Compute Equidistant Cutpoints ("Cut the Metric")
  # We divide the effective range into n_levels equal-width intervals.
  # seq() includes both start and end, so length.out is n_levels + 1
  endpoints <- seq(from = range_lims[1], to = range_lims[2], length.out = n_levels + 1)

  # 4. Construct Breaks for cut()
  # We replace the outer boundaries with Inf/-Inf to ensure no data is lost
  # if it falls slightly outside the theoretical 99.8% range.
  breaks <- endpoints
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf

  # 5. Bin the data
  likert <- cut(
    z_score,
    breaks = breaks,
    labels = 1:n_levels,
    include.lowest = TRUE
  )

  as.integer(as.character(likert))
}
