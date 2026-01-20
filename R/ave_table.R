#' Create an AVE and CR table with factor correlations (convergent and discriminant validity)
#'
#' Fits a CFA model using \pkg{lavaan} and returns a validity table containing
#' Composite Reliability (CR), Average Variance Extracted (AVE), and the latent
#' factor correlation matrix with \eqn{\sqrt{AVE}} on the diagonal.
#'
#' The returned table is commonly used to evaluate:
#' \itemize{
#'   \item Convergent validity via AVE (often AVE >= .50) and CR (often CR >= .70).
#'   \item Discriminant validity via the Fornell Larcker criterion, where
#'     \eqn{\sqrt{AVE}} for each factor is expected to exceed its correlations with other factors.
#' }
#'
#' All numeric entries are formatted in APA style by removing the leading zero
#' when \eqn{|x| < 1} (for example, \code{.810}, \code{.436***}, \code{-.215}).
#'
#' @param data A \code{data.frame} (or tibble) containing observed indicators used in \code{model}.
#' @param model A lavaan model string specifying the measurement model (for example, \code{F1 =~ x1 + x2 + x3}).
#' @param digits_stat Integer, number of digits for CR, AVE, and \eqn{\sqrt{AVE}} on the diagonal.
#' @param digits_corr Integer, number of digits for factor correlations.
#' @param show_stars Logical, whether to append significance stars to correlations based on p values.
#' @param ... Additional arguments passed to \code{lavaan::cfa()}, such as \code{estimator},
#'   \code{missing}, \code{std.lv}, \code{ordered}, and \code{parameterization}.
#'
#' @return A \code{data.frame} with columns \code{Domain}, \code{CR}, \code{AVE}, and numeric
#'   index columns \code{1..k} for the correlation matrix, where \code{k} is the number of factors.
#'
#' @details
#' CR and AVE are computed from the standardized solution (\code{std.all}).
#' For a given factor with standardized loadings \eqn{\lambda_i} and indicator residual variances \eqn{\theta_i}:
#' \deqn{CR = \frac{(\sum \lambda_i)^2}{(\sum \lambda_i)^2 + \sum \theta_i}}
#' \deqn{AVE = \frac{\sum \lambda_i^2}{\sum \lambda_i^2 + \sum \theta_i}}
#'
#' If an item residual variance is not available in the standardized solution, it will be treated as \code{NA}
#' and CR and AVE for that factor will become \code{NA}.
#'
#' @examples
#' if (requireNamespace("lavaan", quietly = TRUE)) {
#'   set.seed(123)
#'   model_pop <- '
#'     F1 =~ 0.8*x1 + 0.7*x2 + 0.6*x3
#'     F2 =~ 0.8*y1 + 0.7*y2 + 0.6*y3
#'     F1 ~~ 0.4*F2
#'   '
#'   dat <- lavaan::simulateData(model_pop, sample.nobs = 300)
#'
#'   model_fit <- '
#'     F1 =~ x1 + x2 + x3
#'     F2 =~ y1 + y2 + y3
#'   '
#'
#'   ave_table(dat, model_fit, digits_stat = 3, digits_corr = 3)
#' }
#'
#' @export
ave_table <- function(data,
                      model,
                      digits_stat = 3,
                      digits_corr = 3,
                      show_stars  = TRUE,
                      ...) {

  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("Install lavaan to use ave_table().")
  }

  fit <- lavaan::cfa(model = model, data = data, ...)

  pe <- lavaan::parameterEstimates(fit, standardized = TRUE)

  lv_names <- lavaan::lavNames(fit, type = "lv")
  if (length(lv_names) == 0) stop("No latent variables found in the fitted model.")

  k <- length(lv_names)

  # APA numeric formatter: remove leading zero when |x| < 1
  fmt_apa <- function(x, digits) {
    if (is.na(x)) return("")
    s <- sprintf(paste0("%.", digits, "f"), x)
    sub("^(-?)0\\.", "\\1.", s)
  }

  get_stars <- function(p_val) {
    if (!show_stars) return("")
    if (is.na(p_val))  return("")
    if (p_val < 0.001) return("***")
    if (p_val < 0.01)  return("**")
    if (p_val < 0.05)  return("*")
    ""
  }

  load_pe <- pe[pe$op == "=~" & pe$lhs %in% lv_names, , drop = FALSE]
  res_pe  <- pe[pe$op == "~~" & pe$lhs == pe$rhs, , drop = FALSE]

  CR  <- rep(NA_real_, k)
  AVE <- rep(NA_real_, k)

  for (i in seq_len(k)) {
    f <- lv_names[i]
    rows <- load_pe[load_pe$lhs == f, , drop = FALSE]
    if (nrow(rows) == 0) next

    inds   <- rows$rhs
    lambda <- rows$std.all

    theta_rows <- res_pe[res_pe$lhs %in% inds, , drop = FALSE]
    theta <- theta_rows$std.all
    names(theta) <- theta_rows$lhs
    theta <- theta[inds]

    if (any(is.na(lambda)) || any(is.na(theta))) {
      CR[i]  <- NA_real_
      AVE[i] <- NA_real_
      next
    }

    sum_lambda <- sum(lambda)
    sum_theta  <- sum(theta)

    CR[i]  <- (sum_lambda^2) / ((sum_lambda^2) + sum_theta)
    AVE[i] <- sum(lambda^2) / (sum(lambda^2) + sum_theta)
  }

  sqrt_ave <- sqrt(AVE)

  mat <- matrix("", nrow = k, ncol = k)

  for (i in seq_len(k)) {
    mat[i, i] <- fmt_apa(sqrt_ave[i], digits_stat)
  }

  corr_pe <- pe[
    pe$op == "~~" &
      pe$lhs %in% lv_names &
      pe$rhs %in% lv_names &
      pe$lhs != pe$rhs,
    ,
    drop = FALSE
  ]

  for (i in seq_len(k)) {
    for (j in seq_len(k)) {
      if (i > j) {
        fi <- lv_names[i]
        fj <- lv_names[j]

        idx <- which((corr_pe$lhs == fi & corr_pe$rhs == fj) | (corr_pe$lhs == fj & corr_pe$rhs == fi))
        if (length(idx) == 0) next

        r_val <- corr_pe$std.all[idx[1]]
        p_val <- corr_pe$pvalue[idx[1]]

        if (is.na(r_val)) next

        r_str <- fmt_apa(r_val, digits_corr)
        s_str <- get_stars(p_val)
        mat[i, j] <- paste0(r_str, s_str)
      }
    }
  }

  out <- data.frame(
    Domain = paste0(seq_len(k), ". ", lv_names),
    CR     = vapply(CR,  fmt_apa, character(1), digits = digits_stat),
    AVE    = vapply(AVE, fmt_apa, character(1), digits = digits_stat),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  mat_df <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(mat_df) <- as.character(seq_len(k))

  cbind(out, mat_df)
}
