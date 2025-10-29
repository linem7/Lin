#' Corrected test_args for crosstable
#'
#' A customized \pkg{crosstable} \code{test_args} function designed to apply
#' *Welch’s unequal variance tests* and *Yates’ correction* when creating
#' descriptive summary tables.
#'
#' Specifically:
#' \itemize{
#'   \item For two-group numeric variables: uses the \strong{Welch t-test} (\code{var.equal = FALSE});
#'   \item For multi-group numeric variables: uses \strong{Welch ANOVA} (\code{oneway.test(..., var.equal = FALSE)});
#'   \item For 2×2 categorical tables: applies \strong{Pearson’s Chi-squared test with Yates’ correction};
#'   \item For larger categorical tables: applies \strong{standard Pearson’s Chi-squared test (no correction)}.
#' }
#'
#' The output printed in the \code{test} column includes:
#' \enumerate{
#'   \item the test statistic (t / F / χ²);
#'   \item degrees of freedom (including fractional df for Welch tests);
#'   \item p-value;
#'   \item the method name.
#' }
#'
#' @return A list object that can be passed to \code{crosstable()} via the
#' \code{test_args} argument.
#' @seealso \code{\link[crosstable]{crosstable_test_args}}, \code{\link[crosstable]{crosstable}}
#' @export
#' @importFrom crosstable crosstable_test_args plim
#' @importFrom stats t.test oneway.test chisq.test
#'
#' @examples
#' \dontrun{
#' library(crosstable)
#' library(flextable)
#'
#' # Example dataset
#' set.seed(123)
#' dat <- data.frame(
#'   group = factor(rep(c("A", "B"), each = 30)),
#'   score = c(rnorm(30, 5, 1), rnorm(30, 6, 1.5)),
#'   gender = factor(sample(c("M", "F"), 60, replace = TRUE))
#' )
#'
#' # Use the custom corrected test_args
#' tab <- crosstable(
#'   dat,
#'   cols = c(score, gender),
#'   by = group,
#'   test = TRUE,
#'   test_args = Lin::crosstable_args()
#' )
#'
#' flextable::as_flextable(tab)
#' }
crosstable_args <- function(){
  crosstable::crosstable_test_args(
    # ---- Numeric variable vs grouping ----
    test_summarize = function(x, g){
      x <- as.numeric(x)
      ng <- length(table(g))
      if (ng <= 1) return(list(p.value=NULL, method=NULL))

      if (ng == 2) {
        tt <- stats::t.test(x ~ g, var.equal = FALSE)  # Welch t
        list(
          p.value   = tt$p.value,
          method    = "Welch t-test",
          statistic = unname(tt$statistic),
          parameter = unname(tt$parameter),   # df (may be fractional)
          stat_name = "t"
        )
      } else {
        w1 <- stats::oneway.test(x ~ g, var.equal = FALSE)  # Welch ANOVA
        list(
          p.value   = w1$p.value,
          method    = "Welch ANOVA",
          statistic = unname(w1$statistic),
          parameter = unname(w1$parameter),   # c(df1, df2)
          stat_name = "F"
        )
      }
    },

    # ---- Categorical vs categorical ----
    test_tabular = function(x, y){
      tab <- table(x, y)
      if (all(dim(tab) == c(2, 2))) {
        cs <- suppressWarnings(stats::chisq.test(tab, correct = TRUE))
        mth <- "Pearson's Chi-squared with Yates correction"
      } else {
        cs <- suppressWarnings(stats::chisq.test(tab, correct = FALSE))
        mth <- "Pearson's Chi-squared test"
      }
      list(
        p.value   = cs$p.value,
        method    = mth,
        statistic = unname(cs$statistic),
        parameter = unname(cs$parameter),   # df
        stat_name = "\u03C7\u00B2"          # χ²
      )
    },

    # ---- Display format ----
    test_display = function(test, digits = 4, method = TRUE){
      if (all(sapply(test, is.null))) return("No test")
      ptxt <- crosstable::plim(test$p.value, digits = digits)

      stat_lab <- if (!is.null(test$stat_name)) test$stat_name else "stat"
      stat_txt <- if (!is.null(test$statistic))
        paste0(stat_lab, " = ", sprintf("%.2f", test$statistic)) else ""

      fmt_df <- function(v){
        if (is.null(v)) return("")
        if (length(v) == 1) {
          if (abs(v - round(v)) < 1e-6) paste0("df = ", round(v))
          else paste0("df = ", sprintf("%.1f", v))
        } else {
          p1 <- v[1]; p2 <- v[2]
          d1 <- if (abs(p1 - round(p1)) < 1e-6) as.character(round(p1)) else sprintf("%.1f", p1)
          d2 <- if (abs(p2 - round(p2)) < 1e-6) as.character(round(p2)) else sprintf("%.1f", p2)
          paste0("df1 = ", d1, ", df2 = ", d2)
        }
      }
      df_txt <- fmt_df(test$parameter)
      meth_txt <- if (method && !is.null(test$method)) paste0(" (", test$method, ")") else ""

      paste0(stat_txt,
             if (nzchar(stat_txt) && nzchar(df_txt)) ", " else "",
             df_txt,
             if (nzchar(stat_txt) || nzchar(df_txt)) ", " else "",
             "p = ", ptxt,
             meth_txt)
    },

    show_method = TRUE
  )
}
