#' Item Analysis for Classical Test Theory
#'
#' Computes a variety of item-level statistics for scales (e.g., means, SDs,
#' skewness, kurtosis, extreme-group discrimination, Cronbach's alpha if item deleted,
#' and corrected item-total correlations with significance stars), plus overall
#' reliability.
#'
#' @param data A data frame containing item response columns.
#' @param ... <[`tidy-select`][dplyr::dplyr_tidy_select]> Column selection
#'   supporting all tidyselect syntax. You can use:
#'   \itemize{
#'     \item Column ranges: \code{item1:item5}
#'     \item Specific columns: \code{c(item1, item3, item5)}
#'     \item Selection helpers: \code{starts_with("A")}, \code{contains("scale")}
#'     \item Custom pattern helper: \code{pattern("Q{i}_rev", 1:5)}
#'   }
#'   These forms can be combined freely within a single call.
#' @param digits Integer; number of decimal places for skewness, kurtosis, CR, CITC,
#'   and alpha-if-deleted values. Means and SDs are always formatted to two decimals.
#'   Default is 3.
#' @param total Logical; if \code{TRUE}, appends a "Total" column with the overall
#'   Cronbach's alpha in the first row and blanks elsewhere. Default is \code{TRUE}.
#'
#' @return A data frame with one row per item and the following columns:
#' \item{Item}{Item name.}
#' \item{Mean}{Item mean score (2 decimal places).}
#' \item{SD}{Standard deviation of the item (2 decimal places).}
#' \item{Skew}{Skewness of the item distribution.}
#' \item{Kurt}{Kurtosis of the item distribution.}
#' \item{CR}{t-value from high vs. low group discrimination, with significance stars.}
#' \item{Alpha.if" deleted}{Cronbach's alpha when the item is removed.}
#' \item{CITC}{Corrected item-total correlation, with significance stars.}
#' \item{Total}{(if \code{total = TRUE}) Overall alpha in the first row; blanks in remaining rows.}
#'
#' @importFrom dplyr select mutate left_join across filter pull everything case_when
#' @importFrom tidyselect eval_select
#' @importFrom rlang expr
#' @importFrom tibble rownames_to_column
#' @importFrom psych describe
#' @importFrom stats t.test cor.test quantile
#' @importFrom car leveneTest
#'
#' @details
#' Column selection is powered by \code{tidyselect}, providing the same flexibility
#' available in \code{dplyr::select()}. The custom helper \code{pattern()} allows
#' glue-style name generation (e.g., \code{pattern("item{i}", 1:5)} resolves to
#' \code{item1, item2, ..., item5}) and can be used alongside any other tidyselect
#' expression.
#'
#' To assess internal consistency, the function computes overall reliability using
#' \code{psych::alpha()} (with \code{check.keys = TRUE}) and extracts the
#' \code{alpha.drop} table to obtain Cronbach's alpha when each item is removed.
#'
#' Descriptive statistics, including mean, standard deviation, skewness, and kurtosis,
#' are calculated with \code{psych::describe()}. Means and SDs are formatted to two
#' decimal places, while skewness and kurtosis are rounded according to the
#' \code{digits} argument.
#'
#' Extreme-group discrimination (CR) is assessed by splitting respondents into high
#' (top 27\%) and low (bottom 27\%) groups based on total scores computed as row
#' means across selected items. Levene's test (\code{car::leveneTest()}) determines
#' whether equal-variance or Welch's t-test is applied. The resulting t-values are
#' annotated with significance stars (* p < .05, ** p < .01, *** p < .001).
#'
#' The corrected item-total correlation (CITC) is computed by correlating each item
#' with the sum of all remaining items, tested via \code{stats::cor.test()}, and
#' formatted with the same significance star convention.
#'
#' All results are merged into a single data frame. If \code{total = TRUE}, the
#' overall Cronbach's alpha is displayed in a "Total" column in the first row.
#'
#' @examples
#' data(good_rel)
#' data(poor_rel)
#'
#' # --- Column range ---
#' item_analysis(good_rel, item1:item5)
#'
#' # --- Specific columns ---
#' item_analysis(good_rel, c(item1, item3, item5))
#'
#' # --- Selection helper ---
#' item_analysis(good_rel, starts_with("item"))
#'
#' # --- Pattern helper (glue-style) ---
#' item_analysis(good_rel, pattern("item{i}", 1:5))
#'
#' # --- Combine multiple selectors ---
#' # item_analysis(my_data, starts_with("qp"), starts_with("hc"))
#'
#' # --- Format with apa() ---
#' item_analysis(good_rel, item1:item5) %>% apa()
#'
#' # --- Poor-reliability example ---
#' item_analysis(poor_rel, pattern("item{i}", 1:5))
#'
#' @export

item_analysis <- function(data,
                          ...,
                          digits = 3,
                          total  = TRUE) {

  #── Dependencies
  if (!requireNamespace("dplyr", quietly = TRUE))      stop("Install dplyr")
  if (!requireNamespace("psych", quietly = TRUE))       stop("Install psych")
  if (!requireNamespace("tidyselect", quietly = TRUE))  stop("Install tidyselect")
  if (!requireNamespace("tibble", quietly = TRUE))      stop("Install tibble")

  #── 1. Resolve column selection via tidyselect
  col_pos <- tidyselect::eval_select(rlang::expr(c(...)), data = data)

  if (length(col_pos) == 0) {
    stop("No columns matched. Please check your selection.")
  }

  item_names     <- names(col_pos)
  selected_items <- data[, col_pos, drop = FALSE]

  #── 2. Compute total alpha and get alpha.drop
  alpha_res <- suppressWarnings(psych::alpha(selected_items, check.keys = TRUE))
  ad <- alpha_res$alpha.drop
  if (is.data.frame(ad)) {
    if ("raw_alpha" %in% colnames(ad)) {
      drop_vals <- ad[["raw_alpha"]]
    } else {
      drop_vals <- ad[[1]]
    }
  } else {
    drop_vals <- unlist(ad)
  }
  names(drop_vals) <- item_names

  #── 3. Prepare for CR: compute TotalScore & high/low groups
  tmp <- selected_items %>%
    dplyr::mutate(
      TotalScore = rowMeans(dplyr::across(dplyr::everything()), na.rm = TRUE)
    )

  n_valid <- sum(!is.na(tmp$TotalScore))
  rk <- rank(tmp$TotalScore, na.last = "keep", ties.method = "average")

  tmp <- tmp %>%
    dplyr::mutate(
      PctRank = rk / n_valid * 100,
      PerformanceGroup = dplyr::case_when(
        rk * 100 <= 27 * n_valid ~ "Low",
        rk * 100 >= 73 * n_valid ~ "High",
        TRUE ~ NA_character_
      )
    )

  #── 4. Descriptives
  desc_df <- psych::describe(selected_items) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Item") %>%
    dplyr::transmute(
      Item,
      Mean = sprintf("%.2f", mean),
      SD   = sprintf("%.2f", sd),
      Skew = sprintf(paste0("%.", digits, "f"), skew),
      Kurt = sprintf(paste0("%.", digits, "f"), kurtosis)
    )

  #── 5. CR via t.test + stars
  cr_df <- do.call(rbind, lapply(item_names, function(it) {
    high <- tmp %>% dplyr::filter(PerformanceGroup == "High") %>% dplyr::pull(.data[[it]])
    low  <- tmp %>% dplyr::filter(PerformanceGroup == "Low")  %>% dplyr::pull(.data[[it]])

    lev_df <- data.frame(
      value = c(high, low),
      group = factor(rep(c("High","Low"), c(length(high), length(low))))
    )
    p_lev   <- car::leveneTest(value ~ group, data = lev_df, center = mean)[1, "Pr(>F)"]
    eq_var  <- (p_lev > 0.05)

    tt   <- stats::t.test(high, low, var.equal = eq_var)
    stars <- if (tt$p.value < 0.001) "***" else
      if (tt$p.value < 0.01)  "**"  else
        if (tt$p.value < 0.05)  "*"   else ""
    t_fmt <- sprintf(paste0("%.", digits, "f"), unname(tt$statistic))
    data.frame(Item = it, CR = paste0(t_fmt, stars), stringsAsFactors = FALSE)
  }))

  #── 6. Alpha if deleted
  alpha_drop_df <- data.frame(
    Item               = item_names,
    `Alpha if deleted` = sprintf(paste0("%.", digits, "f"), drop_vals),
    stringsAsFactors   = FALSE
  )

  #── 7. CITC
  citc_df <- do.call(rbind, lapply(item_names, function(it) {
    others  <- setdiff(item_names, it)
    sum_oth <- rowSums(selected_items[, others, drop = FALSE], na.rm = TRUE)
    ct      <- stats::cor.test(selected_items[[it]], sum_oth)
    stars   <- if (ct$p.value < 0.001) "***" else
      if (ct$p.value < 0.01)  "**"  else
        if (ct$p.value < 0.05)  "*"   else ""
    r_fmt   <- sprintf(paste0("%.", digits, "f"), unname(ct$estimate))
    data.frame(Item = it, CITC = paste0(r_fmt, stars), stringsAsFactors = FALSE)
  }))

  #── 8. Merge everything
  final_df <- desc_df %>%
    dplyr::left_join(cr_df,         by = "Item") %>%
    dplyr::left_join(alpha_drop_df, by = "Item") %>%
    dplyr::left_join(citc_df,       by = "Item")

  #── 9. Add Total-alpha column if requested
  if (total) {
    tot <- sprintf(paste0("%.", digits, "f"), alpha_res$total$raw_alpha)
    final_df$Total <- c(tot, rep("", nrow(final_df) - 1))
  }

  #── 10. Clean up any NAs to blanks
  final_df[is.na(final_df)] <- ""

  #── Return
  final_df
}
