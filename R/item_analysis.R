#' Item Analysis for Classical Test Theory
#'
#' Computes a variety of item‐level statistics for scales (e.g., means, SDs,
#' skewness, kurtosis, extreme‐group discrimination, Cronbach’s α if item deleted,
#' and corrected item–total correlations with significance stars), plus overall
#' reliability.
#'
#' @param data A data frame containing item response columns.
#' @param pattern A glue‐style pattern for item names (e.g. `"item{i}"`).
#'   The placeholder `i` will be filled by each value in `indices`.
#' @param indices An integer vector of indices to substitute into `pattern`.
#'   For example, `pattern = "item{i}", indices = 1:5` selects `"item1"` … `"item5"`.
#' @param digits Integer; number of decimal places for skewness, kurtosis, CR, CITC,
#'   and alpha‐if‐deleted values. Means and SDs are always formatted to two decimals.
#'   Default is 3.
#' @param total Logical; if `TRUE`, prepends a “Total” column with the overall
#'   Cronbach’s α in the first row and blanks elsewhere. Default is `TRUE`.
#' @param ... Additional arguments passed to `psych::alpha()` (e.g.
#'   `check.keys = TRUE` to auto‐flip negatively keyed items).
#'
#' @return A data frame with one row per item and the following columns:
#' \item{Mean}{Item mean score.}
#' \item{SD}{Standard deviation of the item.}
#' \item{Skew}{Skewness of the item distribution.}
#' \item{Kurt}{Kurtosis of the item distribution.}
#' \item{CR}{t-value from high vs. low group discrimination, with significance stars.}
#' \item{Alpha if deleted}{Cronbach’s α when the item is removed.}
#' \item{CITC}{Corrected item–total correlation, with significance stars.}
#' \item{Total (optional)}{Overall α in the first row; blanks in remaining rows.}
#'
#' @importFrom dplyr select mutate left_join
#' @importFrom tibble rownames_to_column
#' @importFrom psych describe
#' @importFrom glue glue
#'
#' @details
#' The function constructs item names using \code{glue::glue()} with the provided \code{pattern} and \code{indices}, then selects those columns from the input \code{data}.
#'
#' To assess internal consistency, it computes overall reliability using \code{psych::alpha()} and extracts the \code{alpha.drop} table to obtain Cronbach’s α if each item were deleted. Total scores are calculated by summing across all selected items. Based on these scores, the function defines high and low performance groups using the 27th and 73rd percentiles, respectively.
#'
#' Descriptive statistics, including mean, standard deviation, skewness, and kurtosis, are calculated with \code{psych::describe()}. Means and SDs are formatted to two decimal places, while skewness and kurtosis are rounded based on the \code{digits} argument.
#'
#' Discrimination (CR) is assessed through two-sample t-tests (\code{t.test()} with \code{var.equal = TRUE}) comparing high and low groups. t-values are formatted and annotated with significance stars. The corrected item–total correlation (CITC) is computed by correlating each item with the total score excluding that item, and formatted similarly with stars based on p-values.
#'
#' Finally, all results—including descriptives, CR, α if deleted, CITC, and optionally the overall α—are merged into a tidy data frame. And the overall α can be displayed in a “Total” column in the first row if requested.
#' @examples
#' data(good_rel)
#' data(poor_rel)
#'
#' # Good‐reliability example
#' # good_rel: a data.frame where items 1–5 all load positively on one factor
#' item_analysis(good_rel, pattern = "item{i}", indices = 1:5)
#'
#' # Use Lin::apa() format the final output in html
#' item_analysis(good_rel, pattern = "item{i}", indices = 1:5) %>% apa()#
#'
#' # Poor‐reliability example
#' # poor_rel: a data.frame with random noise or mixed floor/ceiling effects
#' item_analysis(poor_rel, pattern = "item{i}", indices = 1:5, check.keys = TRUE)
#'
#' @export

item_analysis <- function(data,
                          pattern,
                          indices,
                          digits = 3,
                          total  = TRUE,
                          ...) {
  #── Dependencies
  if (!requireNamespace("dplyr", quietly = TRUE))  stop("Install dplyr")
  if (!requireNamespace("psych", quietly = TRUE))  stop("Install psych")
  if (!requireNamespace("glue", quietly = TRUE))   stop("Install glue")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Install tibble")

  #── 1. Build and select item columns
  item_names     <- as.character(glue::glue(pattern, i = indices))
  selected_items <- data %>% dplyr::select(dplyr::all_of(item_names))

  #── 2. Compute total α (with any extra args) and get alpha.drop
  psych_alpha <- psych::alpha
  alpha_res       <- psych_alpha(selected_items, ...)
  ad <- alpha_res$alpha.drop
  if (is.data.frame(ad)) {
    # prefer a column named “raw_alpha” if it exists
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
      TotalScore = rowSums(across(everything()), na.rm = TRUE),
      PerformanceGroup = dplyr::case_when(
        TotalScore <= quantile(TotalScore, 0.27, na.rm = TRUE) ~ "Low",
        TotalScore >= quantile(TotalScore, 0.73, na.rm = TRUE) ~ "High",
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
    high <- tmp %>% filter(PerformanceGroup == "High") %>% pull(.data[[it]])
    low  <- tmp %>% filter(PerformanceGroup == "Low")  %>% pull(.data[[it]])

    # prepare for Levene’s test
    lev_df <- data.frame(
      value = c(high, low),
      group = factor(rep(c("High","Low"), c(length(high), length(low))))
    )
    p_lev   <- car::leveneTest(value ~ group, data = lev_df)[1, "Pr(>F)"]
    eq_var  <- (p_lev > 0.05)

    # run the matching t‐test
    tt   <- t.test(high, low, var.equal = eq_var)
    stars <- if (tt$p.value < 0.001) "***" else
      if (tt$p.value < 0.01)  "**"  else
        if (tt$p.value < 0.05)  "*"   else ""
    t_fmt <- sprintf(paste0("%.", digits, "f"), unname(tt$statistic))
    data.frame(Item = it, CR = paste0(t_fmt, stars), stringsAsFactors = FALSE)
  }))

  #── 6. Alpha if deleted (from psych::alpha)
  alpha_drop_df <- data.frame(
    Item               = item_names,
    `Alpha if deleted` = sprintf(paste0("%.", digits, "f"), drop_vals),
    stringsAsFactors   = FALSE
  )

  #── 7. CITC: correlate each item with sum of the others + stars
  citc_df <- do.call(rbind, lapply(item_names, function(it) {
    others  <- setdiff(item_names, it)
    sum_oth <- rowSums(selected_items[, others, drop = FALSE], na.rm = TRUE)
    ct      <- cor.test(selected_items[[it]], sum_oth)
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

  #── 9. Add Total‐alpha column if requested
  if (total) {
    tot <- sprintf(paste0("%.", digits, "f"), alpha_res$total$raw_alpha)
    final_df$Total <- c(tot, rep("", nrow(final_df) - 1))
  }

  #── 10. Clean up any NAs → blanks
  final_df[is.na(final_df)] <- ""

  #── Return
  final_df
}
