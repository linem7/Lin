#' Plot the distribution of one or more questionnaire items
#'
#' `item_plot()` reshapes selected variables from **wide** to **long**
#' format and draws a vertically–stacked set of histograms, plus an
#' optional “Total” panel that shows the row-mean when multiple items
#' are supplied.  Facet strips are printed on the right (à la
#' \code{switch = "y"}), and a clean grey–grid style is obtained via
#' \code{\link[ggplot2]{theme_bw}()}.
#'
#' @param data    A data frame in **wide** format.
#' @param pattern A character string used to locate item columns.
#'                If it contains `"{i}"` it will be expanded with
#'                \code{glue::glue(pattern, i = indices)}.
#' @param indices Optional numeric or character vector that either
#'                indexes variables directly or supplies the values
#'                for `{i}` in \code{pattern}.
#' @param digits  Number of decimal places (reserved for future printed
#'                statistics).  Default is 2.
#' @param total   Logical.  If `TRUE` (default) and \eqn{n > 1} items
#'                are selected, a “Total” panel is appended showing the
#'                mean of the items for each respondent.
#' @param binwidth Width of the histogram bins.  Defaults to 1, which
#'                aligns nicely with Likert-type scales; change as
#'                needed for continuous data.
#'
#' @return A \pkg{ggplot2} object (use `print()` or just type its name
#'         to display the figure).
#' @export
#'
#' @examples
#' ## Example 1 – Likert items -------------------------------------------------
#' set.seed(1)
#' likert_dat <- data.frame(
#'   var1 = sample(1:5, 120, TRUE),
#'   var2 = sample(1:5, 120, TRUE),
#'   var3 = sample(1:5, 120, TRUE)
#' )
#' item_plot(likert_dat, "var{i}", 1:3)
#'
#' ## Example 2 – Continuous items (change `binwidth`) -------------------------
#' cont_dat <- mtcars[ , c("mpg", "disp", "hp")]
#' item_plot(cont_dat, pattern = "", indices = 1:3, binwidth = 5, total = FALSE)
item_plot <- function(data,
                      pattern  = NULL,
                      indices  = NULL,
                      digits   = 2,
                      total    = TRUE,
                      binwidth = 1) {

  ## --- 1. Resolve item names -------------------------------------------------
  if (!is.null(pattern) && grepl("\\{i\\}", pattern) && !is.null(indices)) {
    item_names <- as.character(glue::glue(pattern, i = indices))
  } else if (!is.null(pattern)) {
    item_names <- grep(pattern, names(data), value = TRUE)
  } else if (!is.null(indices)) {
    item_names <- names(data)[indices]
  } else {
    stop("Provide either 'pattern' or 'indices' (or both with '{i}' in pattern).")
  }

  if (length(item_names) == 0) {
    stop("No variables matched the supplied pattern/indices.")
  }

  ## --- 2. Create a trimmed data set -----------------------------------------
  selected_items <- dplyr::select(data, dplyr::all_of(item_names))

  ## --- 3. Wide → long --------------------------------------------------------
  long_df <- tidyr::pivot_longer(selected_items,
                                 tidyselect::everything(),
                                 names_to  = "item",
                                 values_to = "score")

  ## --- 4. Add total-score panel if requested --------------------------------
  if (total && length(item_names) > 1) {
    tot <- dplyr::mutate(selected_items,
                         Total = rowMeans(dplyr::across(tidyselect::everything()),
                                          na.rm = TRUE))
    long_df <- dplyr::bind_rows(
      long_df,
      dplyr::select(tot, Total) %>%
        dplyr::rename(score = Total) %>%
        dplyr::mutate(item = "Total")
    )
  }

  ## --- 5. Plot ---------------------------------------------------------------
  ggplot2::ggplot(long_df, ggplot2::aes(x = score)) +
    ggplot2::geom_histogram(color    = "black",
                            binwidth = binwidth,
                            boundary = 0.5) +
    ggplot2::facet_grid(item ~ .) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.01)))
}
