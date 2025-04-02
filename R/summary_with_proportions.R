#' Summary with proportions
#' Title Fit indices and proportion of each class in LCA analysis
#'
#' @param res object returned by \code{readModels}.
#'
#' @returns a data.frame
#' @export
#'
#' @examples
#' summary_with_proportions(res_pre) %>%
#'       mutate(Title = str_replace(Title, " model;", "")) %>%
#'       kable(digits = 2, caption = "Table x. Model fit summary.") %>%
#'       kable_classic(html_font = "Cambria")
#' @importFrom bruceR cc

summary_with_proportions <- function(res) {
  # Retrieve the initial summary table
  summary_data <- SummaryTable(res, keepCols = cc("Title, AIC, BIC, aBIC, Entropy, T11_LMR_PValue, BLRT_PValue")) %>%
    rename(LMR_PValue = T11_LMR_PValue)

  # Extract and format the class proportions for each model
  proportions_data <- map(res, ~ {
    proportions <- .x[["class_counts"]][["mostLikely"]][["proportion"]]
    proportions <- sort(proportions, decreasing = TRUE)
    proportion_string <- paste(formatC(proportions * 100, format = "f", digits = 1), collapse = "/")
    return(proportion_string)
  })

  # Combine the proportion data with the summary table
  summary_data$"Proportion (%)" <- unlist(proportions_data)

  return(summary_data)
}
