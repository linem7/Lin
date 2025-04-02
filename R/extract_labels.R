#' Extract Variable Labels from a Data Frame
#'
#' This function retrieves the variable labels from each column of a data frame.
#' It expects that each column of the data frame has a "label" attribute (as used in imported datasets or assigned manually).
#'
#' @param df A data frame (or tibble) whose columns have \code{label} attributes.
#' @return A tibble with two columns: \code{variable} (the column name) and \code{label} (the label value).
#' @examples
#' # Create a sample data frame with labels
#' df <- data.frame(id = 1:3, score = c(5, 4, 8))
#' attr(df$id, "label") <- "ID Number"
#' attr(df$score, "label") <- "Score"
#' extract_labels(df)
#' @export
extract_labels <- function(data) {
  data %>%
    map(~ attr(.x, "label")) %>%
    enframe(name = "Variable", value = "Label") %>%
    filter(!is.na(Label))  # Keep only variables that have labels
}
