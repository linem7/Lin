#' Extract Variable Labels from a Data Frame
#'
#' This function retrieves the variable labels from each column of a data frame.
#' It expects that each column of the data frame has a "label" attribute (as used in imported datasets or assigned manually).
#'
#' @param df A data frame (or tibble) whose columns have \code{label} attributes.
#' @return A tibble with two columns: \code{variable} (the column name) and \code{label} (the label value).
#' @importFrom tibble enframe
#' @importFrom dplyr filter
#' @importFrom purrr map_chr
#' @examples
#' # Create a sample data frame with labels
#' df <- data.frame(id = 1:3, score = c(5, 4, 8))
#' attr(df$id, "label") <- "ID Number"
#' attr(df$score, "label") <- "Score"
#' extract_labels(df)
#'
#' @export
extract_labels <- function(data) {
  data %>%
    purrr::map_chr(~ attr(.x, "label") %||% NA_character_) %>%
    tibble::enframe(name = "Variable", value = "Label") %>%
    dplyr::filter(!is.na(Label))
}
