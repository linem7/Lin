#' APA a table
#'
#' Using kable and kable_classic to create a apa format table
#'
#' @param x data.frame
#' @param title title of table
#' @param notes footnote of table
#' @param digits digits of the value within the table
#' @return a data.frame
#'
#' @examples
#' # Format the first few rows of mtcars as an APA style table
#' apa(head(mtcars, 5), title = "Descriptive Statistics for mtcars Sample")
#' @importFrom kableExtra footnote
#'
#'@export
apa <- function(x, title = NULL, notes = NULL, digits = 3){
  x %>% knitr::kable(digits = digits, caption = title) %>%
    kableExtra::kable_classic(html_font = "Cambria", position = "left") %>%
    kableExtra::footnote(general = notes, footnote_as_chunk = T)
}
