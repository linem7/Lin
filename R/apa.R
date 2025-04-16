#' Create an APA-Style Table
#'
#' Formats a data frame as an APA-style table using \code{kable} and \code{kable_classic}.
#' Optionally adds a title and footnote, and controls the number of displayed digits.
#'
#' @param x A \code{data.frame}. The data to be formatted as a table.
#' @param title \code{character}. Optional. The title (caption) of the table.
#' @param notes \code{character}. Optional. Footnote text to display below the table.
#' @param digits \code{integer}. Number of decimal places to display for numeric columns. Default is 3.
#'
#' @return An object of class \code{kableExtra}, representing the formatted table. Typically rendered in RMarkdown or HTML output.
#'
#' @details
#' This function uses \code{knitr::kable()} to create a table and \code{kableExtra::kable_classic()} to apply APA-style formatting.
#' A footnote can be added using \code{kableExtra::footnote()}.
#'
#' @examples
#' # Example 1: Format the first five rows of mtcars as an APA-style table
#' apa(
#'   head(mtcars, 5),
#'   title = "Descriptive Statistics for mtcars Sample",
#'   notes = "Note. Data are from the mtcars dataset.",
#'   digits = 2
#' )
#'
#' # Example 2: Format iris summary statistics
#' apa(
#'   summary(iris),
#'   title = "Summary Statistics for Iris Dataset",
#'   notes = "Note. Classic Fisher's iris data.",
#'   digits = 1
#' )
#'
#' @importFrom knitr kable
#' @importFrom kableExtra kable_classic footnote
#' @export
apa <- function(x, title = NULL, notes = NULL, digits = 3){
  x %>%
    kable(digits = digits, caption = title) %>%
    kable_classic(html_font = "Cambria", position = "left") %>%
    footnote(general = notes, footnote_as_chunk = TRUE)
}
