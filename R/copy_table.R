#' Copy the Last Printed Data Frame to Clipboard
#'
#' Copies the most recently evaluated data frame or matrix
#' to the system clipboard in tab-separated format, suitable
#' for pasting into spreadsheet applications.
#'
#' @param row.names Logical. Should row names be included? Default is FALSE.
#' @param col.names Logical. Should column names be included? Default is TRUE.
#' @param sep Field separator string. Default is tab character.
#' @param ... Additional arguments passed to \code{\link{write.table}}.
#'
#' @return Invisibly returns the copied object.
#'
#' @examples
#' \dontrun{
#' head(mtcars)
#' copy_table()
#' }
#'
#' @export
copy_table <- function(row.names = FALSE, col.names = TRUE, sep = "\t", ...) {
  obj <- .Last.value
  if (is.matrix(obj)) obj <- as.data.frame(obj)
  if (!is.data.frame(obj)) stop("Last value is not a data.frame or matrix.")

  os <- Sys.info()["sysname"]
  if (os == "Windows") {
    con <- "clipboard"
  } else if (os == "Darwin") {
    con <- pipe("pbcopy", "w")
  } else {
    con <- pipe("xclip -selection clipboard", "w")
  }

  write.table(obj, con, sep = sep, row.names = row.names,
              col.names = col.names, quote = FALSE, ...)
  if (inherits(con, "connection")) close(con)

  invisible(obj)
}
