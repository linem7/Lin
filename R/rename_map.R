#' Bulk-rename Columns
#'
#' Rename columns in a \code{data.frame} or tibble by replacing specified old prefixes
#' with new prefixes.  Only column names that start with an \code{old} prefix will be
#' modified; all other names remain unchanged.
#'
#' @param data A \code{data.frame} or tibble whose column names will be modified.
#' @param map  A named character vector: use \strong{new_name = "old_name"} (old name should be quoted)
#'
#' @return A \code{data.frame} (or tibble) with columns renamed according to \code{map}.
#'
#' @examples
#' # Create a sample data.frame
#' df <- data.frame(
#'   measureA1 = rnorm(5),
#'   measureA2 = rnorm(5),
#'   factorX1  = runif(5),
#'   factorX2  = runif(5),
#'   Timestamp = Sys.time() + 1:5,
#'   stringsAsFactors = FALSE
#' )
#'
#' # Define mapping: old → new prefixes
#' prefix_map <- c(
#'   measureA = "mA",
#'   factorX  = "fX"
#' )
#'
#' # Apply rename_map()
#' df_renamed <- rename_map(df, prefix_map)
#' names(df_renamed)
#' #> [1] "mA1"      "mA2"      "fX1"      "fX2"      "Timestamp"
#'
#' @export
rename_map <- function(data, map) {
  nm <- names(data)
  for (old_pref in names(map)) {
    new_pref <- map[[old_pref]]
    nm <- sub(paste0("^", old_pref), new_pref, nm)
  }
  names(data) <- nm
  data
}
