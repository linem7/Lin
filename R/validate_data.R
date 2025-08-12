#' Validate Data Against Allowed Values
#'
#' @description
#' Checks a data frame for values that fall outside a specified set of allowed values.
#' Returns the locations and values of any invalid entries found.
#'
#' @param data A data frame to validate
#' @param allowed_values A numeric vector or single value specifying which values are considered valid
#' @param allow_missing Logical. If TRUE (default), NA values are considered valid.
#'   If FALSE, NA values are reported as invalid.
#'
#' @return A data frame with columns \code{column}, \code{row}, and \code{value} containing
#'   the locations and values of invalid entries. Returns \code{invisible(NULL)} if all
#'   values are valid.
#'
#' @details
#' This function iterates through each column of the input data frame and identifies
#' values that are not in the set of allowed values. It's particularly useful for
#' validating binary data (0/1), categorical data with known levels, or checking
#' for unexpected values in processed datasets.
#'
#' @examples
#' # Create sample data
#' set.seed(42)
#' df <- data.frame(
#'   treatment = sample(c(0, 1), 20, replace = TRUE),
#'   response = sample(c(0, 1), 20, replace = TRUE),
#'   group = sample(c(1, 2, 3), 20, replace = TRUE)
#' )
#'
#' # Introduce some invalid values
#' df[5, "treatment"] <- 2
#' df[10, "response"] <- -1
#' df[15, "group"] <- 99
#' df[3, "treatment"] <- NA
#'
#' # Example 1: Check for binary values (0, 1) allowing missing values
#' invalid_binary <- validate_data(df[, c("treatment", "response")], c(0, 1))
#' invalid_binary
#'
#' # Example 2: Check for binary values, treating NA as invalid
#' invalid_strict <- validate_data(df[, c("treatment", "response")], c(0, 1),
#'                                  allow_missing = FALSE)
#' invalid_strict
#'
#' # Example 3: Check categorical variable for expected levels
#' invalid_groups <- validate_data(df[, "group", drop = FALSE], c(1, 2, 3))
#' invalid_groups
#'
#' # Example 4: Check if all values are valid (returns NULL silently)
#' clean_df <- data.frame(
#'   a = c(0, 1, 0, 1, NA),
#'   b = c(1, 1, 0, 0, 1)
#' )
#' result <- validate_data(clean_df, c(0, 1))
#' is.null(result)  # TRUE
#'
#' # Example 5: Use in data validation pipeline
#' \dontrun{
#' # Load data
#' my_data <- read.csv("data.csv")
#'
#' # Validate binary columns
#' binary_cols <- c("is_treatment", "has_disease", "is_male")
#' invalid <- validate_data(my_data[binary_cols], c(0, 1))
#'
#' if (!is.null(invalid)) {
#'   warning(paste("Found", nrow(invalid), "invalid values in binary columns"))
#'   print(invalid)
#' }
#' }
#'
#' @export
#' @seealso
#' \code{\link{is.na}} for checking missing values,
#' \code{\link{which}} for finding positions of values
#'
#' @importFrom stats na.omit
validate_data <- function(data, allowed_values, allow_missing = TRUE) {
  # Convert single number to vector if needed
  if(length(allowed_values) == 1) {
    allowed_values <- c(allowed_values)
  }

  invalid_df <- data.frame()

  for(col_name in names(data)) {
    col_data <- data[[col_name]]

    if(allow_missing) {
      # Check for values that are not in allowed_values and not NA
      invalid_idx <- which(!is.na(col_data) & !(col_data %in% allowed_values))
    } else {
      # Check for values that are not in allowed_values (including NA as invalid)
      invalid_idx <- which(!(col_data %in% allowed_values))

      # Also find NA locations
      na_idx <- which(is.na(col_data))
      if(length(na_idx) > 0) {
        temp_df <- data.frame(
          column = col_name,
          row = na_idx,
          value = "NA",
          stringsAsFactors = FALSE
        )
        invalid_df <- rbind(invalid_df, temp_df)
      }
    }

    if(length(invalid_idx) > 0) {
      temp_df <- data.frame(
        column = col_name,
        row = invalid_idx,
        value = as.character(col_data[invalid_idx]),
        stringsAsFactors = FALSE
      )
      invalid_df <- rbind(invalid_df, temp_df)
    }
  }

  if(nrow(invalid_df) == 0) {
    return(invisible(NULL))
  } else {
    return(invalid_df)
  }
}
