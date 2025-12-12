#' Generate Variable Names from a Pattern
#'
#' A helper function to generate a vector of variable names by expanding a
#' pattern string with a numeric sequence.
#'
#' @param str A character string containing the placeholder \code{"{i}"}.
#' @param n A numeric vector of indices to substitute into the placeholder.
#'
#' @return A character vector of expanded variable names.
#' @export
pattern <- function(str, n) {
  if (!grepl("\\{i\\}", str)) {
    stop("Pattern must contain '{i}' placeholder")
  }
  sapply(n, function(i) gsub("\\{i\\}", as.character(i), str))
}

#' Compute New Variables by Aggregating Multiple Columns
#'
#' Creates new variables by applying aggregation functions (e.g., \code{mean}, \code{sum})
#' to multiple columns. It leverages \code{tidyselect} syntax for flexible variable
#' selection, including ranges, specific names, and a custom \code{pattern()} helper.
#'
#' @param .data A data frame or tibble containing the source variables.
#' @param ... Named arguments defining new variables to create. The value of each
#'   argument should be a selection of columns using:
#'   \itemize{
#'     \item \strong{Standard Tidyselect:} \code{var1:var5}, \code{starts_with("item")}, \code{c(var1, var2)}.
#'     \item \strong{Pattern Helper:} \code{pattern("item{i}", 1:5)} for constructing names dynamically.
#'   }
#' @param .method The aggregation function to apply (e.g., \code{mean}, \code{sum}).
#'   Defaults to \code{mean}.
#' @param .na.rm Logical. Should missing values be removed during aggregation?
#'   Defaults to \code{FALSE}.
#'
#' @details
#' \strong{Variable Selection:}
#' This function supports full \code{dplyr}-style selection. You can use:
#' \itemize{
#'   \item \code{item1:item5} (ranges)
#'   \item \code{c(item1, item3, item5)} (specific columns)
#'   \item \code{starts_with("A")}, \code{contains("scale")} (helpers)
#'   \item \code{pattern("Q{i}_rev", 1:5)} (custom pattern helper)
#' }
#'
#' \strong{Performance Optimization:}
#' If \code{.method} is set to \code{mean} or \code{sum}, the function automatically
#' switches to the highly optimized \code{rowMeans()} or \code{rowSums()} base R functions.
#' For other functions (e.g., \code{sd}, \code{median}), it falls back to \code{apply()},
#' which may be slower on large datasets.
#'
#' \strong{Error Handling:}
#' \itemize{
#'   \item If a selection matches no columns, the new variable is created and filled with \code{NA},
#'   and a warning is issued.
#'   \item If non-numeric columns are selected, an error will occur during computation.
#' }
#'
#' @return A data frame identical to \code{.data} with the new computed variables appended.
#'
#' @examples
#' library(dplyr)
#'
#' # Create sample data
#' df <- data.frame(
#'   id = 1:5,
#'   A1 = c(1, 2, 3, 4, 5),
#'   A2 = c(2, 2, 3, 4, 4),
#'   A3 = c(1, 1, 1, 1, 1),
#'   B_1 = c(5, 4, 3, 2, 1),
#'   B_2 = c(5, 5, 4, 4, 3),
#'   chk_1 = c(1, 0, 1, 0, 1),
#'   chk_2 = c(0, 1, 0, 1, 0)
#' )
#'
#' # 1. Basic Range Selection (Standard Tidyselect)
#' result <- df %>%
#'   compute(
#'     scale_A = A1:A3,
#'     .method = mean
#'   )
#'
#' # 2. Using the pattern() helper for complex names
#' result <- df %>%
#'   compute(
#'     scale_B = pattern("B_{i}", 1:2),
#'     .method = sum
#'   )
#'
#' # 3. Combining methods + na.rm
#' result <- df %>%
#'   compute(
#'     scale_A = A1:A3,
#'     scale_B = pattern("B_{i}", 1:2),
#'     check_score = starts_with("chk"),
#'     .method = mean,
#'     .na.rm = TRUE
#'   )
#' @importFrom rlang enquos
#' @importFrom tidyselect eval_select
#' @export
compute <- function(.data, ..., .method = mean, .na.rm = FALSE) {

  # 1. Capture the '...' arguments (the new variables to create)
  #    This allows us to handle input like: scale_a = item1:item5
  dots <- enquos(...)

  # 2. Initialize result as the original data
  result <- .data

  # 3. Check if the method supports rowMeans/rowSums optimization
  #    (This makes the function much faster than standard apply)
  use_optimized <- FALSE
  if (identical(.method, mean)) {
    op_func <- rowMeans
    use_optimized <- TRUE
  } else if (identical(.method, sum)) {
    op_func <- rowSums
    use_optimized <- TRUE
  }

  # 4. Iterate through each new variable definition
  for (new_var_name in names(dots)) {

    # Extract the logic provided by the user (e.g., item1:item5)
    current_quo <- dots[[new_var_name]]

    # MAGICAL STEP: use tidyselect to find column positions
    # This automatically handles:
    # - Names: c("a", "b")
    # - Ranges: item1:item5
    # - Helpers: starts_with("x")
    # - Our helper: use_pattern("x{i}", 1:3)
    tryCatch({
      col_pos <- tidyselect::eval_select(current_quo, data = result)
    }, error = function(e) {
      stop(paste("Error selecting variables for:", new_var_name, "\n", e$message))
    })

    # If no columns selected, fill with NA and warn
    if (length(col_pos) == 0) {
      warning(paste("No columns matched for", new_var_name))
      result[[new_var_name]] <- NA_real_
      next
    }

    # Subset the data for calculation
    # We use base R subsetting [ , col_pos] for speed
    selected_data <- result[, col_pos, drop = FALSE]

    # 5. Calculate the score
    if (use_optimized) {
      # Fast path for Mean/Sum
      result[[new_var_name]] <- op_func(selected_data, na.rm = .na.rm)
    } else {
      # Slower path for custom functions (e.g., sd, max)
      # We verify if the function accepts na.rm to avoid errors
      if ("na.rm" %in% names(formals(.method))) {
        result[[new_var_name]] <- apply(selected_data, 1, .method, na.rm = .na.rm)
      } else {
        result[[new_var_name]] <- apply(selected_data, 1, .method)
      }
    }
  }

  return(result)
}
