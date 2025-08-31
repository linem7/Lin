#' Compute New Variables by Aggregating Multiple Columns
#'
#' Creates new variables by applying aggregation functions (mean, sum, etc.) to
#' multiple columns using flexible variable selection methods. Supports three
#' approaches: pattern-based indexing, explicit variable vectors, and variable ranges.
#' Designed to work seamlessly with dplyr pipelines.
#'
#' @param .data A data frame or tibble containing the source variables
#' @param ... Variable specifications followed by aggregation function and options:
#'   \itemize{
#'     \item **Named variable specifications** using one of three methods:
#'       \itemize{
#'         \item \strong{Pattern + Index:} \code{var_name = "pattern{i}", indices}
#'               (e.g., \code{score1 = "item{i}", 1:5})
#'         \item \strong{Explicit Vector:} \code{var_name = c("var1", "var2", ...)}
#'               (e.g., \code{score2 = c("item1", "item3", "item5")})
#'         \item \strong{Variable Range:} \code{var_name = "start_var:end_var"}
#'               (e.g., \code{score3 = "item1:item5"})
#'       }
#'     \item **Aggregation function** as unnamed argument (e.g., \code{mean}, \code{sum}, \code{median})
#'     \item **Optional \code{na.rm}** logical argument for handling missing values
#'   }
#'
#' @details
#' **Variable Selection Methods:**
#'
#' \strong{Pattern Method:} Uses a pattern string containing \code{{i}} placeholder
#' followed by numeric indices. The pattern \code{"item{i}"} with indices \code{1:5}
#' expands to \code{c("item1", "item2", "item3", "item4", "item5")}.
#'
#' \strong{Explicit Vector Method:} Directly specify column names as a character vector.
#' Useful for non-sequential or custom variable selections.
#'
#' \strong{Variable Range Method:} Use colon notation \code{"start:end"} for sequential
#' variables with the same prefix and numeric suffixes. Automatically handles
#' zero-padded formats (e.g., \code{"q01:q05"}).
#'
#' **Performance Optimization:** For \code{mean} and \code{sum} functions with multiple
#' columns, the function uses optimized \code{rowMeans()} and \code{rowSums()}
#' implementations. Other functions use \code{apply()} row-wise computation.
#'
#' **Error Handling:** Missing columns generate warnings but don't stop execution.
#' All selected columns must be numeric. If no valid columns are found for a
#' specification, the result variable is filled with \code{NA}.
#'
#' @return A data frame identical to \code{.data} with additional computed variables
#'   appended as new columns. Original columns remain unchanged.
#'
#' @seealso \code{\link[base]{mean}}, \code{\link[base]{sum}},
#'   \code{\link[base]{rowMeans}}, \code{\link[base]{rowSums}},
#'   \code{\link[dplyr]{mutate}}
#'
#' @examples
#' library(dplyr)
#'
#' # Create sample data
#' set.seed(123)
#' data <- data.frame(
#'   id = 1:100,
#'   item1 = rnorm(100, 3, 1), item2 = rnorm(100, 3, 1), item3 = rnorm(100, 3, 1),
#'   item4 = rnorm(100, 3, 1), item5 = rnorm(100, 3, 1), item6 = rnorm(100, 4, 1),
#'   item7 = rnorm(100, 4, 1), item8 = rnorm(100, 4, 1), item9 = rnorm(100, 4, 1),
#'   item10 = rnorm(100, 4, 1), item11 = rnorm(100, 5, 1), item12 = rnorm(100, 5, 1),
#'   item13 = rnorm(100, 5, 1), item14 = rnorm(100, 5, 1), item15 = rnorm(100, 5, 1)
#' )
#'
#' # Add missing values for testing
#' data[1:5, c("item1", "item6", "item11")] <- NA
#'
#' # Method 1: Pattern + Index specification
#' data_means <- data %>%
#'   compute(
#'     scale1_mean = "item{i}", 1:5,
#'     scale2_mean = "item{i}", 6:10,
#'     scale3_mean = "item{i}", 11:15,
#'     mean, na.rm = TRUE
#'   )
#'
#' # Method 2: Explicit variable vector
#' data_custom <- data %>%
#'   compute(
#'     custom_score = c("item1", "item3", "item5", "item7"),
#'     mean, na.rm = TRUE
#'   )
#'
#' # Method 3: Variable range specification
#' data_ranges <- data %>%
#'   compute(
#'     early_items = "item1:item5",
#'     late_items = "item11:item15",
#'     mean, na.rm = TRUE
#'   )
#'
#' # Mixed methods in single call
#' data_mixed <- data %>%
#'   compute(
#'     pattern_mean = "item{i}", 1:3,              # Pattern method
#'     explicit_mean = c("item4", "item6", "item8"), # Explicit method
#'     range_sum = "item10:item12",                 # Range method
#'     mean, na.rm = TRUE
#'   )
#'
#' # Using different aggregation functions
#' data_sums <- data %>%
#'   compute(
#'     total1 = "item1:item5",
#'     total2 = c("item6", "item8", "item10"),
#'     sum, na.rm = FALSE
#'   )
#'
#' # Other statistical functions
#' data_stats <- data %>%
#'   compute(
#'     median_score = "item1:item5",
#'     median, na.rm = TRUE
#'   ) %>%
#'   compute(
#'     sd_score = c("item6", "item7", "item8"),
#'     sd, na.rm = TRUE
#'   )
#'
#' # Zero-padded variable names
#' padded_data <- data.frame(
#'   id = 1:50,
#'   q01 = rnorm(50), q02 = rnorm(50), q03 = rnorm(50),
#'   q04 = rnorm(50), q05 = rnorm(50)
#' )
#'
#' padded_result <- padded_data %>%
#'   compute(
#'     questionnaire_mean = "q01:q05",
#'     mean, na.rm = TRUE
#'   )
#'
#' # Non-sequential selections with pattern method
#' data_nonseq <- data %>%
#'   compute(
#'     odd_items = "item{i}", c(1, 3, 5, 7, 9),
#'     even_items = "item{i}", c(2, 4, 6, 8, 10),
#'     mean, na.rm = TRUE
#'   )
#'
#' @export
compute <- function(.data, ...) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for this function")
  }

  # Helper function to parse varrange (e.g., "item1:item5")
  parse_varrange <- function(range_str) {
    if (!is.character(range_str) || length(range_str) != 1) {
      stop("varrange must be a single character string like 'item1:item5'")
    }

    parts <- strsplit(range_str, ":", fixed = TRUE)[[1]]
    if (length(parts) != 2) {
      stop("varrange must be in format 'start_var:end_var'")
    }

    start_var <- trimws(parts[1])
    end_var <- trimws(parts[2])

    # Extract prefix and numeric parts using regex
    extract_parts <- function(var_name) {
      # Match: (prefix)(digits at end)
      match <- regexec("^(.*?)(\\d+)$", var_name)
      result <- regmatches(var_name, match)[[1]]

      if (length(result) != 3) {
        stop(paste("Cannot parse variable name:", var_name,
                   "- must end with numbers (e.g., 'item1')"))
      }

      list(
        prefix = result[2],
        num_str = result[3],
        num = as.integer(result[3])
      )
    }

    start_info <- extract_parts(start_var)
    end_info <- extract_parts(end_var)

    # Validate same prefix
    if (start_info$prefix != end_info$prefix) {
      stop("varrange endpoints must have the same prefix (e.g., 'item1:item5')")
    }

    # Validate numeric sequence
    if (start_info$num > end_info$num) {
      stop("Start number must be <= end number in varrange")
    }

    # Preserve zero-padding from the original format
    padding <- max(nchar(start_info$num_str), nchar(end_info$num_str))

    # Generate sequence
    num_sequence <- seq(start_info$num, end_info$num)
    formatted_nums <- sprintf(paste0("%0", padding, "d"), num_sequence)

    return(paste0(start_info$prefix, formatted_nums))
  }

  # Helper function to expand pattern with indices
  expand_pattern <- function(pattern, indices) {
    if (!is.character(pattern) || length(pattern) != 1) {
      stop("Pattern must be a single character string containing '{i}'")
    }
    if (!grepl("\\{i\\}", pattern)) {
      stop("Pattern must contain '{i}' placeholder")
    }

    sapply(indices, function(i) gsub("\\{i\\}", as.character(i), pattern))
  }

  # Helper function to apply aggregation function efficiently
  apply_aggregation <- function(data, col_names, func, na_rm) {
    selected_data <- data[, col_names, drop = FALSE]

    # Validate numeric columns
    if (!all(sapply(selected_data, is.numeric))) {
      stop("All selected columns must be numeric")
    }

    # Check if function supports na.rm
    func_formals <- names(formals(func))
    supports_na_rm <- "na.rm" %in% func_formals

    # Use optimized functions when possible
    if (ncol(selected_data) > 1) {
      if (identical(func, mean) && supports_na_rm) {
        return(rowMeans(selected_data, na.rm = na_rm))
      } else if (identical(func, sum) && supports_na_rm) {
        return(rowSums(selected_data, na.rm = na_rm))
      }
    }

    # General case using apply
    if (supports_na_rm) {
      return(apply(selected_data, 1, func, na.rm = na_rm))
    } else {
      if (na_rm) {
        warning("Function does not accept na.rm parameter; ignoring na.rm")
      }
      return(apply(selected_data, 1, func))
    }
  }

  # Parse arguments
  dots <- list(...)
  call_names <- names(dots)

  # Find the aggregation function (first unnamed function)
  func_idx <- NULL
  func_obj <- NULL

  for (i in seq_along(dots)) {
    if (is.null(call_names[i]) || call_names[i] == "") {
      arg_val <- dots[[i]]

      # Check if it's a function
      if (is.function(arg_val)) {
        func_obj <- arg_val
        func_idx <- i
        break
      } else if (is.name(arg_val) && exists(as.character(arg_val))) {
        potential_func <- get(as.character(arg_val))
        if (is.function(potential_func)) {
          func_obj <- potential_func
          func_idx <- i
          break
        }
      }
    }
  }

  if (is.null(func_idx)) {
    stop("Please provide a function (e.g., mean, sum) as an unnamed argument")
  }

  # Extract na.rm if provided
  na_rm_idx <- which(call_names == "na.rm")
  na_rm <- if (length(na_rm_idx) > 0) dots[[na_rm_idx[1]]] else FALSE

  # Remove function and na.rm from variable specifications
  remove_idx <- c(func_idx, na_rm_idx)
  var_specs <- dots[-remove_idx]
  var_names <- names(var_specs)[-remove_idx]

  # Process each variable specification
  result <- .data
  i <- 1

  while (i <= length(var_specs)) {
    var_name <- var_names[i]
    if (is.null(var_name) || var_name == "") {
      stop("Each variable specification must be named (e.g., var1 = ...)")
    }

    current_spec <- var_specs[[i]]
    col_names <- NULL

    # Determine specification type and generate column names
    if (is.character(current_spec)) {
      if (length(current_spec) == 1) {
        if (grepl("\\{i\\}", current_spec)) {
          # Pattern + index method: "item{i}", 1:5
          if (i + 1 > length(var_specs)) {
            stop(paste("Pattern", current_spec, "for", var_name,
                       "must be followed by indices"))
          }

          indices <- var_specs[[i + 1]]
          if (!is.numeric(indices)) {
            stop(paste("Indices for", var_name, "must be numeric"))
          }

          col_names <- expand_pattern(current_spec, indices)
          i <- i + 2  # Skip both pattern and indices

        } else if (grepl(":", current_spec)) {
          # Varrange method: "item1:item5"
          col_names <- parse_varrange(current_spec)
          i <- i + 1

        } else {
          # Single variable name
          col_names <- current_spec
          i <- i + 1
        }
      } else {
        # Vars method: c("item1", "item2", "item3")
        col_names <- current_spec
        i <- i + 1
      }
    } else {
      stop(paste("Invalid specification for", var_name,
                 "- must be character pattern, varrange, or vars"))
    }

    # Validate and compute
    existing_cols <- intersect(col_names, names(.data))
    missing_cols <- setdiff(col_names, names(.data))

    if (length(missing_cols) > 0) {
      warning(paste("For", var_name, "- missing columns:",
                    paste(missing_cols, collapse = ", ")))
    }

    if (length(existing_cols) == 0) {
      warning(paste("No valid columns found for", var_name, "- filling with NA"))
      result[[var_name]] <- NA_real_
    } else {
      result[[var_name]] <- apply_aggregation(.data, existing_cols, func_obj, na_rm)
    }
  }

  return(result)
}
