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
#'         \item \strong{Variable Range (Position-Based):} \code{var_name = "start_var:end_var"}
#'               Selects all columns from start_var to end_var based on their positions
#'               in the data frame, regardless of naming patterns. Supports reverse order.
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
#' \strong{Variable Range Method (Position-Based):} Use colon notation \code{"start:end"}
#' to select all columns between the positions of the start and end variables (inclusive).
#' The start and end variables do not need to share prefixes or follow any naming pattern.
#' Selection is based purely on column positions in the data frame. If start appears
#' after end, selection proceeds in reverse order.
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
#' @examples
#' library(dplyr)
#'
#' # Create sample data with mixed naming patterns
#' set.seed(123)
#' data <- data.frame(
#'   id = 1:100,
#'   item1 = rnorm(100, 3, 1),
#'   item2 = rnorm(100, 3, 1),
#'   height = rnorm(100, 170, 10),
#'   weight = rnorm(100, 70, 15),
#'   score_a = rnorm(100, 85, 10),
#'   score_b = rnorm(100, 88, 12),
#'   final_grade = rnorm(100, 90, 8),
#'   temperature = rnorm(100, 25, 5)
#' )
#'
#' # Position-based range selection (NEW FUNCTIONALITY)
#' # Works regardless of naming patterns
#' data_ranges <- data %>%
#'   compute(
#'     physical_measures = "height:score_b",    # Columns 4-6: height, weight, score_a, score_b
#'     mixed_range = "item1:weight",           # Columns 2-4: item1, item2, height, weight
#'     end_range = "score_a:temperature",      # Columns 5-8: score_a, score_b, final_grade, temperature
#'     mean, na.rm = TRUE
#'   )
#'
#' # Traditional pattern method still works
#' data_pattern <- data %>%
#'   compute(
#'     item_mean = "item{i}", 1:2,
#'     mean, na.rm = TRUE
#'   )
#'
#' # Explicit vector method
#' data_explicit <- data %>%
#'   compute(
#'     custom_score = c("height", "final_grade", "temperature"),
#'     mean, na.rm = TRUE
#'   )
#'
#' # Reverse order range selection
#' data_reverse <- data %>%
#'   compute(
#'     reverse_range = "temperature:score_a",  # Selects in reverse: temperature, final_grade, score_b, score_a
#'     mean, na.rm = TRUE
#'   )
#'
#' @export
compute <- function(.data, ...) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("dplyr package is required for this function")
  }

  # Helper function to parse varrange based on column positions
  parse_varrange <- function(range_str, data_cols) {
    if (!is.character(range_str) || length(range_str) != 1) {
      stop("varrange must be a single character string like 'start_var:end_var'")
    }

    parts <- strsplit(range_str, ":", fixed = TRUE)[[1]]
    if (length(parts) != 2) {
      stop("varrange must be in format 'start_var:end_var'")
    }

    start_var <- trimws(parts[1])
    end_var <- trimws(parts[2])

    # Check if both variables exist in the data
    missing <- setdiff(c(start_var, end_var), data_cols)
    if (length(missing) > 0) {
      stop(paste0("Variables not found in data: ", paste(missing, collapse = ", ")))
    }

    # Find positions of start and end variables
    start_pos <- match(start_var, data_cols)
    end_pos <- match(end_var, data_cols)

    # Create sequence (handles both forward and reverse order)
    step <- if (start_pos <= end_pos) 1 else -1
    positions <- seq(start_pos, end_pos, by = step)

    return(data_cols[positions])
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
    non_numeric <- !sapply(selected_data, is.numeric)
    if (any(non_numeric)) {
      non_numeric_names <- names(selected_data)[non_numeric]
      stop(paste("Non-numeric columns found:", paste(non_numeric_names, collapse = ", "),
                 "- all selected columns must be numeric"))
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
  var_names <- names(var_specs)  # Fixed: don't subset again

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

        } else if (grepl(":", current_spec, fixed = TRUE)) {
          # Position-based varrange method: "start_var:end_var"
          col_names <- parse_varrange(current_spec, names(.data))
          i <- i + 1

        } else {
          # Single variable name
          col_names <- current_spec
          i <- i + 1
        }
      } else {
        # Explicit vector of variable names
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
