#' Correlation table
#'
#' Function for creating correlation matrix including Mean and SD.
#'
#' @param data A data.frame
#' @param vars variables included in correlation analysis, use \code{select} from dplyr to choose columns
#' @param digits_stat digits in the mean and sd analysis
#' @param digits_corr digits in the correlation coefficients
#' @param stats_position should be mean and standard deviation placed in the first two columns or the bottom of the table
#'
#' @returns a data.frame
#' @export
#'
#' @examples
#' cor_data <- select(data, gender, age, ses, var1, var2, var3)
#'
apa_corr_matrix <- function(data, vars, digits_stat = 2, digits_corr = 3, stats_position = c("columns", "rows")) {
  # Match the argument for stats_position
  stats_position <- match.arg(stats_position)

  # Load necessary package for tidy selection
  library(dplyr)

  # Subset the data using tidyselect syntax
  data_sel <- dplyr::select(data, {{ vars }})
  var_names <- names(data_sel)
  n <- length(var_names)

  # Helper function to determine significance stars
  get_stars <- function(p_val) {
    if (is.na(p_val)) return("")
    if (p_val < 0.001) return("***")
    if (p_val < 0.01)  return("**")
    if (p_val < 0.05)  return("*")
    return("")
  }

  # Compute means and standard deviations, formatting with sprintf
  means <- sapply(data_sel, mean, na.rm = TRUE)
  sds   <- sapply(data_sel, sd, na.rm = TRUE)

  # Initialize matrices for correlations and p-values
  corr_matrix <- matrix(NA, n, n)
  p_matrix    <- matrix(NA, n, n)

  # Compute correlations and corresponding p-values (only for lower triangle)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i > j) {
        test <- cor.test(data_sel[[i]], data_sel[[j]])
        corr_matrix[i, j] <- test$estimate
        p_matrix[i, j]    <- test$p.value
      }
    }
  }

  # Create a formatted correlation matrix as a character matrix
  formatted_corr <- matrix("", n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i > j) {
        r_val <- sprintf(paste0("%.", digits_corr, "f"), corr_matrix[i, j])
        stars <- get_stars(p_matrix[i, j])
        formatted_corr[i, j] <- paste0(r_val, stars)
      } else if (i == j) {
        formatted_corr[i, j] <- "-"  # Diagonal is replaced with a dash
      } else {
        formatted_corr[i, j] <- ""
      }
    }
  }

  # Build the final output depending on stats_position
  if (stats_position == "columns") {
    # Build output with descriptive statistics as the first two columns.
    corr_df <- as.data.frame(formatted_corr, stringsAsFactors = FALSE)
    colnames(corr_df) <- as.character(seq_len(n))

    stats_df <- data.frame(
      M  = sprintf(paste0("%.", digits_stat, "f"), means),
      SD = sprintf(paste0("%.", digits_stat, "f"), sds),
      row.names = var_names,
      stringsAsFactors = FALSE
    )

    final_df <- cbind(stats_df, corr_df)
    rownames(final_df) <- paste0(seq_len(n), ". ", var_names)

  } else if (stats_position == "rows") {
    # Build output with descriptive statistics as the last two rows.
    # Create a matrix with (n + 2) rows: n rows for variables and 2 extra for M and SD;
    # and (n + 1) columns: 1 for variable labels and n for correlations.
    final_table <- matrix("", nrow = n + 2, ncol = n + 1)
    colnames(final_table) <- c("Variable", as.character(seq_len(n)))

    # Fill in the variable rows with labels and correlation values.
    for (i in seq_len(n)) {
      final_table[i, 1] <- paste0(i, ". ", var_names[i])
      for (j in seq_len(n)) {
        final_table[i, j + 1] <- formatted_corr[i, j]
      }
    }

    # Append the mean (M) and standard deviation (SD) as the last two rows.
    final_table[n + 1, 1] <- "M"
    final_table[n + 2, 1] <- "SD"
    for (j in seq_len(n)) {
      final_table[n + 1, j + 1] <- sprintf(paste0("%.", digits_stat, "f"), means[j])
      final_table[n + 2, j + 1] <- sprintf(paste0("%.", digits_stat, "f"), sds[j])
    }

    final_df <- as.data.frame(final_table, stringsAsFactors = FALSE)
  }

  return(final_df)
}
