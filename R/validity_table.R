# Helper function to determine estimator family
get_estimator_family <- function(fit) {
  est <- tryCatch(lavaan::lavInspect(fit, "options")$estimator,
                  error = function(e) NA_character_)
  robust_ests <- c("MLM", "MLMVS", "MLMV", "MLR", "WLSMV", "ULSMV", "WLSM", "DWLS")
  if (is.na(est)) return("unknown")
  if (est %in% robust_ests) "robust" else "ml"
}

# Helper function to select best available fit measure
get_best_measure <- function(meas, name, family) {
  candidates <- switch(name,
                       "chisq" = if (family == "robust") c("chisq.robust", "chisq.scaled", "chisq")
                       else c("chisq.scaled", "chisq"),
                       "df" = if (family == "robust") c("df.robust", "df.scaled", "df")
                       else c("df.scaled", "df"),
                       "cfi" = if (family == "robust") c("cfi.robust", "cfi.scaled", "cfi")
                       else c("cfi.scaled", "cfi"),
                       "tli" = if (family == "robust") c("tli.robust", "tli.scaled", "tli")
                       else c("tli.scaled", "tli"),
                       "rmsea" = if (family == "robust") c("rmsea.robust", "rmsea.scaled", "rmsea")
                       else c("rmsea.scaled", "rmsea"),
                       "srmr" = c("srmr"),
                       name
  )

  hit <- candidates[candidates %in% names(meas)]
  used <- if (length(hit)) hit[1] else name
  val <- if (used %in% names(meas)) as.numeric(meas[[used]]) else NA_real_
  list(value = val, used = used)
}

# Helper function to format column names for publication
format_column_name <- function(x) {
  x <- tolower(x)
  switch(x,
         "chisq" = "χ²",
         "df" = "df",
         "cfi" = "CFI",
         "tli" = "TLI",
         "rmsea" = "RMSEA",
         "srmr" = "SRMR",
         toupper(x)
  )
}

#' Generate Comparative Tables of lavaan Model Fit Indices
#'
#' @description
#' Creates publication-ready tables comparing fit indices across multiple
#' structural equation models fitted with lavaan. Accepts either pre-fitted
#' lavaan objects or model specifications that will be automatically fitted
#' using confirmatory factor analysis (CFA). Automatically selects appropriate
#' fit measures based on estimator type and provides delta (change) columns
#' for model comparison.
#'
#' @param ... Either fitted lavaan objects (from \code{cfa()}, \code{sem()}, etc.)
#'   or character strings containing lavaan model syntax. Can also accept a
#'   single list of models. Model names are automatically extracted from
#'   variable names or can be specified via named arguments.
#' @param data A data frame containing the variables for model fitting. Required
#'   when providing model specifications as character strings; ignored when
#'   providing pre-fitted lavaan objects.
#' @param model_names Character vector of custom model names. If \code{NULL}
#'   (default), names are extracted from function call arguments or variable names.
#' @param indices Character vector specifying which fit indices to include.
#'   Default is \code{c("chisq", "df", "cfi", "tli", "rmsea", "srmr")}.
#'   A χ²/df column is automatically added when possible.
#' @param changes Character vector specifying which indices to compute delta
#'   (change) values for. Default is \code{c("cfi", "rmsea")}.
#' @param ref Character string specifying the reference model for delta
#'   calculations. Options are:
#'   \itemize{
#'     \item \code{"first"} (default): All models compared to the first model
#'     \item \code{"last"}: Sequential comparison (each model vs. previous model)
#'   }
#' @param digits Integer specifying the number of decimal places for rounding.
#'   Default is 3.
#'
#' @return A data frame containing model comparison results with columns:
#'   \itemize{
#'     \item \code{Model}: Model names (extracted or provided)
#'     \item \code{Estimator}: Estimation method used (ML, MLR, etc.)
#'     \item Fit indices as specified in \code{indices} parameter
#'     \item \code{χ²/df}: Chi-square to degrees of freedom ratio
#'     \item Delta columns (e.g., \code{ΔCFI}, \code{ΔRMSEA}) showing changes
#'           relative to reference model
#'   }
#'
#'   The returned object includes these attributes:
#'   \itemize{
#'     \item \code{measure_names_used}: Named list showing which specific fit
#'           measures were used for each model (e.g., "cfi.robust" vs "cfi")
#'     \item \code{estimator_families}: Character vector indicating estimator
#'           family ("robust", "ml", or "unknown") for each model
#'     \item \code{original_indices}: The \code{indices} argument used
#'     \item \code{original_changes}: The \code{changes} argument used
#'     \item \code{reference_type}: The \code{ref} argument used
#'   }
#'
#' @details
#' The function automatically selects the most appropriate fit measures based on
#' the estimator used. For robust estimators (MLR, MLM, WLSMV, etc.), it
#' prioritizes robust or scaled versions of fit indices when available.
#'
#' **Delta Calculation Logic:**
#' \itemize{
#'   \item \code{ref = "first"}: Δ values show difference from first model
#'         (Model 1 gets NA, others show Model_i - Model_1)
#'   \item \code{ref = "last"}: Δ values show difference from immediately
#'         previous model (Model 1 gets NA, Model 2 shows Model_2 - Model_1,
#'         Model 3 shows Model_3 - Model_2, etc.)
#' }
#'
#' **Model Name Extraction:**
#' Names are determined in this priority order:
#' \enumerate{
#'   \item Explicit \code{model_names} argument
#'   \item Named arguments (e.g., \code{validity_table("My Model" = fit)})
#'   \item Variable names from function call (e.g., \code{fm_t1} becomes "fm_t1")
#'   \item Fallback to "Model i"
#' }
#'
#' @examples
#' \donttest{
#' library(lavaan)
#' data("HolzingerSwineford1939")
#'
#' # Define model specifications
#' model_3f <- '
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#' '
#'
#' model_3f_orth <- '
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#'   visual  ~~ 0*textual
#'   visual  ~~ 0*speed
#'   textual ~~ 0*speed
#' '
#'
#' model_1f <- '
#'   general =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
#' '
#'
#' # Example 1: Using pre-fitted lavaan objects
#' # Variable names (fit_3f, fit_1f, etc.) become model names
#' fit_3f <- cfa(model_3f, data = HolzingerSwineford1939)
#' fit_3f_orth <- cfa(model_3f_orth, data = HolzingerSwineford1939)
#' fit_1f <- cfa(model_1f, data = HolzingerSwineford1939)
#'
#' result1 <- validity_table(fit_3f, fit_3f_orth, fit_1f, digits = 3)
#' print(result1)
#'
#' # Example 2: Using model specifications (automatic fitting)
#' result2 <- validity_table(model_3f, model_3f_orth, model_1f,
#'                          data = HolzingerSwineford1939,
#'                          model_names = c("3F Correlated", "3F Orthogonal", "1F"),
#'                          digits = 3)
#' print(result2)
#'
#' # Example 3: Sequential comparison (each model vs. previous)
#' result3 <- validity_table(fit_3f, fit_3f_orth, fit_1f,
#'                          ref = "last",  # Sequential comparison
#'                          digits = 3)
#' print(result3)
#'
#' # Example 4: Using named arguments for custom model labels
#' result4 <- validity_table(
#'   "Three Factor" = fit_3f,
#'   "Orthogonal" = fit_3f_orth,
#'   "Single Factor" = fit_1f,
#'   digits = 3
#' )
#' print(result4)
#'
#' # Example 5: Custom fit indices and changes
#' result5 <- validity_table(fit_3f, fit_3f_orth, fit_1f,
#'                          indices = c("chisq", "df", "cfi", "tli", "rmsea", "srmr", "aic"),
#'                          changes = c("cfi", "rmsea", "aic"),
#'                          digits = 3)
#' print(result5)
#'
#' # Example 6: Measurement invariance testing
#' # (Sequential comparison is particularly useful here)
#' configural <- '
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#' '
#'
#' metric <- '
#'   visual  =~ a*x1 + b*x2 + c*x3
#'   textual =~ d*x4 + e*x5 + f*x6
#'   speed   =~ g*x7 + h*x8 + i*x9
#' '
#'
#' fit_config <- cfa(configural, data = HolzingerSwineford1939)
#' fit_metric <- cfa(metric, data = HolzingerSwineford1939)
#'
#' invariance_test <- validity_table(fit_config, fit_metric,
#'                                  model_names = c("Configural", "Metric"),
#'                                  ref = "last",  # Compare metric to configural
#'                                  digits = 3)
#' print(invariance_test)
#'
#' # Inspect which specific fit measures were used
#' attr(invariance_test, "measure_names_used")
#' }
#'
#' @seealso
#' \code{\link[lavaan]{cfa}}, \code{\link[lavaan]{sem}},
#' \code{\link[lavaan]{fitmeasures}}, \code{\link[lavaan]{lavInspect}}
#'
#' @export
validity_table <- function(...,
data = NULL,
model_names = NULL,
indices = c("chisq","df","cfi","tli","rmsea","srmr"),
changes = c("cfi","rmsea"),
ref = c("first", "last"),
digits = 3,
show_estimator_info = FALSE) {

  # Validate ref argument
  ref <- match.arg(ref)

  # Capture model arguments and their names - this is the key fix
  call_args <- match.call(expand.dots = FALSE)
  model_args <- call_args[["..."]]
  dots <- list(...)

  # Handle case where a single list of models is passed
  if (length(dots) == 1 && is.list(dots[[1]]) &&
      all(sapply(dots[[1]], function(x) inherits(x, "lavaan") || is.character(x)))) {
    inner_list <- dots[[1]]
    dots <- inner_list
    # Use names from the list if available
    if (!is.null(names(inner_list)) && any(nzchar(names(inner_list)))) {
      model_args <- as.list(names(inner_list))
    } else {
      model_args <- as.list(paste("Model", seq_along(inner_list)))
    }
  }

  if (length(dots) == 0) {
    stop("Provide at least one lavaan model object or model specification.")
  }

  # Determine if inputs are model specifications or fitted objects
  are_specs <- all(sapply(dots, is.character))
  are_fitted <- all(sapply(dots, inherits, "lavaan"))

  if (!are_specs && !are_fitted) {
    stop("All inputs must be either character strings (model specifications) or fitted lavaan objects.")
  }

  # Fit models if specifications are provided
  if (are_specs) {
    if (is.null(data)) {
      stop("When providing model specifications, 'data' argument is required.")
    }

    models <- vector("list", length(dots))
    for (i in seq_along(dots)) {
      spec <- dots[[i]]
      tryCatch({
        models[[i]] <- lavaan::cfa(model = spec, data = data)
      }, error = function(e) {
        stop(sprintf("Error fitting model %d: %s", i, e$message), call. = FALSE)
      })
    }
  } else {
    models <- dots
  }

  # Validate all models are lavaan objects
  if (!all(sapply(models, inherits, "lavaan"))) {
    stop("All inputs must result in valid lavaan model fits.")
  }

  # Extract model names - improved approach
  if (is.null(model_names)) {
    # Extract names from the actual call arguments
    extracted_names <- character(length(model_args))

    for (i in seq_along(model_args)) {
      arg <- model_args[[i]]
      if (is.symbol(arg)) {
        # Direct variable reference like fm_t1
        extracted_names[i] <- as.character(arg)
      } else if (is.character(arg)) {
        # String literal
        extracted_names[i] <- arg
      } else if (is.call(arg)) {
        # Complex expression - use deparse
        extracted_names[i] <- paste(deparse(arg), collapse = "")
      } else {
        extracted_names[i] <- paste("Model", i)
      }
    }

    # Handle named arguments
    arg_names <- names(model_args)
    if (!is.null(arg_names)) {
      named_positions <- which(nzchar(arg_names))
      extracted_names[named_positions] <- arg_names[named_positions]
    }

    model_names <- extracted_names

  } else if (length(model_names) != length(dots)) {
    stop("`model_names` must match the number of models.")
  }

  # Ensure we have chisq & df for χ²/df calculation
  all_idx <- unique(c(indices, "chisq", "df"))

  results <- vector("list", length(models))
  used_meas <- vector("list", length(models))
  est_info <- character(length(models))

  for (i in seq_along(models)) {
    mod <- models[[i]]
    est <- lavaan::lavInspect(mod, "options")$estimator
    est_info[i] <- est
    family <- get_estimator_family(mod)
    meas <- lavaan::fitmeasures(mod)

    # Collect requested measures
    row <- list()
    used_names <- character(length(all_idx) + 1)
    names(used_names) <- c(all_idx, "chisq_df")

    for (nm in all_idx) {
      res <- get_best_measure(meas, nm, family)
      row[[nm]] <- res$value
      used_names[nm] <- res$used
    }

    # Compute χ²/df
    if (!is.na(row$chisq) && !is.na(row$df) && row$df > 0) {
      row$chisq_df <- row$chisq / row$df
      used_names["chisq_df"] <- paste0(used_names["chisq"], "/", used_names["df"])
    } else {
      row$chisq_df <- NA_real_
      used_names["chisq_df"] <- NA_character_
    }

    results[[i]] <- as.data.frame(row)
    used_meas[[i]] <- used_names
  }

  # Build data.frame with fallback for bind_rows
  if (requireNamespace("dplyr", quietly = TRUE)) {
    fit_df <- dplyr::bind_rows(results)
  } else {
    fit_df <- do.call(rbind, results)
    rownames(fit_df) <- NULL
  }

  fit_df$Model <- model_names
  fit_df$Estimator <- est_info

  # Order columns: Model, Estimator, chisq, df, chisq_df, <others>
  final_raw <- indices
  if (all(c("chisq", "df") %in% indices)) {
    df_pos <- which(indices == "df")
    final_raw <- append(indices, values = "chisq_df", after = df_pos)
  } else {
    final_raw <- unique(c(indices, "chisq_df"))
  }
  fit_df <- fit_df[, c("Model", "Estimator", final_raw)]

  # Rename to publication-style column names
  new_names <- sapply(names(fit_df), function(x) {
    if (x %in% c("Model", "Estimator")) return(x)
    if (x == "chisq_df") return("χ²/df")
    format_column_name(x)
  }, USE.NAMES = FALSE)
  names(fit_df) <- new_names

  # Add Δ-columns for model comparisons based on reference type
  if (nrow(fit_df) > 1) {
    for (m in changes) {
      fm <- format_column_name(m)
      col <- paste0("Δ", fm)
      if (fm %in% names(fit_df)) {
        # Convert to numeric to handle potential character formatting
        numeric_vals <- suppressWarnings(as.numeric(fit_df[[fm]]))
        delta_values <- rep(NA_real_, length(numeric_vals))

        if (ref == "first") {
          # Compare all models to the first model
          base_value <- numeric_vals[1]
          if (!is.na(base_value)) {
            delta_values[-1] <- numeric_vals[-1] - base_value
            # First model gets NA (it's the reference)
            delta_values[1] <- NA
          }
        } else if (ref == "last") {
          # Sequential comparison: each model vs the immediately previous one
          for (i in 2:length(numeric_vals)) {
            if (!is.na(numeric_vals[i]) && !is.na(numeric_vals[i-1])) {
              delta_values[i] <- numeric_vals[i] - numeric_vals[i-1]
            }
          }
          # First model gets NA (no previous model to compare to)
          delta_values[1] <- NA
        }

        fit_df[[col]] <- delta_values
      }
    }
  }

  # Apply rounding to all numeric columns
  is_num <- sapply(fit_df, is.numeric)
  fit_df[is_num] <- lapply(fit_df[is_num], function(x) round(x, digits))

  # Format χ² and χ²/df with fixed decimal places for consistent display
  fmt_fixed <- function(x) {
    if (all(is.na(x))) return(x)
    sprintf(paste0("%.", digits, "f"), x)
  }

  if ("χ²" %in% names(fit_df)) {
    fit_df[["χ²"]] <- fmt_fixed(as.numeric(fit_df[["χ²"]]))
  }
  if ("χ²/df" %in% names(fit_df)) {
    fit_df[["χ²/df"]] <- fmt_fixed(as.numeric(fit_df[["χ²/df"]]))
  }

  # Keep df as integer for cleaner display
  if ("df" %in% names(fit_df)) {
    fit_df[["df"]] <- as.integer(round(as.numeric(fit_df[["df"]]), 0))
  }

  # Store metadata as attributes
  attr(fit_df, "measure_names_used") <- setNames(used_meas, model_names)
  attr(fit_df, "estimator_families") <- sapply(models, get_estimator_family)
  attr(fit_df, "original_indices") <- indices
  attr(fit_df, "original_changes") <- changes
  attr(fit_df, "reference_type") <- ref

  return(fit_df)
}
