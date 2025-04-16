#' Generate a Regression Coefficient Table with Fit Indices
#'
#' @description
#' This function generates a formatted table displaying regression coefficients and model fit indices
#' for one or multiple model objects. These model objects are typically created by the
#' \pkg{MpluAutomation} package, for example via the \code{readModels()} function or as lists returned
#' by the \code{extractModelSummaries()} function. The output table is formatted using
#' \code{knitr::kable} and enhanced with \code{kableExtra} functions.
#'
#' @param ... One or more model objects or nested lists of model objects from the
#' \pkg{MpluAutomation} package.
#' @param stat Character vector specifying which statistics to include. The default is
#' \code{c("est", "se", "p")}. Valid options include \code{"est"}, \code{"se"}, \code{"p"}, and \code{"z"}.
#' Setting \code{stat = "all"} will include all four statistics: \code{"est"}, \code{"se"}, \code{"z"}, and \code{"p"}.
#' @param indices Additional names of model fit indices to include.
#' @param digits_coeff Integer specifying the number of decimal places for coefficient values (default: 3).
#' @param digits_fit Integer specifying the number of decimal places for fit indices (default: 3).
#' @param standardized Logical; if \code{TRUE}, standardized coefficients (from \code{stdyx.standardized}) are used,
#' otherwise unstandardized coefficients are used (default: \code{FALSE}).
#' @param params Character vector indicating which parameter types to include (default: \code{c("regression")}).
#' @param model_names Optional character vector providing names for the models. If \code{NULL}, names are extracted
#' from the model's summary information (e.g., from \code{Filename}) or default names are assigned.
#' @param stars Logical; if \code{TRUE}, significance stars are appended to coefficient estimates based on p-value thresholds.
#' @param title Character string for the table caption (default: \code{"Table X. Unstandardized Regression Coefficients and Model Fit Indices."}).
#' @param notes Character string note to be added as a table footnote (default: \code{"Unstandardized coefficients are reported."}).
#' @param exclude Character vector of substrings; any rows containing these substrings in the right-hand side (RHS)
#' of the model equation will be excluded.
#' @param order_by Character value (\code{"lhs"} or \code{"rhs"}) indicating whether to order the rows by the left-hand
#' side or right-hand side of the equation (default: \code{"lhs"}).
#'
#' @details
#' The function is designed to work with model objects produced by the \pkg{MpluAutomation} package. It supports:
#'
#' \itemize{
#'   \item{\strong{Single model objects:} Directly pass a model created by \code{readModels()}.}
#'   \item{\strong{Complex model lists:} Pass a list of models, such as those returned by \code{extractModelSummaries()}.}
#'   \item{\strong{Combined input:} Mix individual model objects and nested lists of models. The function flattens the input
#'         for consistent processing.}
#' }
#'
#' Internally, helper functions such as \code{fmt_num()}, \code{classify_param()}, and
#' \code{create_formula_parts()} (which should be imported or defined in your environment) are used to format numbers,
#' classify parameters, and construct equation parts respectively.
#'
#' @return A formatted table (using \code{knitr::kable} and \code{kableExtra}) displaying the regression coefficients,
#' corresponding significance levels (with optional star markings), and model fit indices.
#'
#' @seealso \code{\link{readModels}}, \code{\link{extractModelSummaries}}, \code{knitr::kable}, \code{kableExtra::kable_classic}
#'
#' @examples
#' \dontrun{
#' # Example 1: Single model object from MpluAutomation
#' model_single <- readModels("model_file.out")
#' coeff_table(model_single)
#'
#' # Example 2: Complex model list returned by extractModelSummaries
#' model_list <- extractModelSummaries("models_directory")
#' coeff_table(model_list)
#'
#' # Example 3: Combined input with an individual model and a nested model list
#' model_extra <- readModels("another_model.out")
#' models_combined <- list(model_set = extractModelSummaries("other_models_folder"))
#' coeff_table(model_extra, models_combined)
#' }
#'
#' @import knitr
#' @importFrom dplyr bind_rows
#' @importFrom kableExtra kable_classic add_header_above row_spec footnote
#'
#' @export

coeff_table <- function(...,
                        stat = c("est", "se", "p"),
                        indices = character(0),
                        digits_coeff = 3,
                        digits_fit = 3,
                        standardized = FALSE,
                        params = c("regression"),
                        model_names = NULL,
                        stars = FALSE,
                        title = "Table X. Unstandardized Regression Coefficients and Model Fit Indices.",
                        notes = "Unstandardized coefficients are reported.",
                        exclude = character(0),      # Rows with these substrings in RHS will be filtered out.
                        order_by = c("lhs", "rhs")     # Choose ordering by left-hand side or right-hand side.
) {

  # If stat is set to "all", include all statistics: est, se, z and p.
  if (length(stat) == 1 && stat == "all") {
    stat <- c("est", "se", "z", "p")
  }

  ## --- Helper for star assignment ---
  add_stars <- function(est, p) {
    stars_chr <- ifelse(
      is.na(p), "",
      ifelse(p < 0.001, "***",
             ifelse(p < 0.01,  "**",
                    ifelse(p < 0.05,  "*", ""))))
    paste0(est, stars_chr)
  }

  ## --- Helper to recursively flatten a complex list of models ---
  flatten_models <- function(x, parent_name = NULL) {
    flat <- list()
    if (is.list(x)) {
      if ("summaries" %in% names(x)) {
        nm <- parent_name
        if (is.null(nm)) nm <- "model"
        result <- list(x)
        names(result) <- nm
        return(result)
      } else {
        for (i in seq_along(x)) {
          nm <- names(x)[i]
          if (is.null(nm) || nm == "") nm <- paste0("Model", i)
          flat <- c(flat, flatten_models(x[[i]], nm))
        }
      }
    }
    return(flat)
  }

  ## --- STEP 1: Flatten input objects into a list of models ---
  args <- list(...)
  modelList <- list()
  for (obj in args) {
    modelList <- c(modelList, flatten_models(obj))
  }
  if (length(modelList) == 0) {
    stop("No valid model objects were found in the input.")
  }

  ## --- Determine model names ---
  if (is.null(model_names)) {
    model_names <- sapply(modelList, function(m) {
      if (!is.null(m$summaries$Filename) && nzchar(m$summaries$Filename)) {
        return(m$summaries$Filename)
      } else {
        return(NA_character_)
      }
    })
    missing <- which(is.na(model_names) | model_names == "")
    if (length(missing) > 0) {
      model_names[missing] <- paste0("Model ", missing)
    }
  }
  if (length(model_names) != length(modelList)) {
    model_names <- paste0("Model ", seq_len(length(modelList)))
  }

  num_models <- length(modelList)
  model_tables <- list()
  fit_list <- list()

  for (i in seq_len(num_models)) {
    mod_flat <- modelList[[i]]
    if (standardized) {
      if (!is.null(mod_flat$parameters$stdyx.standardized)) {
        coeffs <- mod_flat$parameters$stdyx.standardized
      } else {
        stop(sprintf("Standardized results requested but not found for model %d. Processing stopped.", i))
      }
    } else {
      coeffs <- mod_flat$parameters$unstandardized
    }

    required_cols <- c("paramHeader", "param", "est", "se", "pval", "est_se")
    if (!all(required_cols %in% colnames(coeffs))) {
      stop("Coefficient data missing one or more required columns: ",
           paste(required_cols, collapse = ", "))
    }
    if (!"p" %in% colnames(coeffs) && "pval" %in% colnames(coeffs)) {
      names(coeffs)[names(coeffs) == "pval"] <- "p"
    }

    coeffs$ptype <- mapply(classify_param, coeffs$paramHeader, coeffs$param)
    coeffs <- coeffs[coeffs$ptype %in% params, ]

    formula_parts <- t(mapply(create_formula_parts, coeffs$paramHeader, coeffs$param))
    formula_parts <- as.data.frame(formula_parts, stringsAsFactors = FALSE)
    coeffs <- cbind(formula_parts, coeffs)

    ## --- Exclusion Filtering ---
    if (length(exclude) > 0) {
      # Create a regex pattern that matches any of the exclude strings.
      pattern <- paste(exclude, collapse = "|")
      coeffs <- coeffs[ !grepl(pattern, coeffs$rhs, ignore.case = TRUE), ]
    }

    p_orig <- as.numeric(coeffs$p)
    if ("est" %in% stat)  coeffs$est <- fmt_num(coeffs$est, digits_coeff)
    if ("se" %in% stat)   coeffs$se  <- fmt_num(coeffs$se, digits_coeff)
    if ("z" %in% stat)    coeffs$z   <- fmt_num(coeffs$est_se, digits_coeff)
    if ("p" %in% stat) {
      coeffs$p <- ifelse(p_orig < 0.001,
                         "<0.001",
                         fmt_num(p_orig, digits_coeff))
    }
    if (stars && "est" %in% stat) {
      coeffs$est <- add_stars(coeffs$est, p_orig)
    }

    for (s in stat) {
      original_name <- s
      new_name <- paste0("M", i, "_", s)
      names(coeffs)[names(coeffs) == original_name] <- new_name
    }

    expected_coef_cols <- c("lhs", "op", "rhs", paste0("M", i, "_", stat))
    if (nrow(coeffs) == 0) {
      coeffs <- data.frame(matrix(ncol = length(expected_coef_cols), nrow = 0),
                           stringsAsFactors = FALSE)
      colnames(coeffs) <- expected_coef_cols
    } else {
      coeffs <- coeffs[, expected_coef_cols, drop = FALSE]
    }
    model_tables[[i]] <- coeffs

    ## --- Process fit indices ---
    summ <- mod_flat$summaries
    default_indices <- c("ChiSqM_Value", "ChiSqM_DF", "χ²/df", "CFI", "TLI", "RMSEA_Estimate", "SRMR")
    all_indices <- unique(c(default_indices, indices))

    rename_map <- c(
      "ChiSqM_Value"   = "χ²",
      "ChiSqM_DF"      = "df",
      "RMSEA_Estimate" = "RMSEA"
    )
    idx_labels <- ifelse(all_indices %in% names(rename_map),
                         rename_map[all_indices],
                         all_indices)

    fit_vals <- sapply(all_indices, function(idx) {
      if (idx == "ChiSqM_DF") {
        if ("ChiSqM_DF" %in% names(summ))
          as.character(as.integer(summ[["ChiSqM_DF"]]))
        else NA_character_
      } else if (idx == "Parameters") {
        if ("Parameters" %in% names(summ))
          as.character(as.integer(summ[["Parameters"]]))
        else NA_character_
      } else if (idx == "χ²/df") {
        if ("ChiSqM_Value" %in% names(summ) && "ChiSqM_DF" %in% names(summ) && summ[["ChiSqM_DF"]] != 0)
          fmt_num(summ[["ChiSqM_Value"]] / summ[["ChiSqM_DF"]], digits_fit)
        else NA_character_
      } else {
        if (idx %in% names(summ)) fmt_num(summ[[idx]], digits_fit) else NA_character_
      }
    }, USE.NAMES = FALSE)

    fit_row <- data.frame(lhs = idx_labels,
                          op  = rep("", length(idx_labels)),
                          rhs = rep("", length(idx_labels)),
                          stringsAsFactors = FALSE)
    for (s in stat) {
      fit_row[[s]] <- if (s == stat[1]) fit_vals else NA_character_
    }
    for (s in stat) {
      new_name <- paste0("M", i, "_", s)
      names(fit_row)[names(fit_row) == s] <- new_name
    }
    fit_list[[i]] <- fit_row
  }  # end loop over models

  ## --- Merge coefficient tables ---
  merged_coeff <- Reduce(function(x, y) full_join(x, y, by = c("lhs", "op", "rhs")), model_tables)

  ## --- Ordering the merged coefficients ---
  # Define a fixed operator order.
  op_order <- c("=~" = 1, "~" = 2, "~ 1" = 3, "~~" = 4, ":=" = 5)
  merged_coeff$order <- op_order[merged_coeff$op]
  # Use the new argument to decide how to order.
  order_by <- match.arg(order_by)
  if(order_by == "lhs"){
    merged_coeff <- merged_coeff[order(merged_coeff$order, merged_coeff$lhs), ]
  } else {  # order_by == "rhs"
    merged_coeff <- merged_coeff[order(merged_coeff$order, merged_coeff$rhs), ]
  }
  merged_coeff$order <- NULL

  ## --- Merge fit indices rows ---
  merged_fit <- Reduce(function(x, y) full_join(x, y, by = c("lhs", "op", "rhs")), fit_list)

  final_table <- dplyr::bind_rows(merged_coeff, merged_fit)
  row.names(final_table) <- NULL
  final_table[is.na(final_table)] <- ""

  ## --- Force final_table to have all expected columns ---
  expected_cols <- c("lhs", "op", "rhs",
                     unlist(lapply(seq_len(num_models), function(i)
                       paste0("M", i, "_", stat))))
  for (col in expected_cols) {
    if (!(col %in% colnames(final_table))) {
      final_table[[col]] <- ""
    }
  }
  final_table <- final_table[, expected_cols, drop = FALSE]

  ## --- Set Uniform Column Names for Display ---
  if (standardized) {
    mapping <- c(est = "β", se = "SE", z = "z", p = "p")
  } else {
    mapping <- c(est = "b", se = "SE", z = "z", p = "p")
  }
  lower_stat <- sapply(stat, function(s) {
    if (s %in% names(mapping)) mapping[[s]] else s
  }, USE.NAMES = FALSE)
  model_stats_header <- unlist(lapply(seq_len(num_models), function(x) lower_stat))
  new_colnames <- c("LHS", "Operator", "RHS", model_stats_header)

  ## --- Build Header Groups for kableExtra ---
  header_groups <- c(" " = 3)
  for (i in seq_along(model_names)) {
    header_groups[ model_names[i] ] <- length(stat)
  }

  group_row_index <- nrow(merged_coeff)

  knitr::kable(final_table,
               caption = title,
               col.names = new_colnames) %>%
    kableExtra::kable_classic(html_font = "Cambria", full_width = FALSE, position = "left") %>%
    kableExtra::add_header_above(header_groups) %>%
    kableExtra::row_spec(group_row_index, extra_css = "border-bottom: 1px solid;") %>%
    kableExtra::footnote(general = notes, footnote_as_chunk = TRUE)
}
