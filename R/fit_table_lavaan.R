# Helper functions (internal)
# 1. Detect which family of fit‐indices to use, based on estimator
#' @noRd
get_estimator_family <- function(lavaan_model) {
  est <- toupper(lavInspect(lavaan_model, "options")$estimator)

  # Estimators that typically produce SCALED fit indices (e.g., for categorical data)
  scaled_estimators <- c("WLSM", "WLSMV", "DWLS", "ULS", "ULSM", "ULSMV", "MLM", "MLMV", "MLMVS")

  # Estimators that typically produce ROBUST standard errors/chi-square (e.g., MLR)
  robust_estimators <- c("MLR")

  if (est %in% scaled_estimators) {
    "scaled"
  } else if (est %in% robust_estimators) {
    "robust"
  } else {
    # Default for ML, FIML, etc.
    "standard"
  }
}

# 2. Pick the best available version of each index (scaled → robust → standard)
#' @noRd
get_best_measure <- function(all_meas, name, family) {
  if (family=="scaled") {
    candidates <- c(paste0(name,".scaled"),
                    paste0(name,".robust"),
                    name)
  } else if (family=="robust") {
    candidates <- c(paste0(name,".robust"),
                    paste0(name,".scaled"),
                    name)
  } else {
    candidates <- c(name,
                    paste0(name,".scaled"),
                    paste0(name,".robust"))
  }
  for (cand in candidates) {
    if (cand %in% names(all_meas) && !is.na(all_meas[cand])) {
      return(list(value=as.numeric(all_meas[cand]), used=cand))
    }
  }
  return(list(value=NA_real_, used=NA_character_))
}

# 3. Map raw measure names → publication‐style column names
#' @noRd
format_column_name <- function(x) {
  name_map <- c(
    chisq = "χ²",
    df    = "df",
    cfi   = "CFI",
    tli   = "TLI",
    rmsea = "RMSEA",
    srmr  = "SRMR",
    aic   = "AIC",
    bic   = "BIC",
    nnfi  = "NNFI",
    gfi   = "GFI",
    agfi  = "AGFI"
  )
  if (x %in% names(name_map)) {
    name_map[[x]]
  } else {
    toupper(x)
  }
}

#' Compare Multiple Lavaan Models with Publication-Ready Fit Statistics
#'
#' Creates a comprehensive comparison table for multiple lavaan model objects,
#' automatically selecting the most appropriate fit indices based on each model's
#' estimator family and providing publication-ready formatting with change statistics.
#'
#' @param ... Lavaan model objects (from \code{cfa()}, \code{sem()}, etc.) to compare,
#'   or a single list containing multiple lavaan model objects
#' @param model_names Character vector of names for the models. If \code{NULL},
#'   defaults to "Model 1", "Model 2", etc. Must match the number of models provided
#' @param indices Character vector specifying which fit indices to include.
#'   Available options: "chisq", "df", "cfi", "tli", "rmsea", "srmr", "aic", "bic",
#'   "nnfi", "gfi", "agfi". Default: \code{c("chisq","df","cfi","tli","rmsea","srmr")}
#' @param changes Character vector specifying which indices to calculate change scores
#'   for (Δ columns) relative to the first model. Default: \code{c("cfi","rmsea")}
#' @param digits Integer specifying number of decimal places for rounding and display.
#'   Default: 3. Chi-square and chi-square/df values use fixed decimal formatting
#' @param show_estimator_info Logical indicating whether to display detailed estimator
#'   information in console summary. Default: \code{FALSE}
#'
#' @return A data.frame with model comparison results including:
#' \itemize{
#'   \item Model names and estimator information
#'   \item Selected fit indices with appropriate scaling (standard/robust/scaled)
#'   \item Chi-square to degrees of freedom ratio (χ²/df) when both chisq and df are included
#'   \item Change scores (Δ) relative to the first model for specified indices
#' }
#'
#' The returned object includes attributes containing metadata about which specific
#' index versions were used for each model.
#'
#' @details
#' \strong{Estimator Family Detection:}
#' The function automatically detects estimator families and selects appropriate indices:
#' \itemize{
#'   \item \strong{Scaled family} (WLSM, WLSMV, DWLS, ULS, ULSM, ULSMV, MLM, MLMV):
#'         Prefers .scaled → .robust → standard versions
#'   \item \strong{Robust family} (MLR): Prefers .robust → .scaled → standard versions
#'   \item \strong{Standard family} (ML, etc.): Prefers standard → .scaled → .robust versions
#' }
#'
#' \strong{Special Formatting:}
#' Chi-square (χ²) and chi-square/df (χ²/df) values are formatted as character strings
#' with fixed decimal places to ensure consistent display regardless of magnitude.
#' Other indices remain numeric for further calculations.
#'
#' \strong{Change Interpretation:}
#' \itemize{
#'   \item For CFI/TLI/NNFI/GFI/AGFI: Positive changes indicate improvement
#'   \item For RMSEA/SRMR: Negative changes indicate improvement
#'   \item For AIC/BIC: Negative changes indicate improvement
#'   \item For chi-square: Negative changes indicate improvement (better fit)
#' }
#'
#' @examples
#' library(lavaan)
#'
#' # Define CFA models
#' model_1f <- 'general =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9'
#' model_3f <- '
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#' '
#'
#' # Fit models
#' fit_1f <- cfa(model_1f, data = HolzingerSwineford1939)
#' fit_3f <- cfa(model_3f, data = HolzingerSwineford1939)
#'
#' # Basic comparison
#' comparison <- fit_table_lavaan(
#'   fit_1f, fit_3f,
#'   model_names = c("One Factor", "Three Factors"),
#'   digits = 3
#' )
#' print(comparison)
#'
#' # Different change indice
#' comparison <- fit_table_lavaan(
#'   fit_1f, fit_3f,
#'   model_names = c("One Factor", "Three Factors"),
#'   digits = 3,
#'   changes = "srmr"
#' )
#' print(comparison)
#'
#' # Show detailed index selection
#' show_measure_details(comparison)
#'
#' @seealso \code{\link{show_measure_details}} for detailed index selection information
#' @importFrom lavaan fitmeasures lavInspect
#' @importFrom dplyr bind_rows
#' @export
fit_table_lavaan <- function(...,
                             model_names     = NULL,
                             indices         = c("chisq","df","cfi","tli","rmsea","srmr"),
                             changes         = c("cfi","rmsea"),
                             digits          = 3,
                             show_estimator_info = FALSE) {
  # capture models
  models <- list(...)
  if (length(models)==1 &&
      is.list(models[[1]]) &&
      all(sapply(models[[1]], inherits, "lavaan"))) {
    models <- models[[1]]
  }
  if (length(models)==0) {
    stop("Provide at least one lavaan model object.")
  }
  if (!all(sapply(models, inherits, "lavaan"))) {
    stop("All inputs must be lavaan model fits (cfa(), sem(), etc.).")
  }
  # model names
  if (is.null(model_names)) {
    model_names <- paste("Model", seq_along(models))
  } else if (length(model_names) != length(models)) {
    stop("`model_names` must match the number of models.")
  }
  # ensure we have chisq & df for χ²/df
  all_idx <- unique(c(indices, "chisq","df"))

  results       <- vector("list", length(models))
  used_meas     <- vector("list", length(models))
  est_info      <- character(length(models))

  for (i in seq_along(models)) {
    mod      <- models[[i]]
    est      <- lavInspect(mod, "options")$estimator
    est_info[i] <- est
    family   <- get_estimator_family(mod)
    meas     <- fitmeasures(mod)

    # collect requested measures
    row        <- list()
    used_names <- character(length(all_idx)+1)
    names(used_names) <- c(all_idx,"chisq_df")

    for (nm in all_idx) {
      res <- get_best_measure(meas, nm, family)
      row[[nm]]        <- res$value
      used_names[nm]   <- res$used
    }
    # compute χ²/df
    if (!is.na(row$chisq) && !is.na(row$df) && row$df>0) {
      row$chisq_df       <- row$chisq / row$df
      used_names["chisq_df"] <- paste0(
        used_names["chisq"],"/",used_names["df"]
      )
    } else {
      row$chisq_df       <- NA_real_
      used_names["chisq_df"] <- NA_character_
    }

    results[[i]]   <- as.data.frame(row)
    used_meas[[i]] <- used_names
  }

  # build data.frame
  fit_df <- bind_rows(results)
  fit_df$Model     <- model_names
  fit_df$Estimator <- est_info

  #  ─── order columns: Model, Estimator, chisq, df, chisq_df, <others> ─────
  final_raw <- indices
  if (all(c("chisq","df") %in% indices)) {
    df_pos <- which(indices=="df")
    # insert chisq_df AFTER df
    final_raw <- append(
      indices,
      values = "chisq_df",
      after = df_pos
    )
  } else {
    final_raw <- unique(c(indices, "chisq_df"))
  }
  fit_df <- fit_df[, c("Model","Estimator", final_raw)]

  # rename to publication‐style
  new_names <- sapply(names(fit_df), function(x) {
    if (x %in% c("Model","Estimator")) return(x)
    if (x=="chisq_df")       return("χ²/df")
    format_column_name(x)
  }, USE.NAMES=FALSE)
  names(fit_df) <- new_names

  #  ─── add Δ‐columns ─────────────────────────────────────────────────────
  if (nrow(fit_df)>1) {
    for (m in changes) {
      fm  <- format_column_name(m)
      col <- paste0("Δ", fm)
      if (fm %in% names(fit_df)) {
        base <- fit_df[[fm]][1]
        fit_df[[col]] <- c(NA, fit_df[[fm]][-1] - base)
      }
    }
  }

  #  ─── SINGLE COMPREHENSIVE ROUNDING AND DISPLAY FORMATTING ───────────────
  # Round all numeric columns first
  is_num <- sapply(fit_df, is.numeric)
  fit_df[is_num] <- lapply(fit_df[is_num], function(x) round(x, digits))

  # Format χ² and χ²/df with fixed decimal places for proper display
  # This ensures tibble printing shows exactly `digits` decimal places
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

  # store metadata
  attr(fit_df, "measure_names_used") <- setNames(used_meas, model_names)
  attr(fit_df, "estimator_families") <- sapply(models, get_estimator_family)
  attr(fit_df, "original_indices")   <- indices
  attr(fit_df, "original_changes")   <- changes

  #  ─── print concise summary ──────────────────────────────────────────────
  if (nrow(fit_df)>1) {
    cat("=== MODEL COMPARISON SUMMARY ===\n")
    cat(sprintf("(Changes relative to %s)\n\n", model_names[1]))
    if (show_estimator_info) {
      cats <- paste0(
        "- ", model_names, ": ", est_info,
        " (", attr(fit_df,"estimator_families"), ")\n"
      )
      cat("Estimator info:\n", paste(cats, collapse=""), "\n")
    }
    for (i in 2:nrow(fit_df)) {
      cat(model_names[i], ":\n")
      for (m in changes) {
        fm  <- format_column_name(m)
        dc  <- paste0("Δ", fm)
        if (dc %in% names(fit_df) && !is.na(fit_df[[dc]][i])) {
          val <- fit_df[[dc]][i]
          dir <- if (m=="cfi") {
            if (val>0) "improved" else "declined"
          } else if (m=="rmsea") {
            if (val<0) "improved" else "increased"
          } else {
            ""
          }
          cat(sprintf(paste0("  %s: %+.", digits, "f %s\n"), dc, val, dir))
        }
      }
      cat("\n")
    }
  }

  return(fit_df)
}

#' Display Detailed Information About Fit Index Selection
#'
#' Shows exactly which versions of fit indices (standard, robust, or scaled)
#' were automatically selected by \code{fit_table_lavaan} for each model based
#' on the estimator used. This transparency function helps users understand
#' the automatic index selection process.
#'
#' @param fit_table_result The result object returned by \code{fit_table_lavaan()}
#'
#' @return Invisibly returns \code{NULL}. The function is called for its side effect
#'   of printing detailed information about index usage to the console
#'
#' @details
#' For each model in the comparison, this function displays:
#' \itemize{
#'   \item The automatically detected estimator family (standard, robust, or scaled)
#'   \item The specific version of each fit index that was used (e.g., "cfi.scaled", "rmsea.robust")
#'   \item The mapping from requested generic index names to actual lavaan measure names
#' }
#'
#' This information is particularly useful when comparing models with different
#' estimators to ensure that appropriate index versions are being compared.
#'
#' @examples
#' \dontrun{
#' # After running fit_table_lavaan
#' comparison <- fit_table_lavaan(model1, model2)
#' show_measure_details(comparison)
#' }
#'
#' @seealso \code{\link{fit_table_lavaan}}
#' @export
show_measure_details <- function(fit_table_result) {
  info     <- attr(fit_table_result, "measure_names_used")
  families <- attr(fit_table_result, "estimator_families")
  if (is.null(info)) {
    cat("No measure‐use metadata found. Re-run with fit_table_lavaan(...).\n")
    return(invisible(NULL))
  }
  cat("=== DETAILED MEASURE INFORMATION ===\n\n")
  for (i in seq_along(info)) {
    model <- names(info)[i]
    fam   <- families[i]
    cat(sprintf("%s (family = %s):\n", model, fam))
    used  <- info[[i]]
    for (orig in names(used)) {
      act <- used[orig]
      if (is.na(act)) next
      disp_orig <- if (orig=="chisq_df") "χ²/df" else format_column_name(orig)
      cat(sprintf("  %s → %s\n", disp_orig, act))
    }
    cat("\n")
  }
}
