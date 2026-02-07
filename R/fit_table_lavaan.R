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

#' Compare Multiple CFA Models with Publication-Ready Fit Statistics
#'
#' Fits one or more CFA models using \code{lavaan::cfa()} and returns a comparison
#' table of common fit indices with publication-oriented formatting. Change statistics
#' (Δ columns) are computed relative to the first model.
#'
#' @param ... Lavaan model syntax strings (character). You can pass multiple model
#'   strings directly, or a single list of model strings.
#' @param data A data.frame used to fit all models via \code{lavaan::cfa()}.
#' @param fit_args A list of additional arguments passed to \code{lavaan::cfa()},
#'   such as \code{estimator}，\code{missing}，\code{std.lv}，\code{ordered}，\code{meanstructure}.
#'   Default is \code{list()}.
#' @param model_names Character vector of names for the models. If \code{NULL},
#'   defaults to "Model 1", "Model 2", etc. Must match the number of models provided.
#' @param indices Character vector specifying which fit indices to include.
#'   Common options include "chisq", "df", "cfi", "tli", "rmsea", "srmr", "aic", "bic".
#'   Default: \code{c("chisq","df","cfi","tli","rmsea","srmr")}.
#' @param changes Character vector specifying which indices to calculate change scores
#'   for (Δ columns) relative to the first model. Default: \code{c("cfi","rmsea")}.
#' @param digits Integer specifying the number of decimal places for rounding and display.
#'   Default: 3. \eqn{\chi^2} and \eqn{\chi^2/df} are formatted with fixed decimal places.
#' @param show_estimator_info Logical indicating whether to include an "Estimator" column
#'   in the output table. Default: \code{FALSE}.
#' @param show_summary Logical indicating whether to print a concise console summary of
#'   change statistics. Default: \code{FALSE}.
#'
#' @return A data.frame containing model comparison results, including:
#' \itemize{
#'   \item Model names (and optional estimator column)
#'   \item Requested fit indices
#'   \item \eqn{\chi^2/df} when both "chisq" and "df" are available
#'   \item Δ columns for indices listed in \code{changes}, relative to the first model
#' }
#'
#' The returned object includes attributes:
#' \itemize{
#'   \item \code{measure_names_used}: the exact fitmeasures names used for each model
#'   \item \code{estimator_families}: estimator family labels for each fitted model
#'   \item \code{fit_args}: the \code{fit_args} list actually used
#'   \item \code{fits}: the fitted lavaan objects (named by \code{model_names})
#' }
#'
#' @details
#' All input models are fitted inside the function using \code{lavaan::cfa()} with the
#' same \code{data} and \code{fit_args}. Fit indices are extracted via
#' \code{lavaan::fitmeasures()} and arranged into a single comparison table. When
#' multiple models are provided, Δ columns are computed as the difference between each
#' model and the first model for the indices listed in \code{changes}.
#'
#' @examples
#' library(lavaan)
#'
#' model_1f <- 'general =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9'
#' model_3f <- '
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#' '
#'
#' comparison <- fit_table_lavaan(
#'   model_1f, model_3f,
#'   data = HolzingerSwineford1939,
#'   model_names = c("One Factor", "Three Factors"),
#'   fit_args = list(estimator = "MLR", missing = "fiml"),
#'   digits = 3
#' )
#' print(comparison)
#'
#' @seealso \code{\link{show_measure_details}} for detailed index selection information.
#' @importFrom lavaan cfa fitmeasures lavInspect
#' @importFrom dplyr bind_rows
#' @export
fit_table_lavaan <- function(...,
                             data               = NULL,
                             fit_args           = list(),
                             model_names        = NULL,
                             indices            = c("chisq","df","cfi","tli","rmsea","srmr"),
                             changes            = c("cfi","rmsea"),
                             digits             = 3,
                             show_estimator_info = FALSE,
                             show_summary       = FALSE) {

  models <- list(...)

  if (length(models) == 1 && is.list(models[[1]]) && !is.character(models[[1]])) {
    models <- models[[1]]
  }

  if (length(models) == 0) stop("Provide at least one lavaan model syntax (character).")
  if (!all(sapply(models, is.character))) stop("All inputs must be lavaan model syntax strings (character).")
  if (is.null(data)) stop("`data` must be provided for fitting via lavaan::cfa().")

  if (is.null(model_names)) {
    model_names <- paste("Model", seq_along(models))
  } else if (length(model_names) != length(models)) {
    stop("`model_names` must match the number of models.")
  }

  # fit inside (always cfa)
  fits <- vector("list", length(models))
  for (i in seq_along(models)) {
    spec <- models[[i]]
    if (length(spec) > 1) spec <- paste(spec, collapse = "\n")

    args <- c(list(model = spec, data = data), fit_args)

    fits[[i]] <- tryCatch(
      do.call(lavaan::cfa, args),
      error = function(e) {
        stop(sprintf("Model fitting failed for %s: %s", model_names[i], e$message), call. = FALSE)
      }
    )
  }

  all_idx <- unique(c(indices, "chisq","df"))

  results   <- vector("list", length(fits))
  used_meas <- vector("list", length(fits))
  est_info  <- character(length(fits))

  for (i in seq_along(fits)) {
    mod         <- fits[[i]]
    est         <- lavInspect(mod, "options")$estimator
    est_info[i] <- est
    family      <- get_estimator_family(mod)
    meas        <- fitmeasures(mod)

    row        <- list()
    used_names <- character(length(all_idx) + 1)
    names(used_names) <- c(all_idx, "chisq_df")

    for (nm in all_idx) {
      res <- get_best_measure(meas, nm, family)
      row[[nm]]      <- res$value
      used_names[nm] <- res$used
    }

    if (!is.na(row$chisq) && !is.na(row$df) && row$df > 0) {
      row$chisq_df <- row$chisq / row$df
      used_names["chisq_df"] <- paste0(used_names["chisq"], "/", used_names["df"])
    } else {
      row$chisq_df <- NA_real_
      used_names["chisq_df"] <- NA_character_
    }

    results[[i]]   <- as.data.frame(row)
    used_meas[[i]] <- used_names
  }

  fit_df <- dplyr::bind_rows(results)
  fit_df$Model <- model_names
  if (show_estimator_info) fit_df$Estimator <- est_info

  final_raw <- indices
  if (all(c("chisq","df") %in% indices)) {
    df_pos <- which(indices == "df")
    final_raw <- append(indices, values = "chisq_df", after = df_pos)
  } else {
    final_raw <- unique(c(indices, "chisq_df"))
  }

  cols <- c("Model", final_raw)
  if (show_estimator_info) cols <- c("Model", "Estimator", final_raw)
  fit_df <- fit_df[, cols, drop = FALSE]

  new_names <- vapply(names(fit_df), function(x) {
    if (x %in% c("Model","Estimator")) return(x)
    if (x == "chisq_df") return("χ²/df")
    format_column_name(x)
  }, FUN.VALUE = character(1))
  names(fit_df) <- new_names

  if (nrow(fit_df) > 1) {
    for (m in changes) {
      fm  <- format_column_name(m)
      col <- paste0("Δ", fm)
      if (fm %in% names(fit_df)) {
        base <- fit_df[[fm]][1]
        fit_df[[col]] <- c(NA, fit_df[[fm]][-1] - base)
      }
    }
  }

  is_num <- vapply(fit_df, is.numeric, logical(1))
  fit_df[is_num] <- lapply(fit_df[is_num], function(x) round(x, digits))

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
  if ("df" %in% names(fit_df)) {
    fit_df[["df"]] <- as.integer(round(as.numeric(fit_df[["df"]]), 0))
  }

  attr(fit_df, "measure_names_used") <- stats::setNames(used_meas, model_names)
  attr(fit_df, "estimator_families") <- vapply(fits, get_estimator_family, character(1))
  attr(fit_df, "original_indices")   <- indices
  attr(fit_df, "original_changes")   <- changes
  attr(fit_df, "fit_args")           <- fit_args
  attr(fit_df, "fits")               <- stats::setNames(fits, model_names)

  if (show_summary && nrow(fit_df) > 1) {
    cat("=== MODEL COMPARISON SUMMARY ===\n")
    cat(sprintf("(Changes relative to %s)\n\n", model_names[1]))

    if (show_estimator_info) {
      cats <- paste0(
        "- ", model_names, ": ", est_info,
        " (", attr(fit_df, "estimator_families"), ")\n"
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
          dir <- if (m == "cfi") {
            if (val > 0) "improved" else "declined"
          } else if (m == "rmsea") {
            if (val < 0) "improved" else "increased"
          } else ""
          cat(sprintf(paste0("  %s: %+.", digits, "f %s\n"), dc, val, dir))
        }
      }
      cat("\n")
    }
  }

  fit_df
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
