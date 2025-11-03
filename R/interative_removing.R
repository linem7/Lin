#' Internal: compute reliability for (potential) second-order factor models
#'
#' Computes reliability indices for constructs specified in `var_list` by
#' fitting CFA models using lavaan. For second-order constructs, it reports
#' per first-order dimension reliability (omega-like, based on loadings)
#' and an overall item-level Cronbach's alpha by pooling all items across
#' first-order factors. For single-factor constructs, it reports an omega-like
#' reliability based on standardized loadings.
#'
#' Thresholds used downstream: first-order dimensions >= 0.70; pooled
#' second-order alpha >= 0.80.
#'
#' @param data data.frame of observed variables.
#' @param var_list named list of lavaan measurement model syntax (one per construct).
#' @param ... ignored (placeholder to be flexible with upstream callers).
#' @return A named list per construct with fields: `type`, `pass`,
#'   `alpha_overall`, `n_dimensions` (for second-order), and `dimensions`.
#' @keywords internal
calculate_reliability_alpha <- function(data, var_list, ...) {
  results <- list()

  for (construct_name in names(var_list)) {
    model_syntax <- var_list[[construct_name]]

    # Fit the model
    fit <- tryCatch({
      suppressWarnings(cfa(model_syntax, data = data, std.lv = TRUE))
    }, error = function(e) {
      warning("Fit ", construct_name, " failed: ", e$message)
      return(NULL)
    })

    if (is.null(fit)) {
      results[[construct_name]] <- list(
        type = "Unknown",
        pass = FALSE,
        alpha_overall = NA,
        dimensions = NULL,
        error = "Model fitting failed"
      )
      next
    }

    # Extract parameter estimates
    pe <- parameterEstimates(fit, standardized = TRUE)

    # Identify model structure
    # Get all latent variables
    all_latent <- unique(c(pe$lhs[pe$op == "=~"], pe$rhs[pe$op == "=~"]))

    # Get all observed variables (present in data)
    all_observed <- colnames(data)

    # Distinguish latent vs observed
    latent_vars <- all_latent[!(all_latent %in% all_observed)]

    # Latent->latent relations (second-order)
    higher_order_relations <- pe[pe$op == "=~" &
                                   pe$lhs %in% latent_vars &
                                   pe$rhs %in% latent_vars, ]

    is_hierarchical <- nrow(higher_order_relations) > 0

    if (is_hierarchical) {
      # === Second-order factor model ===

      # Name of higher-order factor
      higher_factor <- unique(higher_order_relations$lhs)[1]

      # First-order factor names
      first_order_factors <- unique(higher_order_relations$rhs)

      # Compute reliability per first-order dimension
      dim_results <- list()

      for (first_factor in first_order_factors) {
        # Items loaded by this first-order factor
        item_relations <- pe[pe$op == "=~" &
                               pe$lhs == first_factor &
                               pe$rhs %in% all_observed, ]

        if (nrow(item_relations) > 0) {
          item_loadings <- item_relations$std.all

          if (length(item_loadings) > 1) {
            sum_loadings <- sum(item_loadings)
            sum_loadings_sq <- sum(item_loadings^2)
            error_var <- length(item_loadings) - sum_loadings_sq
            omega_dim <- sum_loadings^2 / (sum_loadings^2 + error_var)

            dim_results[[first_factor]] <- list(
              alpha = omega_dim,
              n_items = length(item_loadings)
            )
          } else {
            # Single-item dimension
            dim_results[[first_factor]] <- list(
              alpha = NA,
              n_items = length(item_loadings)
            )
          }
        } else {
          cat("  [Warning] ", first_factor, " has no observed items found\n")
          # Keep placeholder even if no observed items, for reporting
          dim_results[[first_factor]] <- list(
            alpha = NA,
            n_items = 0
          )
        }
      }

      # Compute item-level Cronbach's alpha for higher-order factor (combine all items)
      all_items <- character(0)
      for (first_factor in first_order_factors) {
        item_relations <- pe[pe$op == "=~" & pe$lhs == first_factor & pe$rhs %in% all_observed, ]
        if (nrow(item_relations) > 0) {
          all_items <- c(all_items, item_relations$rhs)
        }
      }
      all_items <- unique(all_items)

      alpha_total <- tryCatch({
        if (length(all_items) >= 2) {
          S <- stats::cov(data[, all_items, drop = FALSE], use = "pairwise.complete.obs")
          p <- ncol(S)
          if (p >= 2 && is.finite(sum(S))) {
            p/(p - 1) * (1 - sum(diag(S)) / sum(S))
          } else {
            NA_real_
          }
        } else {
          NA_real_
        }
      }, error = function(e) NA_real_)

      # Pass/fail check
      valid_dim_alphas <- sapply(dim_results, function(x) {
        if (!is.null(x) && !is.na(x$alpha)) x$alpha else -1
      })
      all_dim_pass <- all(valid_dim_alphas >= 0.7)
      overall_pass <- !is.na(alpha_total) && alpha_total >= 0.8

      results[[construct_name]] <- list(
        type = "Second-order factor",
        pass = all_dim_pass && overall_pass,
        alpha_overall = alpha_total,
        n_dimensions = length(first_order_factors),
        dimensions = dim_results
      )

    } else {
      # === Single (first-order) factor model ===

      # Items for the (single) factor
      factor_name <- latent_vars[1]  # assume only one latent variable

      item_relations <- pe[pe$op == "=~" &
                             pe$lhs == factor_name &
                             pe$rhs %in% all_observed, ]

      if (nrow(item_relations) > 0) {
        loadings <- item_relations$std.all

        if (length(loadings) > 1) {
          sum_loadings <- sum(loadings)
          sum_loadings_sq <- sum(loadings^2)
          error_var <- length(loadings) - sum_loadings_sq
          omega <- sum_loadings^2 / (sum_loadings^2 + error_var)

          results[[construct_name]] <- list(
            type = "Single-dimension construct",
            pass = omega >= 0.7,
            alpha_overall = omega,
            n_items = length(loadings),
            dimensions = NULL
          )
        } else {
          results[[construct_name]] <- list(
            type = "Single-dimension construct",
            pass = FALSE,
            alpha_overall = NA,
            n_items = length(loadings),
            dimensions = NULL
          )
        }
      } else {
        results[[construct_name]] <- list(
          type = "Single-dimension construct",
          pass = FALSE,
          alpha_overall = NA,
          n_items = 0,
          dimensions = NULL
        )
      }
    }
  }

  return(results)
}

#' Internal: pretty-print reliability report
#'
#' Prints a human-readable reliability report for the output of
#' `calculate_reliability_alpha()`.
#'
#' @param rel_results list returned by `calculate_reliability_alpha()`.
#' @param criteria_alpha numeric, threshold for first-order dimensions (default 0.70).
#' @param criteria_alpha_overall2 numeric, threshold for pooled second-order alpha (default 0.80).
#' @keywords internal
print_reliability_report <- function(rel_results, criteria_alpha = 0.7, criteria_alpha_overall2 = 0.8) {
  cat("\n========================================\n")
  cat("Reliability Assessment Report\n")
  cat("========================================\n\n")

  for (construct_name in names(rel_results)) {
    result <- rel_results[[construct_name]]

    cat("[", construct_name, "]\n", sep = "")
    cat("  Type: ", result$type, "\n", sep = "")

    if (result$type == "Second-order factor") {
      # Second-order factor report
      cat("  Number of dimensions: ", result$n_dimensions, "\n\n", sep = "")

      if (!is.null(result$dimensions) && length(result$dimensions) > 0) {
        cat("  First-order dimensions:\n")
        dims <- result$dimensions
        nm <- names(dims)
        for (i in seq_along(dims)) {
          dim_info <- dims[[i]]
          dim_name <- if (!is.null(nm) && length(nm) >= i && !is.na(nm[i]) && nzchar(nm[i])) nm[i] else paste0("Dimension ", i)
          alpha_val <- dim_info$alpha
          n_items <- dim_info$n_items

          if (!is.na(alpha_val)) {
            status <- ifelse(alpha_val >= criteria_alpha, "✓", "✗")
            cat(sprintf("    %s: α=%.3f (n=%d) %s\n",
                        dim_name, alpha_val, n_items, status))
          } else {
            cat(sprintf("    %s: α=NA (n=%d) ✗ (not computable)\n",
                        dim_name, n_items))
          }
        }
      } else {
        cat("  First-order dimensions: (no data)\n")
      }

      cat("\n  Second-order factor:\n")
      if (!is.na(result$alpha_overall)) {
        status <- ifelse(result$alpha_overall >= criteria_alpha_overall2, "✓", "✗")
        cat(sprintf("    %s: α=%.3f (based on %d dimensions) %s\n",
                    construct_name, result$alpha_overall, result$n_dimensions, status))
      } else {
        cat(sprintf("    %s: α=NA (based on %d dimensions) ✗ (not computable)\n",
                    construct_name, result$n_dimensions))
      }

      overall_status <- ifelse(result$pass, "✓ All pass", "✗ Not all pass")
      cat("\n  Overall: ", overall_status, "\n\n", sep = "")

    } else if (result$type == "Single-dimension construct") {
      # Single-dimension report
      cat("  Number of items: ", result$n_items, "\n", sep = "")

      if (!is.na(result$alpha_overall)) {
        status <- ifelse(result$alpha_overall >= criteria_alpha, "✓", "✗")
        cat(sprintf("  Cronbach's α: %.3f %s\n\n", result$alpha_overall, status))
      } else {
        cat("  Cronbach's α: NA ✗ (not computable)\n\n")
      }

    } else {
      # Error case
      cat("  Status: evaluation failed\n")
      if (!is.null(result$error)) {
        cat("  Error: ", result$error, "\n", sep = "")
      }
      cat("\n")
    }
  }
}

#' Internal: unified evaluation for reliability, paths, and fit
#'
#' Fits models and evaluates three criteria: reliability, key path significance,
#' and overall fit indices (CFI/RMSEA). Used by `interative_removing()`.
#'
#' @param data data.frame of observed variables.
#' @param model_configs named list; each element has `model` (lavaan syntax)
#'   and optional `key_paths` (character vector of labeled paths to test).
#' @param var_list optional named list of measurement models passed to
#'   `calculate_reliability_alpha()`.
#' @param criteria list of thresholds/flags for `reliability`, `paths`, and `fit`.
#' @return list with sublists `reliability`, `paths`, and `fit`, each having
#'   `pass` (logical) and `details`.
#' @keywords internal
check_all_criteria <- function(data, model_configs, var_list = NULL, criteria) {
  results <- list(
    reliability = list(pass = TRUE, details = NULL),
    paths = list(pass = TRUE, details = NULL),
    fit = list(pass = TRUE, details = NULL)
  )

  # === 1. Reliability ===
  if (!is.null(var_list) && criteria$reliability$check) {
    cat("\n[Reliability check]\n")

    rel_results <- tryCatch({
      calculate_reliability_alpha(data, var_list)
    }, error = function(e) {
      cat("✗ Reliability analysis failed:", e$message, "\n")
      return(NULL)
    }, warning = function(w) {
      suppressWarnings(calculate_reliability_alpha(data, var_list))
    })

    if (!is.null(rel_results) && length(rel_results) > 0) {
      # Print formatted report (dimension threshold uses criteria$reliability$min_alpha; second-order uses 0.8)
      print_reliability_report(rel_results, criteria$reliability$min_alpha, 0.8)

      all_pass <- all(sapply(rel_results, function(x) x$pass))

      results$reliability <- list(
        pass = all_pass,
        details = rel_results
      )
    } else {
      cat("✗ Reliability analysis completely failed\n")
      results$reliability <- list(pass = FALSE, details = NULL)
    }
  }

  # === 2. Path significance ===
  if (criteria$paths$check) {
    cat("\n[Path significance check]\n")

    all_paths_pass <- TRUE
    path_details <- list()

    for (config_name in names(model_configs)) {
      config <- model_configs[[config_name]]

      fit <- tryCatch({
        cfa(config$model, data = data, std.lv = TRUE)
      }, error = function(e) {
        cat("✗ Model", config_name, "fit failed:", e$message, "\n")
        return(NULL)
      }, warning = function(w) {
        tryCatch(cfa(config$model, data = data, std.lv = TRUE),
                 error = function(e) NULL)
      })

      if (is.null(fit)) {
        all_paths_pass <- FALSE
        next
      }

      converged <- tryCatch({
        lavInspect(fit, "converged")
      }, error = function(e) FALSE)

      if (!converged) {
        cat("✗ Model", config_name, "did not converge\n")
        all_paths_pass <- FALSE
        next
      }

      pe <- tryCatch({
        parameterEstimates(fit)
      }, error = function(e) {
        cat("✗ Failed to extract parameter estimates\n")
        return(NULL)
      })

      if (is.null(pe)) {
        all_paths_pass <- FALSE
        next
      }

      if (!is.null(config$key_paths)) {
        for (path_label in config$key_paths) {
          path_row <- pe[pe$label == path_label, ]

          if (nrow(path_row) == 0) {
            cat("⚠️ Path", path_label, "not found\n")
            all_paths_pass <- FALSE
            next
          }

          if (is.na(path_row$pvalue) || is.null(path_row$pvalue)) {
            cat(sprintf("  ✗ %s: β = %.3f, p = NA (not computable)\n",
                        path_label, path_row$est))
            all_paths_pass <- FALSE

            path_details[[path_label]] <- list(
              estimate = path_row$est,
              pvalue = NA,
              significant = FALSE
            )
            next
          }

          is_sig <- path_row$pvalue < 0.05
          status <- ifelse(is_sig, "✓", "✗")

          cat(sprintf("  %s %s: β = %.3f, p = %.3f\n",
                      status, path_label, path_row$est, path_row$pvalue))

          if (!is_sig) all_paths_pass <- FALSE

          path_details[[path_label]] <- list(
            estimate = path_row$est,
            pvalue = path_row$pvalue,
            significant = is_sig
          )
        }
      }
    }

    results$paths <- list(
      pass = all_paths_pass,
      details = path_details
    )
  }

  # === 3. Fit indices ===
  if (criteria$fit$check) {
    cat("\n[Fit indices check]\n")

    all_fit_pass <- TRUE
    fit_details <- list()

    for (config_name in names(model_configs)) {
      config <- model_configs[[config_name]]

      fit <- tryCatch({
        cfa(config$model, data = data, std.lv = TRUE)
      }, error = function(e) NULL,
      warning = function(w) {
        tryCatch(cfa(config$model, data = data, std.lv = TRUE),
                 error = function(e) NULL)
      })

      if (is.null(fit)) {
        all_fit_pass <- FALSE
        next
      }

      converged <- tryCatch(lavInspect(fit, "converged"), error = function(e) FALSE)

      if (!converged) {
        all_fit_pass <- FALSE
        next
      }

      fit_measures <- tryCatch({
        fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr"))
      }, error = function(e) {
        cat("✗ Failed to extract fit indices\n")
        return(NULL)
      })

      if (is.null(fit_measures)) {
        all_fit_pass <- FALSE
        next
      }

      cfi_pass <- fit_measures["cfi"] >= criteria$fit$min_cfi
      rmsea_pass <- fit_measures["rmsea"] <= criteria$fit$max_rmsea

      cat(sprintf("  %s CFI = %.3f (criterion: ≥ %.2f)\n",
                  ifelse(cfi_pass, "✓", "✗"),
                  fit_measures["cfi"], criteria$fit$min_cfi))

      cat(sprintf("  %s RMSEA = %.3f (criterion: ≤ %.2f)\n",
                  ifelse(rmsea_pass, "✓", "✗"),
                  fit_measures["rmsea"], criteria$fit$max_rmsea))

      if (!cfi_pass || !rmsea_pass) all_fit_pass <- FALSE

      fit_details[[config_name]] <- as.list(fit_measures)
    }

    results$fit <- list(
      pass = all_fit_pass,
      details = fit_details
    )
  }

  return(results)
}


#' Internal: identify problematic cases for removal
#'
#' Strategy: prioritize reliability when requested, using a greedy
#' case-deleted Cronbach's alpha improvement (Method 1). When reliability
#' does not apply or Method 1 finds no positive-gain cases, fall back to a
#' generic z-score-based selection (largest absolute z-score sums).
#'
#' @param data data.frame with a temporary `..row_id` column for mapping.
#' @param model_configs named list of model configurations.
#' @param var_list named list of measurement models; required for reliability.
#' @param criteria list of thresholds/flags.
#' @param n_to_remove integer, desired number of cases to remove.
#' @param focus one of `"auto"`, `"reliability"`, `"paths"`, `"fit"`.
#' @return integer vector of row indices (relative to current `data`).
#' @keywords internal
identify_bad_cases <- function(data, model_configs, var_list, criteria, n_to_remove, focus = "auto") {
  cat("\n[Identify problematic cases]\n")

  # Print focus (no longer using influence_stat)
  if (identical(focus, "reliability")) {
    cat("Focus: Reliability (Method 1: remove cases to increase α)\n")
  } else if (identical(focus, "paths")) {
    cat("Focus: Path significance (use generic fallback selection)\n")
  } else if (identical(focus, "fit")) {
    cat("Focus: Model fit (use generic fallback selection)\n")
  } else {
    cat("Focus: Auto\n")
  }

  # Helper: compute Cronbach's alpha from covariance matrix
  alpha_from_cov <- function(S) {
    if (is.null(S) || any(!is.finite(S))) return(NA_real_)
    p <- ncol(S)
    if (is.null(p) || p < 2) return(NA_real_)
    denom <- sum(S)
    if (!is.finite(denom) || denom <= 0) return(NA_real_)
    p/(p - 1) * (1 - sum(diag(S)) / denom)
  }

  # Fallback: pick top cases by sum of absolute Z-scores across all variables
  fallback_pick_by_zsum <- function(df, k) {
    obs_names <- setdiff(colnames(df), "..row_id")
    X <- as.data.frame(df[, obs_names, drop = FALSE])
    suppressWarnings({
      Z <- scale(X)
    })
    # Row sum (absolute), NA treated as 0
    row_score <- rowSums(abs(Z), na.rm = TRUE)
    ord <- order(row_score, decreasing = TRUE)
    head(ord, min(k, length(ord)))
  }

  # Method 1: remove cases to improve alpha (used for reliability focus or auto first)
  try_method1 <- function(df, k) {
    if (is.null(var_list) || !isTRUE(criteria$reliability$check)) return(integer(0))

    obs_all <- setdiff(colnames(df), "..row_id")

    # Collect items per construct and current alpha
    cand <- list()
    for (construct_name in names(var_list)) {
      model_syntax <- var_list[[construct_name]]
      fit <- tryCatch({
        suppressWarnings(cfa(model_syntax, data = df, std.lv = TRUE))
      }, error = function(e) NULL)

      if (is.null(fit)) next
      converged <- tryCatch(lavInspect(fit, "converged"), error = function(e) FALSE)
      if (!converged) next

      pe <- tryCatch(parameterEstimates(fit, standardized = TRUE), error = function(e) NULL)
      if (is.null(pe)) next

      # Collect items
      is_loading <- pe$op == "=~"
      lhs <- pe$lhs[is_loading]
      rhs <- pe$rhs[is_loading]

      # Identify latent vs observed variables
      all_observed <- obs_all
      latent_lhs <- lhs[!(lhs %in% all_observed)]
      latent_rhs <- rhs[!(rhs %in% all_observed)]
      # latent->latent (second-order)
      is_hier <- any((lhs %in% latent_lhs) & (rhs %in% latent_rhs))

      items <- character(0)
      if (is_hier) {
        first_factors <- unique(rhs[(lhs %in% latent_lhs) & (rhs %in% latent_rhs)])
        for (f1 in first_factors) {
          it <- rhs[(lhs == f1) & (rhs %in% all_observed)]
          if (length(it) > 0) items <- c(items, it)
        }
      } else {
        # direct latent->observed items
        items <- rhs[(lhs %in% latent_lhs) & (rhs %in% all_observed)]
      }

      items <- unique(intersect(items, obs_all))
      if (length(items) < 2) next

      S0 <- tryCatch(stats::cov(df[, items, drop = FALSE], use = "pairwise.complete.obs"), error = function(e) NULL)
      a0 <- alpha_from_cov(S0)
      if (is.na(a0)) next

      cand[[construct_name]] <- list(items = items, alpha0 = a0, is_hier = is_hier)
    }

    if (length(cand) == 0) return(integer(0))

    # Choose target construct: prefer second-order with α<0.8; then single with α<0.7; else lowest α
    pick_target <- function(cand) {
      nm <- names(cand)
      df <- data.frame(
        name = nm,
        alpha0 = vapply(cand, function(x) x$alpha0, numeric(1)),
        is_hier = vapply(cand, function(x) x$is_hier, logical(1)),
        stringsAsFactors = FALSE
      )
      idx <- which(df$is_hier & df$alpha0 < 0.8)
      if (length(idx) == 0) idx <- which(!df$is_hier & df$alpha0 < 0.7)
      if (length(idx) == 0) idx <- which.min(df$alpha0)
      df$name[idx[1]]
    }

    target_name <- pick_target(cand)
    target <- cand[[target_name]]
    items <- target$items
    a0 <- target$alpha0

    cat(sprintf("Method 1 target: %s (items=%d, current α=%.3f)\n", target_name, length(items), a0))

    # Compute per-case deleted Δα
    n <- nrow(df)
    deltas <- rep(NA_real_, n)
    X <- as.data.frame(df[, items, drop = FALSE])
    for (i in seq_len(n)) {
      Xi <- X[-i, , drop = FALSE]
      S <- tryCatch(stats::cov(Xi, use = "pairwise.complete.obs"), error = function(e) NULL)
      ai <- alpha_from_cov(S)
      deltas[i] <- if (!is.na(ai)) ai - a0 else NA_real_
    }

    # 选择 Δα 最大的样本（优先正增益），不足则返回空让回退处理
    ok <- which(is.finite(deltas) & deltas > 0)
    if (length(ok) == 0) return(integer(0))
    ord <- ok[order(deltas[ok], decreasing = TRUE)]
    head(ord, min(k, length(ord)))
  }

  picks <- integer(0)

  # Decision: prefer Method 1 (reliability); otherwise fallback
  if (identical(focus, "reliability") || identical(focus, "auto")) {
    picks <- try_method1(data, n_to_remove)
    if (length(picks) > 0) {
      cat("  ✓ Method 1: selected ", length(picks), " case(s)\n", sep = "")
      return(picks)
    } else {
      cat("  ↪ No positive-gain cases found by Method 1; switching to fallback\n")
    }
  }

  # For paths/fit focus or when Method 1 fails, use fallback
  picks <- fallback_pick_by_zsum(data, n_to_remove)
  cat("  ✓ Fallback: selected ", length(picks), " case(s) by Z-score sum\n", sep = "")
  return(picks)
}

#' Removal-only iterative cleaning to meet reliability/paths/fit
#'
#' Iteratively removes problematic cases from a dataset until all criteria
#' are satisfied: (1) reliability (first-order dimensions >= 0.70 and pooled
#' second-order alpha >= 0.80), (2) key path significance (p < 0.05 for
#' labeled paths), and (3) fit indices (CFI >= 0.90; RMSEA <= 0.08).
#'
#' Focus priority per round: Fit > Paths > Reliability. When focusing on
#' reliability (or auto mode), cases are selected via a greedy, case-deleted
#' Cronbach's alpha improvement (Method 1). Otherwise a generic fallback
#' selects cases by the sum of absolute z-scores across variables.
#'
#' This function mutates a temporary `..row_id` to map removals back to
#' original rows; it is removed before returning.
#'
#' @param data data.frame of observed variables.
#' @param model_configs named list with elements containing `model` (lavaan syntax)
#'   and optional `key_paths` (character vector of labeled paths to test).
#' @param var_list optional named list of measurement models (lavaan syntax) keyed by construct name.
#' @param initial_remove_pct numeric in (0,1], initial removal proportion per round (default 0.05).
#' @param pct_increment numeric in (0,1], increment added to removal proportion each round (default 0.05).
#' @param max_iterations integer, maximum number of iterations (default 10).
#' @param max_total_remove_pct numeric in (0,1], cap on cumulative removals relative to initial n (default 0.5).
#' @param seed optional integer seed for reproducibility.
#'
#' @return A list with elements:
#'   - `cleaned_data`: the cleaned data.frame
#'   - `iterations`: number of iterations performed
#'   - `initial_n`: initial sample size
#'   - `final_n`: final sample size
#'   - `total_removed`: number of removed rows
#'   - `removed_row_ids`: original row ids removed
#'   - `removal_proportion`: last round removal proportion
#'   - `final_results`: evaluation results from the last round
#'
#' @section Criteria:
#' - Reliability: first-order dimensions >= 0.70; second-order pooled alpha >= 0.80.
#' - Paths: labeled paths significant at p < 0.05.
#' - Fit: CFI >= 0.90; RMSEA <= 0.08.
#'
#' @section Notes:
#' - Requires lavaan. Do not call `library(lavaan)` inside package code; instead
#'   use NAMESPACE imports. See `@importFrom` below.
#' - The greedy reliability improvement can be compute-intensive for large n;
#'   consider pre-filtering candidates for performance if needed.
#'
#' @seealso calculate_reliability_alpha, check_all_criteria, identify_bad_cases
#'
#' @examples
#'   # Minimal example: simulate data and run cleaning
#'   set.seed(123)
#'   pop_model <- '
#'     X =~ 0.80*f1 + 0.80*f2 + 0.80*f3
#'     f1 =~ 0.8*X1 + 0.82*X2 + 0.79*X3
#'     f2 =~ 0.81*X4 + 0.83*X5 + 0.78*X6
#'     f3 =~ 0.8*X7 + 0.82*X8 + 0.79*X9
#'     M =~ 0.8*M1 + 0.82*M2 + 0.81*M3 + 0.79*M4 + 0.83*M5
#'     Y =~ 0.8*Y1 + 0.81*Y2 + 0.82*Y3 + 0.79*Y4 + 0.8*Y5
#'     M ~ 0.5*X
#'     Y ~ 0.6*M
#'     X ~~ 1*X; f1 ~~ 1*f1; f2 ~~ 1*f2; f3 ~~ 1*f3; M ~~ 1*M; Y ~~ 1*Y
#'   '
#'   dat <- lavaan::simulateData(model = pop_model, sample.nobs = 500)
#'   obs <- c(paste0("X", 1:9), paste0("M", 1:5), paste0("Y", 1:5))
#'   dat <- as.data.frame(dat)[, obs, drop = FALSE]
#'
#'   # Inject 50 outliers using rnorm with inflated SD
#'   set.seed(999)
#'   n_out <- 50
#'   outliers <- as.data.frame(
#'     sapply(dat, function(x) rnorm(n_out, mean = mean(x), sd = sd(x) * 3))
#'   )
#'   dat <- rbind(dat, outliers)
#'
#'   var_list <- list(
#'     "Construct X" = 'X =~ f1 + f2 + f3\nf1 =~el_cfgs <- list(
#'     "Mediation model" = list(
#'       model = 'X =~ f1 + f2 + f3\nf1 =~ X1 + X2 + X3\nf2 =~ X4 + X5 + X6\nf3 =~ X7 + X8 + X9\nM =~ M1 + M2 + M3 + M4 + M5\nY =~ Y1 + Y2 + Y3 + Y4 + Y5\nM ~ path_a*X\nY ~ path_b*M',
#'       key_paths = c("path_a", "path_b")
#'     )
#'   )
#'
#'   res <- interative_removint simulateData
#' @importFrom stats cov
#' @export
interative_removing <- function(data,
                                model_configs,
                                var_list = NULL,
                                initial_remove_pct = 0.05,
                                pct_increment = 0.05,
                                max_iterations = 10,
                                max_total_remove_pct = 0.5,
                                seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Default thresholds (consistent with Core_file.rmd)
  criteria <- list(
    reliability = list(
      min_alpha = 0.7,
      check = TRUE
    ),
    paths = list(
      check = TRUE
    ),
    fit = list(
      min_cfi = 0.90,
      max_rmsea = 0.08,
      check = TRUE
    )
  )

  cat("========================================\n")
  cat("Start removal-only iterative cleaning (interative_removing)\n")
  cat("========================================\n\n")

  # Create original row id for tracking removals
  data$..row_id <- seq_len(nrow(data))
  current_data <- data
  initial_n <- nrow(current_data)

  removal_proportion <- initial_remove_pct
  iteration <- 0
  removed_row_ids <- integer(0)

  # Initial evaluation
  cat("[Step 1] Evaluate initial data\n")
  cat("Initial sample size:", nrow(current_data), "\n\n")
  check_results <- check_all_criteria(
    data = current_data,
    model_configs = model_configs,
    var_list = var_list,
    criteria = criteria
  )

  all_pass <- isTRUE(check_results$reliability$pass) &&
    isTRUE(check_results$paths$pass) &&
    isTRUE(check_results$fit$pass)

  if (all_pass) {
    cat("\n✓✓✓ Initial data already meet all criteria — no removal needed! ✓✓✓\n")
    current_data$..row_id <- NULL
    return(list(
      cleaned_data = current_data,
      iterations = 0,
      initial_n = initial_n,
      final_n = nrow(current_data),
      total_removed = 0,
      removed_row_ids = removed_row_ids,
      removal_proportion = 0,
      final_results = check_results
    ))
  }

  cat("\n✗✗✗ Criteria not met; starting iterative removal...\n\n")

  while (iteration < max_iterations) {
    iteration <- iteration + 1
    cat("========================================\n")
    cat("[Iteration ", iteration, "] Current removal rate: ", removal_proportion * 100, "%\n", sep = "")
    cat("========================================\n")

    # Guard: total removal cap
    total_removed_so_far <- length(removed_row_ids)
    total_removed_ratio <- total_removed_so_far / initial_n
    if (total_removed_ratio >= max_total_remove_pct) {
      cat("⚠️ Total removal proportion reached limit ", max_total_remove_pct * 100, "%, stopping\n\n", sep = "")
      break
    }

    # Target removals this round (based on current n)
    n_to_remove <- ceiling(removal_proportion * nrow(current_data))
    # Do not exceed remaining allowed removals
    remaining_allow <- max(0, floor(max_total_remove_pct * initial_n) - total_removed_so_far)
    if (n_to_remove > remaining_allow) n_to_remove <- remaining_allow
    if (n_to_remove <= 0) {
      cat("⚠️ No removal quota left, stopping\n\n")
      break
    }
    cat("Target removals:", n_to_remove, "\n")

    if (!isTRUE(check_results$fit$pass)) {
      focus <- "fit"
    } else if (!isTRUE(check_results$paths$pass)) {
      focus <- "paths"
    } else if (!isTRUE(check_results$reliability$pass)) {
      focus <- "reliability"
    } else {
      focus <- "auto"
    }
    cat("Focus this round: ", switch(focus,
                                     "reliability" = "Reliability",
                                     "paths" = "Path significance",
                                     "fit" = "Model fit",
                                     "auto" = "Auto"), "\n", sep = "")
    bad_cases_local_idx <- tryCatch({
      identify_bad_cases(
        data = current_data,
        model_configs = model_configs,
        var_list = var_list,
        criteria = criteria,
        n_to_remove = n_to_remove,
        focus = focus
      )
    }, error = function(e) integer(0))
    if (length(bad_cases_local_idx) == 0) {
      cat("⚠️ No cases identified, stopping\n\n")
      break
    }

    # Map local row indices back to original row_id
    new_removed_ids <- current_data$..row_id[bad_cases_local_idx]
    removed_row_ids <- c(removed_row_ids, new_removed_ids)

    # Remove cases
    data_cleaned <- current_data[-bad_cases_local_idx, , drop = FALSE]
    cat("Removed cases:", length(bad_cases_local_idx), "\n")
    cat("Sample size after removal:", nrow(data_cleaned), "\n")

    current_data <- data_cleaned

    # Re-evaluate
    cat("\n[Re-evaluate]\n")
    check_results <- check_all_criteria(
      data = current_data,
      model_configs = model_configs,
      var_list = var_list,
      criteria = criteria
    )

    all_pass <- isTRUE(check_results$reliability$pass) &&
      isTRUE(check_results$paths$pass) &&
      isTRUE(check_results$fit$pass)

    cat(sprintf("\nCurrent status: Reliability[%s] Paths[%s] Fit[%s]\n\n",
                ifelse(isTRUE(check_results$reliability$pass), "✓", "✗"),
                ifelse(isTRUE(check_results$paths$pass), "✓", "✗"),
                ifelse(isTRUE(check_results$fit$pass), "✓", "✗")))

    if (all_pass) {
      cat("\n✓✓✓ All criteria satisfied! Removal completed! ✓✓✓\n")
      break
    }

    # Increase removal rate for next round
    cat("✗ Still not satisfied; increasing removal rate...\n\n")
    removal_proportion <- removal_proportion + pct_increment

    # Guard: stop when sample size is too small
    if (nrow(current_data) < 20) {
      cat("⚠️ Sample size too small (<20), stopping\n\n")
      break
    }
  }

  # Return results
  cat("========================================\n")
  cat("Removal process finished\n")
  cat("========================================\n")
  cat("Initial sample size:", initial_n, "\n")
  cat("Final sample size:", nrow(current_data), "\n")
  cat("Iterations:", iteration, "\n")
  cat("Total removal proportion:", round(length(removed_row_ids) / initial_n * 100, 2), "%\n")

  current_data$..row_id <- NULL
  list(
    cleaned_data = current_data,
    iterations = iteration,
    initial_n = initial_n,
    final_n = nrow(current_data),
    total_removed = length(removed_row_ids),
    removed_row_ids = removed_row_ids,
    removal_proportion = removal_proportion,
    final_results = check_results
  )
}

#
