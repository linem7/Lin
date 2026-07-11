# iterative_removing2: revised removal-only iterative cleaning
#
# Differences from interative_removing():
# - Per-model fit checking: each model_configs entry may set check_fit = FALSE
#   (e.g., total-score structural models); saturated models (df = 0) are
#   auto-skipped with a note.
# - New common method bias (CMB) check: Harman single-factor test on all core
#   items (first unrotated factor variance must not exceed a threshold).
# - min_alpha_single is honored (bug fix).
# - JOINT multi-criteria case selection (no semfindr dependency): every
#   failing criterion contributes a normalized per-case gain vector
#   (reliability -> LOO delta-alpha; paths -> exact dfbeta regression
#   influence on the full proxy equation, supporting controls and A:B
#   interactions; fit -> Mahalanobis distance; cmb -> LOO first-eigenvalue
#   reduction); vectors are summed and harm to currently-significant key
#   paths is subtracted as a protection penalty.
# - Best-state checkpointing via an aggregate distance-to-pass (total_gap2):
#   if the loop ends without full success the best intermediate state is
#   restored, so the result is never worse than the input.
# - Fallback z-score selection only uses variables that appear in models.
# - Convergence is checked everywhere; omega has Heywood-case guards.
# - All thresholds configurable; verbose switch; no fake seed argument.
#
# NOTE (didactic use only): this family of functions removes cases to make
# results meet criteria. It is intended for simulation/teaching demonstrations
# and must not be used to clean real data for substantive analysis.

#' Internal: extract item names per factor from lavaan syntax (no fitting)
#'
#' @param model_syntax lavaan model syntax.
#' @param observed character vector of observed variable names in data.
#' @return list with `latents` (character), `items_by_factor` (named list),
#'   `first_order` (character), `is_hier` (logical).
#' @keywords internal
parse_measurement2 <- function(model_syntax, observed) {
  pt <- tryCatch(lavaan::lavaanify(model_syntax), error = function(e) NULL)
  if (is.null(pt)) return(NULL)
  lo <- pt[pt$op == "=~", , drop = FALSE]
  latents <- unique(lo$lhs)
  items_by_factor <- lapply(latents, function(f) {
    it <- lo$rhs[lo$lhs == f]
    it[it %in% observed]
  })
  names(items_by_factor) <- latents
  # first-order factors = latents whose rhs side appears under another latent
  hier_rows <- lo[lo$rhs %in% latents, , drop = FALSE]
  is_hier <- nrow(hier_rows) > 0
  first_order <- if (is_hier) unique(hier_rows$rhs) else latents
  list(latents = latents, items_by_factor = items_by_factor,
       first_order = first_order, is_hier = is_hier)
}

#' Internal: Cronbach's alpha from a covariance matrix
#' @keywords internal
alpha_from_cov2 <- function(S) {
  if (is.null(S) || any(!is.finite(S))) return(NA_real_)
  p <- ncol(S)
  if (is.null(p) || p < 2) return(NA_real_)
  denom <- sum(S)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  p / (p - 1) * (1 - sum(diag(S)) / denom)
}

#' Internal: omega-like reliability from standardized loadings, with guards
#' @return list(value, note) — value is NA when loadings are unusable.
#' @keywords internal
omega_from_loadings2 <- function(loadings) {
  if (length(loadings) < 2 || any(!is.finite(loadings))) {
    return(list(value = NA_real_, note = "insufficient/non-finite loadings"))
  }
  if (any(abs(loadings) >= 1)) {
    return(list(value = NA_real_, note = "Heywood case (|std loading| >= 1)"))
  }
  if (any(loadings > 0) && any(loadings < 0)) {
    note <- "mixed loading signs (check reverse-scored items)"
  } else {
    note <- NULL
  }
  sl <- sum(abs(loadings))
  err <- length(loadings) - sum(loadings^2)
  if (err <= 0) return(list(value = NA_real_, note = "non-positive error variance"))
  list(value = sl^2 / (sl^2 + err), note = note)
}

#' Internal: reliability per construct (revised)
#'
#' Same purpose as `calculate_reliability_alpha()` but: checks convergence,
#' guards omega against Heywood cases, and honors a separate single-factor
#' threshold.
#'
#' @param data data.frame of observed variables.
#' @param var_list named list of lavaan measurement syntax per construct.
#' @param criteria_alpha_dim threshold for first-order dimensions.
#' @param criteria_alpha_overall2 threshold for pooled second-order alpha.
#' @param criteria_alpha_single threshold for single-factor constructs.
#' @return named list per construct: `type`, `pass`, `alpha_overall`,
#'   `n_dimensions`/`n_items`, `dimensions`, optional `note`/`error`.
#' @keywords internal
calculate_reliability_alpha2 <- function(data, var_list,
                                         criteria_alpha_dim = 0.7,
                                         criteria_alpha_overall2 = 0.8,
                                         criteria_alpha_single = 0.7) {
  results <- list()
  obs <- colnames(data)

  for (construct_name in names(var_list)) {
    model_syntax <- var_list[[construct_name]]

    fit <- tryCatch(suppressWarnings(cfa(model_syntax, data = data, std.lv = TRUE)),
                    error = function(e) NULL)
    converged <- !is.null(fit) &&
      isTRUE(tryCatch(lavInspect(fit, "converged"), error = function(e) FALSE))

    if (!converged) {
      results[[construct_name]] <- list(
        type = "Unknown", pass = FALSE, alpha_overall = NA_real_,
        dimensions = NULL, error = "Model did not fit/converge"
      )
      next
    }

    pe <- tryCatch(parameterEstimates(fit, standardized = TRUE), error = function(e) NULL)
    ms <- parse_measurement2(model_syntax, obs)
    if (is.null(pe) || is.null(ms)) {
      results[[construct_name]] <- list(
        type = "Unknown", pass = FALSE, alpha_overall = NA_real_,
        dimensions = NULL, error = "Could not extract estimates"
      )
      next
    }

    get_std <- function(factor, items) {
      rows <- pe[pe$op == "=~" & pe$lhs == factor & pe$rhs %in% items, , drop = FALSE]
      rows$std.all
    }

    if (ms$is_hier) {
      # === Second-order model: omega per first-order dim + pooled alpha ===
      dim_results <- list()
      notes <- character(0)
      for (f1 in ms$first_order) {
        items <- ms$items_by_factor[[f1]]
        if (length(items) >= 2) {
          om <- omega_from_loadings2(get_std(f1, items))
          if (!is.null(om$note)) notes <- c(notes, paste0(f1, ": ", om$note))
          dim_results[[f1]] <- list(alpha = om$value, n_items = length(items))
        } else {
          dim_results[[f1]] <- list(alpha = NA_real_, n_items = length(items))
        }
      }

      all_items <- unique(unlist(ms$items_by_factor[ms$first_order], use.names = FALSE))
      all_items <- intersect(all_items, obs)
      alpha_total <- if (length(all_items) >= 2) {
        S <- tryCatch(stats::cov(data[, all_items, drop = FALSE],
                                 use = "pairwise.complete.obs"),
                      error = function(e) NULL)
        alpha_from_cov2(S)
      } else NA_real_

      dim_alphas <- vapply(dim_results, function(x) {
        if (is.na(x$alpha)) -Inf else x$alpha
      }, numeric(1))
      all_dim_pass <- length(dim_alphas) > 0 && all(dim_alphas >= criteria_alpha_dim)
      overall_pass <- !is.na(alpha_total) && alpha_total >= criteria_alpha_overall2

      results[[construct_name]] <- list(
        type = "Second-order factor",
        pass = all_dim_pass && overall_pass,
        alpha_overall = alpha_total,
        n_dimensions = length(ms$first_order),
        dimensions = dim_results,
        note = if (length(notes) > 0) paste(notes, collapse = "; ") else NULL
      )
    } else {
      # === Single-factor model ===
      f <- ms$latents[1]
      items <- ms$items_by_factor[[f]]
      if (length(items) >= 2) {
        om <- omega_from_loadings2(get_std(f, items))
        results[[construct_name]] <- list(
          type = "Single-dimension construct",
          pass = !is.na(om$value) && om$value >= criteria_alpha_single,
          alpha_overall = om$value,
          n_items = length(items),
          dimensions = NULL,
          note = om$note
        )
      } else {
        results[[construct_name]] <- list(
          type = "Single-dimension construct",
          pass = FALSE, alpha_overall = NA_real_,
          n_items = length(items), dimensions = NULL
        )
      }
    }
  }
  results
}

#' Internal: Harman single-factor test for common method bias
#'
#' Pools all items of the core constructs and extracts the first unrotated
#' factor; CMB is flagged when its explained variance exceeds `max_variance`.
#'
#' @param data data.frame.
#' @param items character vector of item names to pool.
#' @param max_variance numeric threshold for first-factor variance (default 0.40).
#' @param method "pca" (eigen of correlation matrix, default) or "fa"
#'   (ML single-factor via `stats::factanal`).
#' @return list(pass, first_factor_variance, n_items, method, error).
#' @keywords internal
check_cmb <- function(data, items, max_variance = 0.40, method = c("pca", "fa")) {
  method <- match.arg(method)
  items <- intersect(items, colnames(data))
  if (length(items) < 3) {
    return(list(pass = NA, first_factor_variance = NA_real_,
                n_items = length(items), method = method,
                error = "Fewer than 3 items available for CMB test"))
  }
  X <- data[, items, drop = FALSE]
  X <- X[stats::complete.cases(X), , drop = FALSE]
  if (nrow(X) < length(items) + 1) {
    return(list(pass = NA, first_factor_variance = NA_real_,
                n_items = length(items), method = method,
                error = "Too few complete cases for CMB test"))
  }

  prop1 <- tryCatch({
    if (method == "pca") {
      R <- stats::cor(X)
      ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
      ev[1] / length(items)
    } else {
      fa <- stats::factanal(X, factors = 1)
      sum(fa$loadings[, 1]^2) / length(items)
    }
  }, error = function(e) NA_real_)

  if (!is.finite(prop1)) {
    return(list(pass = NA, first_factor_variance = NA_real_,
                n_items = length(items), method = method,
                error = "CMB computation failed"))
  }
  list(pass = prop1 <= max_variance, first_factor_variance = prop1,
       n_items = length(items), method = method, error = NULL)
}

#' Internal: collect all items appearing in var_list measurement models
#' @keywords internal
all_items_from_var_list2 <- function(var_list, observed) {
  if (is.null(var_list)) return(character(0))
  items <- unlist(lapply(var_list, function(m) {
    ms <- parse_measurement2(m, observed)
    if (is.null(ms)) return(character(0))
    unlist(ms$items_by_factor, use.names = FALSE)
  }), use.names = FALSE)
  unique(items)
}

#' Internal: unified evaluation (revised)
#'
#' Evaluates four criteria: reliability (per construct in `var_list`), key
#' path significance, fit indices, and common method bias. Each model in
#' `model_configs` is fitted ONCE and the fit object is reused for both the
#' path and the fit checks. A model with `check_fit = FALSE`, or a saturated
#' model (df = 0), is excluded from the fit check.
#'
#' @param data data.frame of observed variables (may contain `..row_id`).
#' @param model_configs named list; each element has `model` (lavaan syntax),
#'   optional `key_paths` (labeled paths), optional `check_fit` (logical,
#'   default TRUE).
#' @param var_list optional named list of measurement models.
#' @param criteria list of thresholds/flags (see [iterative_removing2()]).
#' @param verbose logical; print a report.
#' @return list with sublists `reliability`, `paths`, `fit`, `cmb`, each with
#'   `pass` (logical or NA when not checked) and `details`.
#' @keywords internal
check_all_criteria2 <- function(data, model_configs, var_list = NULL,
                                criteria, verbose = TRUE) {
  say <- function(...) if (verbose) cat(...)
  obs <- setdiff(colnames(data), "..row_id")

  results <- list(
    reliability = list(pass = NA, details = NULL),
    paths       = list(pass = NA, details = NULL),
    fit         = list(pass = NA, details = NULL),
    cmb         = list(pass = NA, details = NULL)
  )

  # === Fit all structural models once ===
  fits <- list()
  for (config_name in names(model_configs)) {
    config <- model_configs[[config_name]]
    fit <- tryCatch(suppressWarnings(cfa(config$model, data = data, std.lv = TRUE)),
                    error = function(e) NULL)
    converged <- !is.null(fit) &&
      isTRUE(tryCatch(lavInspect(fit, "converged"), error = function(e) FALSE))
    fits[[config_name]] <- if (converged) fit else NULL
  }

  # === 1. Reliability ===
  if (!is.null(var_list) && isTRUE(criteria$reliability$check)) {
    say("\n[Reliability check]\n")
    min_alpha <- criteria$reliability$min_alpha
    min_single <- if (!is.null(criteria$reliability$min_alpha_single)) {
      criteria$reliability$min_alpha_single
    } else min_alpha
    min_overall2 <- if (!is.null(criteria$reliability$min_alpha_overall2)) {
      criteria$reliability$min_alpha_overall2
    } else 0.8

    rel_results <- tryCatch(
      calculate_reliability_alpha2(
        data[, obs, drop = FALSE], var_list,
        criteria_alpha_dim = min_alpha,
        criteria_alpha_overall2 = min_overall2,
        criteria_alpha_single = min_single
      ),
      error = function(e) NULL
    )

    if (!is.null(rel_results) && length(rel_results) > 0) {
      for (nm in names(rel_results)) {
        r <- rel_results[[nm]]
        if (identical(r$type, "Second-order factor")) {
          say(sprintf("  [%s] second-order: pooled α=%s %s\n", nm,
                      ifelse(is.na(r$alpha_overall), "NA",
                             sprintf("%.3f", r$alpha_overall)),
                      ifelse(isTRUE(r$pass), "✓", "✗")))
          for (dn in names(r$dimensions)) {
            d <- r$dimensions[[dn]]
            say(sprintf("      %s: α=%s (n=%d) %s\n", dn,
                        ifelse(is.na(d$alpha), "NA", sprintf("%.3f", d$alpha)),
                        d$n_items,
                        ifelse(!is.na(d$alpha) && d$alpha >= min_alpha,
                               "✓", "✗")))
          }
        } else {
          say(sprintf("  [%s] %s: α/ω=%s %s\n", nm, r$type,
                      ifelse(is.na(r$alpha_overall), "NA",
                             sprintf("%.3f", r$alpha_overall)),
                      ifelse(isTRUE(r$pass), "✓", "✗")))
        }
        if (!is.null(r$note)) say("      note: ", r$note, "\n", sep = "")
        if (!is.null(r$error)) say("      error: ", r$error, "\n", sep = "")
      }
      results$reliability <- list(
        pass = all(vapply(rel_results, function(x) isTRUE(x$pass), logical(1))),
        details = rel_results
      )
    } else {
      say("  ✗ Reliability analysis failed\n")
      results$reliability <- list(pass = FALSE, details = NULL)
    }
  }

  # === 2. Path significance ===
  if (isTRUE(criteria$paths$check)) {
    say("\n[Path significance check]\n")
    alpha_level <- if (!is.null(criteria$paths$alpha_level)) {
      criteria$paths$alpha_level
    } else 0.05

    all_paths_pass <- TRUE
    path_details <- list()
    any_path_checked <- FALSE

    for (config_name in names(model_configs)) {
      config <- model_configs[[config_name]]
      if (is.null(config$key_paths)) next
      fit <- fits[[config_name]]
      if (is.null(fit)) {
        say("  ✗ Model ", config_name, " failed to fit/converge\n", sep = "")
        all_paths_pass <- FALSE
        next
      }
      pe <- tryCatch(parameterEstimates(fit), error = function(e) NULL)
      if (is.null(pe)) { all_paths_pass <- FALSE; next }

      for (path_spec in config$key_paths) {
        any_path_checked <- TRUE
        # Optional direction prefix: "+lab" requires a POSITIVE significant
        # coefficient, "-lab" a NEGATIVE one; no prefix = direction-agnostic
        expected_sign <- if (startsWith(path_spec, "+")) 1
        else if (startsWith(path_spec, "-")) -1 else 0
        path_label <- sub("^[+-]", "", path_spec)

        path_row <- pe[pe$label == path_label, , drop = FALSE]
        if (nrow(path_row) == 0) {
          say("  ⚠ Path ", path_label, " not found in ", config_name, "\n", sep = "")
          all_paths_pass <- FALSE
          next
        }
        path_row <- path_row[1, ]
        pv <- suppressWarnings(as.numeric(path_row$pvalue))
        zv <- suppressWarnings(as.numeric(path_row$z))
        is_sig <- is.finite(pv) && pv < alpha_level
        sign_ok <- expected_sign == 0 ||
          (is.finite(path_row$est) && sign(path_row$est) == expected_sign)
        path_pass <- is_sig && sign_ok

        dir_note <- if (expected_sign == 0) "" else {
          paste0(", expected ", ifelse(expected_sign > 0, "positive", "negative"),
                 ifelse(sign_ok, "", " [WRONG SIGN]"))
        }
        say(sprintf("  %s %s (%s): est = %.3f, p = %s%s\n",
                    ifelse(path_pass, "✓", "✗"), path_label, config_name,
                    path_row$est,
                    ifelse(is.finite(pv), sprintf("%.3f", pv), "NA"),
                    dir_note))
        if (!path_pass) all_paths_pass <- FALSE
        path_details[[paste(config_name, path_label, sep = ".")]] <- list(
          model = config_name, label = path_label,
          estimate = path_row$est, pvalue = pv, z = zv,
          expected_sign = expected_sign, sign_ok = sign_ok,
          significant = is_sig, pass = path_pass
        )
      }
    }

    results$paths <- list(
      pass = if (any_path_checked || !all_paths_pass) all_paths_pass else NA,
      details = path_details
    )
  }

  # === 3. Fit indices (per-model opt-out + saturated detection) ===
  if (isTRUE(criteria$fit$check)) {
    say("\n[Fit indices check]\n")
    all_fit_pass <- TRUE
    any_fit_checked <- FALSE
    fit_details <- list()

    for (config_name in names(model_configs)) {
      config <- model_configs[[config_name]]
      if (isFALSE(config$check_fit)) {
        say("  - ", config_name, ": fit check disabled for this model\n", sep = "")
        next
      }
      fit <- fits[[config_name]]
      if (is.null(fit)) {
        say("  ✗ Model ", config_name, " failed to fit/converge\n", sep = "")
        all_fit_pass <- FALSE
        next
      }
      fm <- tryCatch(fitMeasures(fit, c("df", "cfi", "tli", "rmsea", "srmr")),
                     error = function(e) NULL)
      if (is.null(fm)) { all_fit_pass <- FALSE; next }
      if (fm[["df"]] == 0) {
        say("  - ", config_name, ": saturated model (df = 0), fit check skipped\n",
            sep = "")
        fit_details[[config_name]] <- c(as.list(fm), list(saturated = TRUE))
        next
      }
      any_fit_checked <- TRUE
      cfi_pass <- is.finite(fm[["cfi"]]) && fm[["cfi"]] >= criteria$fit$min_cfi
      rmsea_pass <- is.finite(fm[["rmsea"]]) && fm[["rmsea"]] <= criteria$fit$max_rmsea
      say(sprintf("  %s %s: CFI = %.3f (≥ %.2f), RMSEA = %.3f (≤ %.2f)\n",
                  ifelse(cfi_pass && rmsea_pass, "✓", "✗"), config_name,
                  fm[["cfi"]], criteria$fit$min_cfi,
                  fm[["rmsea"]], criteria$fit$max_rmsea))
      if (!cfi_pass || !rmsea_pass) all_fit_pass <- FALSE
      fit_details[[config_name]] <- c(as.list(fm), list(saturated = FALSE))
    }

    results$fit <- list(
      pass = if (any_fit_checked || !all_fit_pass) all_fit_pass else NA,
      details = fit_details
    )
  }

  # === 4. Common method bias (Harman single-factor) ===
  if (isTRUE(criteria$cmb$check)) {
    say("\n[Common method bias check (Harman)]\n")
    cmb_items <- if (!is.null(criteria$cmb$items)) {
      criteria$cmb$items
    } else {
      all_items_from_var_list2(var_list, obs)
    }
    max_var <- if (!is.null(criteria$cmb$max_variance_first)) {
      criteria$cmb$max_variance_first
    } else 0.40
    cmb_method <- if (!is.null(criteria$cmb$method)) criteria$cmb$method else "pca"

    cmb <- check_cmb(data[, obs, drop = FALSE], cmb_items,
                     max_variance = max_var, method = cmb_method)
    if (!is.null(cmb$error)) {
      say("  ⚠ ", cmb$error, "\n", sep = "")
    } else {
      say(sprintf("  %s First factor explains %.1f%% of variance (criterion: ≤ %.0f%%, %d items, %s)\n",
                  ifelse(isTRUE(cmb$pass), "✓", "✗"),
                  cmb$first_factor_variance * 100, max_var * 100,
                  cmb$n_items, cmb$method))
    }
    results$cmb <- list(pass = cmb$pass, details = cmb)
  }

  results
}

#' Internal: identify problematic cases via joint multi-criteria scoring
#'
#' Computes one per-case gain vector per FAILING criterion, normalizes each
#' (units are not commensurable), sums them into a joint score, subtracts a
#' penalty for harming currently-significant key paths, and selects the
#' top-k. Gain vectors (all self-contained, no semfindr):
#' - reliability: LOO case-deleted alpha improvement per failing construct,
#'   weighted by that construct's distance below its threshold.
#' - paths: EXACT LOO regression influence (`stats::dfbeta`) of each case on
#'   the target coefficient, computed on the FULL proxy equation of the
#'   dependent variable — so control variables, parallel predictors, and
#'   interaction terms (`A:B`) are partialled out correctly. The dfbeta sign
#'   convention is self-calibrated empirically at run time.
#' - fit: Mahalanobis distance on the failing models' observed variables
#'   (cases beyond the median distance only).
#' - cmb: LOO reduction of the first-eigenvalue proportion (Harman).
#' - fallback: largest |z| sum across MODEL variables when no targeted gains.
#'
#' @param data data.frame with `..row_id` column.
#' @param model_configs named list of model configurations.
#' @param var_list named list of measurement models.
#' @param criteria list of thresholds/flags.
#' @param n_to_remove integer, number of cases to select.
#' @param check_results last evaluation from `check_all_criteria2()` (used to
#'   know which constructs/models/paths are failing).
#' @param verbose logical.
#' @return integer vector of row indices relative to `data`.
#' @keywords internal
identify_bad_cases2 <- function(data, model_configs, var_list, criteria,
                                n_to_remove, check_results,
                                verbose = TRUE) {
  say <- function(...) if (verbose) cat(...)
  obs <- setdiff(colnames(data), "..row_id")
  n <- nrow(data)

  # Restrict every strategy to variables that actually appear in the models
  model_vars <- unique(c(
    all_items_from_var_list2(var_list, obs),
    unlist(lapply(model_configs, function(cf) {
      ms <- parse_measurement2(cf$model, obs)
      pt <- tryCatch(lavaan::lavaanify(cf$model), error = function(e) NULL)
      c(if (!is.null(ms)) unlist(ms$items_by_factor, use.names = FALSE),
        if (!is.null(pt)) intersect(unique(c(pt$lhs, pt$rhs)), obs))
    }), use.names = FALSE)
  ))
  model_vars <- intersect(model_vars, obs)
  num_vars <- model_vars[vapply(data[model_vars], is.numeric, logical(1))]

  # --- Generic helpers for gain-vector based selection ---
  top_k <- function(g, k) {
    if (is.null(g)) return(integer(0))
    g[!is.finite(g)] <- 0
    ok <- which(g > 0)
    if (length(ok) == 0) return(integer(0))
    head(ok[order(g[ok], decreasing = TRUE)], min(k, length(ok)))
  }

  # Proxy getter for one model: observed numeric -> the column itself;
  # latent -> row mean of its items; product term "A:B" -> product of proxies
  make_proxy <- function(model_syntax) {
    ms <- parse_measurement2(model_syntax, obs)
    getv <- function(v) {
      if (v %in% num_vars) return(data[[v]])
      if (!is.null(ms) && v %in% names(ms$items_by_factor)) {
        it <- intersect(ms$items_by_factor[[v]], num_vars)
        if (length(it) >= 1) {
          return(rowMeans(data[, it, drop = FALSE], na.rm = TRUE))
        }
      }
      if (grepl(":", v, fixed = TRUE)) {
        parts <- strsplit(v, ":", fixed = TRUE)[[1]]
        vals <- lapply(parts, getv)
        if (!any(vapply(vals, is.null, logical(1)))) return(Reduce(`*`, vals))
      }
      NULL
    }
    getv
  }

  # --- Reliability: LOO delta-alpha gain vector over failing constructs ---
  # Each construct's gain is weighted by how far its alpha is below the
  # threshold, so constructs in bigger trouble dominate the selection. The
  # per-round quota is modest and gains are recomputed every round, so the
  # outer iteration acts as a batch-level greedy refinement.
  gain_reliability <- function() {
    if (is.null(var_list)) return(NULL)
    rel <- check_results$reliability$details
    if (is.null(rel)) return(NULL)
    failing <- names(rel)[vapply(rel, function(x) !isTRUE(x$pass), logical(1))]
    if (length(failing) == 0) return(NULL)

    thr_single <- if (!is.null(criteria$reliability$min_alpha_single)) {
      criteria$reliability$min_alpha_single
    } else criteria$reliability$min_alpha

    g <- rep(0, n)
    used <- character(0)
    for (nm in failing) {
      ms <- parse_measurement2(var_list[[nm]], obs)
      if (is.null(ms)) next
      items <- unique(unlist(ms$items_by_factor, use.names = FALSE))
      items <- items[items %in% num_vars]
      if (length(items) < 2) next

      X <- data[, items, drop = FALSE]
      S0 <- tryCatch(stats::cov(X, use = "pairwise.complete.obs"),
                     error = function(e) NULL)
      a0 <- alpha_from_cov2(S0)
      if (is.na(a0)) next

      thr <- if (identical(rel[[nm]]$type, "Second-order factor")) {
        if (!is.null(criteria$reliability$min_alpha_overall2)) {
          criteria$reliability$min_alpha_overall2
        } else 0.8
      } else thr_single
      w <- 1 + max(0, thr - a0) * 10  # bigger gap -> up to ~2-3x weight

      for (i in seq_len(n)) {
        S <- tryCatch(stats::cov(X[-i, , drop = FALSE],
                                 use = "pairwise.complete.obs"),
                      error = function(e) NULL)
        ai <- alpha_from_cov2(S)
        if (!is.na(ai) && ai > a0) g[i] <- g[i] + (ai - a0) * w
      }
      used <- c(used, nm)
    }
    if (length(used) == 0 || all(g <= 0)) return(NULL)
    say("  gain[reliability]: LOO Δα on construct(s): ",
        paste(used, collapse = ", "), "\n", sep = "")
    g
  }

  # --- Paths: exact LOO regression influence via dfbeta on proxy equations ---
  #
  # For each labeled path, the FULL regression equation of its dependent
  # variable is reconstructed from the model syntax (all predictors,
  # including control variables and product terms), latents are proxied by
  # item means, and stats::dfbeta() gives the exact per-case change of the
  # target coefficient upon deletion. This handles mediation, moderation
  # (A:B products) and covariate-adjusted models, where a bivariate
  # correlation would misrank cases.
  #
  # Returns a list of per-case influence vectors: for each path,
  # `delta[i]` = b(-i) - b (change of the coefficient when case i is
  # removed), NA for cases not usable in that equation.
  path_influences <- function(labels_by_model) {
    out <- list()
    for (mname in names(labels_by_model)) {
      config <- model_configs[[mname]]
      if (is.null(config)) next
      pt <- tryCatch(lavaan::lavaanify(config$model), error = function(e) NULL)
      if (is.null(pt)) next
      proxy <- make_proxy(config$model)

      for (lab in labels_by_model[[mname]]) {
        row <- pt[pt$label == lab & pt$op == "~", , drop = FALSE]

        if (nrow(row) == 0) {
          # Labeled covariance (~~): LOO influence on the plain correlation
          row2 <- pt[pt$label == lab & pt$op == "~~", , drop = FALSE]
          if (nrow(row2) == 0) next
          row2 <- row2[1, ]
          x <- proxy(row2$rhs); y <- proxy(row2$lhs)
          if (is.null(x) || is.null(y)) next
          cc <- stats::complete.cases(x, y)
          if (sum(cc) < 5) next
          r0 <- suppressWarnings(stats::cor(x[cc], y[cc]))
          if (!is.finite(r0)) next
          delta <- rep(NA_real_, n)
          for (i in which(cc)) {
            keep <- cc; keep[i] <- FALSE
            ri <- suppressWarnings(stats::cor(x[keep], y[keep]))
            if (is.finite(ri)) delta[i] <- ri - r0
          }
          out[[paste(mname, lab, sep = ".")]] <- list(delta = delta, b = r0)
          next
        }

        row <- row[1, ]
        dv <- row$lhs
        # All predictors of this dv in the model (the full equation)
        preds <- unique(pt$rhs[pt$op == "~" & pt$lhs == dv])
        if (!(row$rhs %in% preds)) next

        y <- proxy(dv)
        Xl <- lapply(preds, proxy)
        if (is.null(y) || any(vapply(Xl, is.null, logical(1)))) next
        X <- do.call(cbind, Xl)
        colnames(X) <- paste0("V", seq_along(preds))
        target_col <- paste0("V", match(row$rhs, preds))

        df_reg <- data.frame(.y = y, X)
        cc <- stats::complete.cases(df_reg)
        if (sum(cc) < length(preds) + 3) next
        f <- tryCatch(stats::lm(.y ~ ., data = df_reg[cc, , drop = FALSE]),
                      error = function(e) NULL)
        if (is.null(f)) next
        b <- stats::coef(f)[target_col]
        if (!is.finite(b)) next
        db <- tryCatch(stats::dfbeta(f)[, target_col], error = function(e) NULL)
        if (is.null(db)) next

        # Sign self-calibration: dfbeta's convention (b - b(-i) vs b(-i) - b)
        # is verified empirically on the most influential case, so the code
        # never depends on remembering the convention.
        i0 <- which.max(abs(db))
        b1 <- tryCatch(
          stats::coef(stats::lm(.y ~ ., data = df_reg[cc, , drop = FALSE][-i0, , drop = FALSE]))[target_col],
          error = function(e) NA_real_)
        sgn <- if (is.finite(b1) && (b - b1) * db[i0] > 0) 1 else -1
        # delta = b(-i) - b
        delta_cc <- -sgn * db

        delta <- rep(NA_real_, n)
        delta[which(cc)] <- delta_cc
        out[[paste(mname, lab, sep = ".")]] <- list(
          delta = delta, b = unname(b))
      }
    }
    out
  }

  # Gain vector for failing key paths. The target direction is the path's
  # EXPECTED sign when one was declared (via a "+lab"/"-lab" prefix), so a
  # coefficient with the wrong sign is pushed THROUGH zero toward the
  # required direction; otherwise the current sign of the coefficient is
  # pushed away from zero. Weight grows with the signed-z distance to the
  # critical value (direction declared) or the p-value excess (agnostic).
  gain_paths <- function() {
    pd <- check_results$paths$details
    if (is.null(pd) || length(pd) == 0) return(NULL)
    failing <- Filter(function(x) !isTRUE(x$pass), pd)
    if (length(failing) == 0) return(NULL)
    alpha_level <- if (!is.null(criteria$paths$alpha_level)) {
      criteria$paths$alpha_level
    } else 0.05
    z_crit <- stats::qnorm(1 - alpha_level / 2)

    by_model <- split(vapply(failing, `[[`, "", "label"),
                      vapply(failing, `[[`, "", "model"))
    infl <- path_influences(lapply(by_model, unique))
    if (length(infl) == 0) return(NULL)

    g <- rep(0, n)
    for (key in names(infl)) {
      delta <- infl[[key]]$delta
      b <- infl[[key]]$b
      f <- failing[[which(vapply(failing, function(x)
        paste(x$model, x$label, sep = ".") == key, logical(1)))[1]]]

      es <- if (!is.null(f$expected_sign)) f$expected_sign else 0
      target <- if (es != 0) es else if (b != 0) sign(b) else 0
      raw <- if (target == 0) abs(delta) else target * delta
      raw[!is.finite(raw)] <- 0

      w <- if (es != 0 && is.finite(f$z)) {
        # signed requirement: distance to the critical z in the expected
        # direction (a wrong-signed path is farthest and gets most weight)
        1 + max(0, (z_crit - es * f$z) / z_crit)
      } else if (is.finite(f$pvalue)) {
        1 + max(0, f$pvalue - alpha_level)
      } else 1
      g <- g + pmax(0, raw) * w
    }
    if (all(g <= 0)) return(NULL)
    say("  gain[paths]: dfbeta influence on ", length(infl),
        " failing path(s)\n", sep = "")
    g
  }

  # Protection penalty: harm to key paths that currently PASS (significant
  # and, when a direction was declared, correctly signed). Positive entries =
  # removing this case pushes a passing coefficient back toward zero.
  penalty_protect <- function() {
    pd <- check_results$paths$details
    if (is.null(pd) || length(pd) == 0) return(NULL)
    passing <- Filter(function(x) isTRUE(x$pass), pd)
    if (length(passing) == 0) return(NULL)

    by_model <- split(vapply(passing, `[[`, "", "label"),
                      vapply(passing, `[[`, "", "model"))
    infl <- path_influences(lapply(by_model, unique))
    if (length(infl) == 0) return(NULL)

    p <- rep(0, n)
    for (key in names(infl)) {
      delta <- infl[[key]]$delta
      b <- infl[[key]]$b
      harm <- if (b == 0) rep(0, n) else -sign(b) * delta  # shrink toward 0
      harm[!is.finite(harm)] <- 0
      p <- p + pmax(0, harm)
    }
    if (all(p <= 0)) return(NULL)
    p
  }

  # --- Fit: Mahalanobis distance on failing models' variables ---
  # Multivariate outliers inflate model misfit; normalized MD serves as the
  # fit-improvement gain (exact per-case CFI/RMSEA influence would need n
  # refits per model per round).
  gain_fit <- function() {
    fd <- check_results$fit$details
    if (isTRUE(is.na(check_results$fit$pass)) ||
        isTRUE(check_results$fit$pass)) return(NULL)
    failing_models <- names(model_configs)
    if (!is.null(fd) && length(fd) > 0) {
      bad <- vapply(names(fd), function(nm) {
        d <- fd[[nm]]
        if (isTRUE(d$saturated)) return(FALSE)
        cfi_ok <- is.finite(d$cfi) && d$cfi >= criteria$fit$min_cfi
        rmsea_ok <- is.finite(d$rmsea) && d$rmsea <= criteria$fit$max_rmsea
        !(cfi_ok && rmsea_ok)
      }, logical(1))
      if (any(bad)) failing_models <- names(fd)[bad]
    }

    vars <- unique(unlist(lapply(failing_models, function(nm) {
      ms <- parse_measurement2(model_configs[[nm]]$model, obs)
      if (is.null(ms)) return(character(0))
      unlist(ms$items_by_factor, use.names = FALSE)
    }), use.names = FALSE))
    vars <- intersect(vars, num_vars)
    if (length(vars) < 2) vars <- num_vars
    if (length(vars) < 2) return(NULL)

    X <- data[, vars, drop = FALSE]
    cc <- stats::complete.cases(X)
    if (sum(cc) <= length(vars) + 1) return(NULL)
    md <- rep(0, n)
    md[cc] <- tryCatch({
      Xc <- as.matrix(X[cc, , drop = FALSE])
      stats::mahalanobis(Xc, colMeans(Xc), stats::cov(Xc))
    }, error = function(e) rep(0, sum(cc)))
    if (all(md <= 0)) return(NULL)
    # Only cases beyond the median distance count as fit-improvement targets;
    # normalize to [0, 1] so it is commensurable with other gain vectors
    thr <- stats::median(md[cc])
    g <- pmax(0, md - thr)
    if (max(g) > 0) g <- g / max(g)
    say("  gain[fit]: Mahalanobis distance on ", length(vars),
        " model variable(s)\n", sep = "")
    g
  }

  # --- CMB: LOO first-eigenvalue proportion reduction ---
  gain_cmb <- function() {
    if (isTRUE(is.na(check_results$cmb$pass)) ||
        isTRUE(check_results$cmb$pass)) return(NULL)
    items <- if (!is.null(criteria$cmb$items)) {
      criteria$cmb$items
    } else {
      all_items_from_var_list2(var_list, obs)
    }
    items <- intersect(items, num_vars)
    if (length(items) < 3) return(NULL)

    X <- data[, items, drop = FALSE]
    cc <- stats::complete.cases(X)
    if (sum(cc) < length(items) + 2) return(NULL)
    prop_first <- function(M) {
      R <- suppressWarnings(stats::cor(M))
      if (any(!is.finite(R))) return(NA_real_)
      eigen(R, symmetric = TRUE, only.values = TRUE)$values[1] / ncol(M)
    }
    Xc <- as.matrix(X)
    p0 <- prop_first(Xc[cc, , drop = FALSE])
    if (!is.finite(p0)) return(NULL)

    g <- rep(0, n)
    for (i in which(cc)) {
      keep <- cc; keep[i] <- FALSE
      pi <- prop_first(Xc[keep, , drop = FALSE])
      if (is.finite(pi) && p0 > pi) g[i] <- p0 - pi
    }
    if (all(g <= 0)) return(NULL)
    say("  gain[cmb]: LOO first-factor variance (current = ",
        sprintf("%.1f%%", p0 * 100), ")\n", sep = "")
    g
  }

  # --- Fallback: |z| sum over model variables only ---
  pick_by_zsum <- function(k) {
    if (length(num_vars) == 0) return(integer(0))
    Z <- suppressWarnings(scale(data[, num_vars, drop = FALSE]))
    row_score <- rowSums(abs(Z), na.rm = TRUE)
    ord <- order(row_score, decreasing = TRUE)
    head(ord, min(k, length(ord)))
  }

  # ============ Joint scoring across ALL failing criteria ============
  # Each gain vector is normalized to max 1 (units differ: Δα vs Δβ vs
  # eigen-proportion vs MD), then summed. Cases that help several failing
  # criteria at once rank highest. Harm to currently-significant key paths
  # enters as a subtractive penalty so fixing one criterion does not silently
  # break another (the seesaw problem in multi-criteria cleaning).
  say("\n[Identify problematic cases] joint scoring over failing criteria\n")

  norm1 <- function(g) if (!is.null(g) && max(g) > 0) g / max(g) else g

  gains <- list(
    reliability = norm1(gain_reliability()),
    paths       = norm1(gain_paths()),
    fit         = gain_fit(),   # already normalized
    cmb         = norm1(gain_cmb())
  )
  gains <- Filter(Negate(is.null), gains)

  if (length(gains) > 0) {
    score <- Reduce(`+`, gains)
    pen <- penalty_protect()
    if (!is.null(pen)) {
      pen <- norm1(pen) * 0.5 * length(gains)
      score <- score - pen
      say("  penalty: protecting currently-significant key path(s)\n")
    }
    picks <- top_k(score, n_to_remove)
    if (length(picks) > 0) {
      say("  ✓ Joint selection (", paste(names(gains), collapse = " + "),
          "): ", length(picks), " case(s)\n", sep = "")
      return(picks)
    }
  }

  say("  ↪ No targeted gains found; using |z|-sum fallback\n")
  picks <- pick_by_zsum(n_to_remove)
  say("  ✓ Fallback: ", length(picks), " case(s)\n", sep = "")
  picks
}

#' Internal: aggregate distance-to-pass across all criteria
#'
#' Sums how far each checked criterion currently is from its threshold
#' (alpha shortfalls, p-value excesses over the significance level, CFI/RMSEA
#' shortfalls, CMB variance excess). Zero means all criteria pass. The units
#' are heterogeneous, so the value is a heuristic ordering of states of the
#' SAME problem — used to keep/restore the best intermediate state — not an
#' absolute quality measure.
#'
#' @param check_results output of `check_all_criteria2()`.
#' @param criteria the resolved criteria list.
#' @return non-negative numeric scalar.
#' @keywords internal
total_gap2 <- function(check_results, criteria) {
  gap <- 0

  # Reliability: alpha shortfalls (overall + first-order dimensions)
  rel <- check_results$reliability$details
  if (!is.null(rel)) {
    thr_single <- if (!is.null(criteria$reliability$min_alpha_single)) {
      criteria$reliability$min_alpha_single
    } else criteria$reliability$min_alpha
    thr2 <- if (!is.null(criteria$reliability$min_alpha_overall2)) {
      criteria$reliability$min_alpha_overall2
    } else 0.8
    for (r in rel) {
      if (isTRUE(r$pass)) next
      thr <- if (identical(r$type, "Second-order factor")) thr2 else thr_single
      a <- r$alpha_overall
      gap <- gap + if (is.finite(a)) max(0, thr - a) else thr * 0.5
      if (!is.null(r$dimensions)) {
        for (d in r$dimensions) {
          if (is.finite(d$alpha)) {
            gap <- gap + max(0, criteria$reliability$min_alpha - d$alpha)
          }
        }
      }
    }
  } else if (identical(check_results$reliability$pass, FALSE)) {
    gap <- gap + 1
  }

  # Paths: p-value excess over the significance level; when a direction is
  # required ("+lab"/"-lab"), use the signed-z distance to the critical value
  # instead — it stays monotone while a coefficient is pushed through zero,
  # so sign-flipping progress is visible to the best-state checkpoint.
  pd <- check_results$paths$details
  al <- if (!is.null(criteria$paths$alpha_level)) criteria$paths$alpha_level else 0.05
  if (!is.null(pd) && length(pd) > 0) {
    z_crit <- stats::qnorm(1 - al / 2)
    for (p in pd) {
      passed <- if (!is.null(p$pass)) isTRUE(p$pass) else isTRUE(p$significant)
      if (passed) next
      es <- if (!is.null(p$expected_sign)) p$expected_sign else 0
      if (es != 0 && is.finite(p$z)) {
        gap <- gap + max(0, (z_crit - es * p$z) / z_crit) * 0.25
      } else {
        gap <- gap + if (is.finite(p$pvalue)) max(0, p$pvalue - al) else 1
      }
    }
  } else if (identical(check_results$paths$pass, FALSE)) {
    gap <- gap + 1
  }

  # Fit: CFI/RMSEA shortfalls of non-saturated models
  fd <- check_results$fit$details
  if (!is.null(fd) && length(fd) > 0) {
    for (d in fd) {
      if (isTRUE(d$saturated)) next
      if (is.finite(d$cfi)) gap <- gap + max(0, criteria$fit$min_cfi - d$cfi)
      if (is.finite(d$rmsea)) gap <- gap + max(0, d$rmsea - criteria$fit$max_rmsea)
    }
  } else if (identical(check_results$fit$pass, FALSE)) {
    gap <- gap + 1
  }

  # CMB: first-factor variance excess
  cd <- check_results$cmb$details
  if (identical(check_results$cmb$pass, FALSE) && !is.null(cd) &&
      is.finite(cd$first_factor_variance)) {
    mx <- if (!is.null(criteria$cmb$max_variance_first)) {
      criteria$cmb$max_variance_first
    } else 0.40
    gap <- gap + max(0, cd$first_factor_variance - mx)
  }

  gap
}

#' Removal-only iterative cleaning (revised) — FOR TEACHING/SIMULATION ONLY
#'
#' Iteratively removes cases from a dataset until all requested criteria are
#' met: (1) reliability of each construct in `var_list`, (2) significance of
#' labeled key paths, (3) fit indices of models that opt in, and (4) absence
#' of common method bias (Harman single-factor test). This is an
#' outcome-driven procedure: it exists to DEMONSTRATE to students how case
#' removal can manufacture "better" results. It must not be used to clean
#' real data for substantive analysis.
#'
#' Compared with [interative_removing()]:
#' \itemize{
#'   \item Per-model fit control: set `check_fit = FALSE` in a
#'     `model_configs` entry to exclude it from the fit check (e.g., a
#'     total-score structural model); saturated models (df = 0) are skipped
#'     automatically. Typical usage fits the measurement CFA of each core
#'     construct as its own model (fit checked) plus one total-score
#'     structural model (paths checked, fit optional).
#'   \item New CMB criterion: pooled items of all constructs must not have a
#'     first unrotated factor explaining more than
#'     `criteria$cmb$max_variance_first` (default 40\%).
#'   \item `min_alpha_single` is honored; all thresholds (including the path
#'     significance level) are configurable.
#'   \item JOINT multi-criteria case selection: every failing criterion
#'     contributes a normalized per-case gain vector (LOO delta-alpha for
#'     reliability; exact `dfbeta` regression influence for paths; Mahalanobis
#'     distance for fit; LOO first-eigenvalue reduction for CMB), the vectors
#'     are summed, and harm to currently-significant key paths is subtracted
#'     as a protection penalty — so several simultaneously-failing criteria
#'     (e.g., a mediation model and a moderation model on the same data) are
#'     optimized together instead of see-sawing between single foci.
#'   \item Path influence uses the FULL equation of each labeled path's
#'     dependent variable (all predictors, control variables, and `A:B`
#'     product terms, with latents proxied by item means), so
#'     covariate-adjusted and moderated paths are targeted correctly.
#'   \item Best-state checkpointing: an aggregate gap (distance-to-pass over
#'     all criteria, see [total_gap2()]) is tracked each round; if the loop
#'     ends without full success, the best intermediate state is restored,
#'     so the result is never worse than the untouched input.
#'   \item Each model is fitted once per evaluation round; convergence is
#'     always checked; omega has Heywood-case guards. The fallback |z|-sum
#'     only uses variables appearing in the models.
#'   \item Single-round removals are capped so the sample never drops below
#'     `min_n`; a `verbose` switch silences all output.
#' }
#'
#' @param data data.frame of observed variables.
#' @param model_configs named list; each element is a list with:
#'   \itemize{
#'     \item `model`: lavaan syntax (required)
#'     \item `key_paths`: character vector of labeled paths to test
#'       (optional). A label may carry a direction prefix: `"+cp"` requires
#'       the coefficient `cp` to be POSITIVE and significant, `"-cp"`
#'       NEGATIVE and significant; a bare label only requires significance.
#'       With a declared direction, case selection pushes the coefficient
#'       toward that sign — a wrong-signed path (e.g., currently negative
#'       when `"+cp"` is requested) is driven through zero and on to
#'       significance on the other side, given enough removable cases.
#'     \item `check_fit`: logical, include this model in the fit check
#'       (optional, default TRUE)
#'   }
#' @param var_list optional named list of measurement models (lavaan syntax)
#'   keyed by construct name; required for the reliability check and the
#'   default CMB item pool.
#' @param initial_remove_pct numeric in (0,1], removal proportion in round 1
#'   (default 0.05).
#' @param pct_increment numeric, increment added to the removal proportion
#'   after each unsuccessful round (default 0.05).
#' @param pct_basis how the per-round removal proportion is interpreted:
#'   \itemize{
#'     \item `"initial"` (default): the proportion is a CUMULATIVE target
#'       relative to the initial sample size. Round r removes
#'       `ceiling(prop_r * initial_n) - already_removed` cases. E.g., with
#'       n = 500, 5\% + 5\%: round 1 removes 25 (475 left), round 2 removes
#'       25 more so that the cumulative removal is 10\% of 500 (450 left).
#'     \item `"current"`: the proportion applies to the CURRENT sample each
#'       round (legacy behavior of [interative_removing()]). Same example:
#'       round 2 removes `ceiling(0.10 * 475) = 48` cases (427 left), so the
#'       cumulative removal exceeds the nominal percentage.
#'   }
#' @param max_iterations integer, maximum number of rounds (default 10).
#' @param max_total_remove_pct numeric in (0,1], cap on cumulative removals
#'   relative to the initial n (default 0.5).
#' @param min_n integer, hard floor on the remaining sample size; per-round
#'   removals are capped so `nrow(data)` never drops below it (default 20).
#' @param criteria optional list overriding defaults:
#'   \preformatted{
#'   list(
#'     reliability = list(check = TRUE, min_alpha = 0.7,
#'                        min_alpha_single = 0.7, min_alpha_overall2 = 0.8),
#'     paths       = list(check = TRUE, alpha_level = 0.05),
#'     fit         = list(check = TRUE, min_cfi = 0.90, max_rmsea = 0.08),
#'     cmb         = list(check = FALSE, max_variance_first = 0.40,
#'                        method = "pca", items = NULL)
#'   )
#'   }
#'   `cmb$items` defaults to all items found in `var_list`. Fields not
#'   provided keep their defaults.
#' @param verbose logical, print progress reports (default TRUE).
#'
#' @return list with elements `cleaned_data`, `iterations`, `initial_n`,
#'   `final_n`, `total_removed`, `removed_row_ids` (positions in the original
#'   data), `all_pass` (logical), `final_gap` (aggregate distance-to-pass of
#'   the returned state; 0 when `all_pass`), `reverted_to_iteration` (the
#'   iteration whose state was restored when the final state was worse, NA
#'   otherwise), `final_results` (evaluation of the returned state), and
#'   `history` (per-round data.frame: iteration, n_removed, n_after, gap,
#'   and pass flags).
#'
#' @examples
#' \donttest{
#' # Chained mediation with 4 constructs: X -> M1 -> M2 -> Y
#' set.seed(123)
#' pop_model <- '
#'   X  =~ 0.8*X1 + 0.8*X2 + 0.8*X3
#'   M1 =~ 0.8*M11 + 0.8*M12 + 0.8*M13
#'   M2 =~ 0.8*M21 + 0.8*M22 + 0.8*M23
#'   Y  =~ 0.8*Y1 + 0.8*Y2 + 0.8*Y3
#'   M1 ~ 0.5*X
#'   M2 ~ 0.5*M1
#'   Y  ~ 0.5*M2
#' '
#' dat <- as.data.frame(lavaan::simulateData(pop_model, sample.nobs = 400))
#'
#' # Contaminate with careless responders
#' set.seed(999)
#' noise <- as.data.frame(sapply(dat, function(x)
#'   rnorm(40, mean = mean(x), sd = sd(x) * 3)))
#' dat <- rbind(dat, noise)
#'
#' var_list <- list(
#'   X  = 'X  =~ X1 + X2 + X3',
#'   M1 = 'M1 =~ M11 + M12 + M13',
#'   M2 = 'M2 =~ M21 + M22 + M23',
#'   Y  = 'Y  =~ Y1 + Y2 + Y3'
#' )
#'
#' # Total scores for the structural model
#' dat$X_sum  <- rowMeans(dat[, c("X1", "X2", "X3")])
#' dat$M1_sum <- rowMeans(dat[, c("M11", "M12", "M13")])
#' dat$M2_sum <- rowMeans(dat[, c("M21", "M22", "M23")])
#' dat$Y_sum  <- rowMeans(dat[, c("Y1", "Y2", "Y3")])
#'
#' # 5 models: 4 measurement CFAs (fit checked) + 1 total-score chain
#' # (paths checked, fit not checked — it is saturated anyway)
#' model_configs <- list(
#'   CFA_X  = list(model = var_list$X),
#'   CFA_M1 = list(model = var_list$M1),
#'   CFA_M2 = list(model = var_list$M2),
#'   CFA_Y  = list(model = var_list$Y),
#'   Chain  = list(
#'     model = '
#'       M1_sum ~ a*X_sum
#'       M2_sum ~ b*M1_sum
#'       Y_sum  ~ c*M2_sum
#'     ',
#'     key_paths = c("a", "b", "c"),
#'     check_fit = FALSE
#'   )
#' )
#'
#' res <- iterative_removing2(
#'   data = dat,
#'   model_configs = model_configs,
#'   var_list = var_list,
#'   criteria = list(
#'     reliability = list(check = TRUE, min_alpha = 0.7),
#'     paths = list(check = TRUE),
#'     fit = list(check = TRUE, min_cfi = 0.90, max_rmsea = 0.08),
#'     cmb = list(check = TRUE, max_variance_first = 0.40)
#'   )
#' )
#' res$history
#' }
#'
#' @seealso interative_removing, calculate_reliability_alpha2,
#'   check_all_criteria2, identify_bad_cases2, check_cmb
#' @importFrom stats cov cor complete.cases mahalanobis factanal
#' @export
iterative_removing2 <- function(data,
                                model_configs,
                                var_list = NULL,
                                initial_remove_pct = 0.05,
                                pct_increment = 0.05,
                                max_iterations = 10,
                                max_total_remove_pct = 0.5,
                                min_n = 20,
                                pct_basis = c("initial", "current"),
                                criteria = NULL,
                                verbose = TRUE) {
  say <- function(...) if (verbose) cat(...)
  pct_basis <- match.arg(pct_basis)

  stopifnot(is.data.frame(data), nrow(data) > 0)
  if (initial_remove_pct <= 0 || initial_remove_pct > 1) {
    stop("initial_remove_pct must be in (0, 1]")
  }
  if (max_total_remove_pct <= 0 || max_total_remove_pct > 1) {
    stop("max_total_remove_pct must be in (0, 1]")
  }

  default_criteria <- list(
    reliability = list(check = TRUE, min_alpha = 0.7,
                       min_alpha_single = 0.7, min_alpha_overall2 = 0.8),
    paths       = list(check = TRUE, alpha_level = 0.05),
    fit         = list(check = TRUE, min_cfi = 0.90, max_rmsea = 0.08),
    cmb         = list(check = FALSE, max_variance_first = 0.40,
                       method = "pca", items = NULL)
  )
  criteria <- if (is.null(criteria)) {
    default_criteria
  } else {
    utils::modifyList(default_criteria, criteria)
  }

  say("========================================\n")
  say("Iterative removal (v2) — teaching/simulation use only\n")
  say("========================================\n")

  data$..row_id <- seq_len(nrow(data))
  current_data <- data
  initial_n <- nrow(current_data)
  removal_proportion <- initial_remove_pct
  iteration <- 0
  removed_row_ids <- integer(0)
  history <- list()

  criterion_pass <- function(x) is.na(x$pass) || isTRUE(x$pass)
  all_pass_fun <- function(cr) {
    criterion_pass(cr$reliability) && criterion_pass(cr$paths) &&
      criterion_pass(cr$fit) && criterion_pass(cr$cmb)
  }

  say("\n[Step 1] Evaluate initial data (n = ", initial_n, ")\n", sep = "")
  check_results <- check_all_criteria2(current_data, model_configs, var_list,
                                       criteria, verbose = verbose)

  if (all_pass_fun(check_results)) {
    say("\n✓ Initial data already meet all criteria — no removal needed\n")
    current_data$..row_id <- NULL
    return(list(
      cleaned_data = current_data, iterations = 0, initial_n = initial_n,
      final_n = nrow(current_data), total_removed = 0,
      removed_row_ids = removed_row_ids, all_pass = TRUE,
      final_gap = 0, reverted_to_iteration = NA_integer_,
      final_results = check_results,
      history = data.frame()
    ))
  }
  cur_gap <- total_gap2(check_results, criteria)
  say("\n✗ Criteria not met (aggregate gap = ", sprintf("%.4f", cur_gap),
      "); starting iterative removal\n", sep = "")

  # Best-state checkpoint: with multiple simultaneous criteria a removal
  # round can improve one criterion while worsening another; if the loop
  # ends without full success we restore the best intermediate state instead
  # of the (possibly worse) final one. State 0 is the untouched data, so the
  # function never returns something worse than its input.
  best <- list(gap = cur_gap, data = current_data, removed = removed_row_ids,
               results = check_results, iteration = 0)

  all_pass <- FALSE
  while (iteration < max_iterations) {
    iteration <- iteration + 1
    say("\n========== [Iteration ", iteration, "] removal rate ",
        round(removal_proportion * 100, 1), "% ==========\n", sep = "")

    # Removal quota: cumulative cap AND hard sample-size floor
    total_removed_so_far <- length(removed_row_ids)
    remaining_allow <- min(
      floor(max_total_remove_pct * initial_n) - total_removed_so_far,
      nrow(current_data) - min_n
    )
    if (remaining_allow <= 0) {
      say("⚠ Removal quota exhausted (cap ", max_total_remove_pct * 100,
          "% or min_n = ", min_n, "), stopping\n", sep = "")
      break
    }
    # 1e-9 tolerance guards against FP accumulation in removal_proportion
    # (e.g., 0.05 * 3 = 0.15000000000000002 would otherwise over-ceiling)
    n_quota <- if (pct_basis == "initial") {
      # Cumulative target relative to initial n: remove only the shortfall
      ceiling(removal_proportion * initial_n - 1e-9) - total_removed_so_far
    } else {
      # Legacy: proportion of the current sample each round
      ceiling(removal_proportion * nrow(current_data) - 1e-9)
    }
    if (n_quota <= 0) {
      # Cumulative target already met (e.g., a previous round hit it exactly);
      # raise the target and try again next iteration
      say("Cumulative removal target already met; raising target\n")
      removal_proportion <- removal_proportion + pct_increment
      next
    }
    n_to_remove <- min(n_quota, remaining_allow)
    say("Target removals this round: ", n_to_remove,
        if (pct_basis == "initial") {
          paste0(" (cumulative target ", round(removal_proportion * 100, 1),
                 "% of initial n)")
        } else "",
        "\n", sep = "")

    bad_idx <- tryCatch(
      identify_bad_cases2(current_data, model_configs, var_list, criteria,
                          n_to_remove, check_results,
                          verbose = verbose),
      error = function(e) {
        say("⚠ Case identification failed: ", conditionMessage(e), "\n", sep = "")
        integer(0)
      }
    )
    if (length(bad_idx) == 0) {
      say("⚠ No cases identified, stopping\n")
      break
    }

    removed_row_ids <- c(removed_row_ids, current_data$..row_id[bad_idx])
    current_data <- current_data[-bad_idx, , drop = FALSE]
    say("Removed ", length(bad_idx), " case(s); n = ", nrow(current_data),
        "\n", sep = "")

    say("\n[Re-evaluate]\n")
    check_results <- check_all_criteria2(current_data, model_configs, var_list,
                                         criteria, verbose = verbose)
    all_pass <- all_pass_fun(check_results)
    cur_gap <- total_gap2(check_results, criteria)

    flag <- function(x) if (is.na(x$pass)) "-" else if (isTRUE(x$pass)) "✓" else "✗"
    say(sprintf("\nStatus: Reliability[%s] Paths[%s] Fit[%s] CMB[%s] | gap = %.4f (best %.4f)\n",
                flag(check_results$reliability), flag(check_results$paths),
                flag(check_results$fit), flag(check_results$cmb),
                cur_gap, min(best$gap, cur_gap)))

    if (cur_gap < best$gap) {
      best <- list(gap = cur_gap, data = current_data,
                   removed = removed_row_ids, results = check_results,
                   iteration = iteration)
    }

    history[[iteration]] <- data.frame(
      iteration = iteration,
      n_removed = length(bad_idx), n_after = nrow(current_data),
      gap = cur_gap,
      reliability_pass = check_results$reliability$pass,
      paths_pass = check_results$paths$pass,
      fit_pass = check_results$fit$pass,
      cmb_pass = check_results$cmb$pass
    )

    if (all_pass) {
      say("\n✓ All criteria satisfied\n")
      break
    }
    removal_proportion <- removal_proportion + pct_increment
  }

  # Restore the best intermediate state if the loop ended on a worse one
  reverted_to <- NA_integer_
  if (!all_pass && best$gap < cur_gap) {
    say("\n↩ Final state (gap ", sprintf("%.4f", cur_gap),
        ") is worse than iteration ", best$iteration, " (gap ",
        sprintf("%.4f", best$gap), "); restoring that state\n", sep = "")
    current_data <- best$data
    removed_row_ids <- best$removed
    check_results <- best$results
    cur_gap <- best$gap
    reverted_to <- best$iteration
  }

  say("\n========================================\n")
  say("Finished: ", initial_n, " -> ", nrow(current_data), " (removed ",
      length(removed_row_ids), ", ",
      round(length(removed_row_ids) / initial_n * 100, 2), "%) in ",
      iteration, " iteration(s)\n", sep = "")
  if (!all_pass) {
    say("⚠ Stopped WITHOUT meeting all criteria (remaining gap = ",
        sprintf("%.4f", cur_gap), ")\n", sep = "")
  }

  current_data$..row_id <- NULL
  list(
    cleaned_data = current_data,
    iterations = iteration,
    initial_n = initial_n,
    final_n = nrow(current_data),
    total_removed = length(removed_row_ids),
    removed_row_ids = removed_row_ids,
    all_pass = all_pass,
    final_gap = cur_gap,
    reverted_to_iteration = reverted_to,
    final_results = check_results,
    history = if (length(history) > 0) do.call(rbind, history) else data.frame()
  )
}
