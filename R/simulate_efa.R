#' Simulate Likert Scale Data for Scale Development
#'
#' Generates simulated Likert scale data based on a structural model (lavaan syntax).
#' It allows specific items to be designated as "bad items" designed to be dropped
#' during Item Analysis (Stage 1) or Exploratory Factor Analysis (Stage 2).
#'
#' The function provides control over factor loadings (\code{loading_range}),
#' inter-factor correlations (\code{factor_cor}), and item-level skewness
#' (\code{skew_params}) to produce data that resembles realistic survey responses.
#'
#' @param model_syntax A string containing the factor structure in lavaan syntax
#'   (e.g., \code{'F1 =~ x1 + x2 + x3'}). Each factor definition should occupy
#'   its own line. Items that already contain explicit loading specifications
#'   (e.g., \code{0.8*x1}) will not be overwritten by \code{loading_range}.
#' @param items_to_drop A named list specifying items to be manipulated to fail
#'   at specific stages.
#'   \describe{
#'     \item{\code{stage1}}{Character vector. Items whose signals will be diluted
#'       (15\% signal + 85\% noise) to produce low Corrected Item-Total
#'       Correlations (CITC), making them candidates for deletion in Item Analysis.}
#'     \item{\code{stage2}}{Character vector. In multidimensional models, these
#'       items receive a cross-loading (0.5) on a non-target factor. In
#'       unidimensional models, their signals are diluted (30\% signal + 70\%
#'       noise) to produce low factor loadings in EFA.}
#'   }
#' @param n_obs Integer. Sample size (default = 300).
#' @param likert_points Integer. Number of points on the Likert scale (default = 5).
#' @param factor_cor Controls the correlations among latent factors.
#'   Ignored when the model contains only one factor. Accepts three forms:
#'   \describe{
#'     \item{Single numeric}{A scalar (default = 0.3) applied to every pair of
#'       factors. For example, \code{factor_cor = 0.5} sets all pairwise
#'       correlations to 0.5.}
#'     \item{Symmetric matrix}{A \eqn{k \times k} correlation matrix where
#'       \eqn{k} equals the number of factors. Row and column names, if present,
#'       must match the factor names in \code{model_syntax}. Diagonal elements
#'       should be 1.}
#'     \item{\code{NULL}}{No factor covariance syntax is appended. lavaan will
#'       use its own defaults (typically free estimation, which under
#'       \code{standardized = TRUE} may default to zero).}
#'   }
#' @param loading_range A numeric vector of length 2 specifying the lower and
#'   upper bounds for randomly generated factor loadings (default =
#'   \code{c(0.70, 0.85)}). Each item without an explicit loading in
#'   \code{model_syntax} receives a loading drawn from
#'   \code{Uniform(loading_range[1], loading_range[2])}. Set to \code{NULL}
#'   to skip loading injection and use lavaan defaults.
#' @param skew_params A named list controlling skewness for realism.
#'   \describe{
#'     \item{\code{good}}{Skewness applied to normal items (default = -0.5).}
#'     \item{\code{bad}}{Absolute skewness applied to Stage 1 bad items, with
#'       sign randomly flipped (default = 2.0).}
#'     \item{\code{jitter}}{Maximum random deviation added to the skewness of
#'       normal items (default = 0.2).}
#'   }
#' @param ia_thresholds A named list of thresholds for Item Analysis.
#'   \describe{
#'     \item{\code{citc}}{Corrected Item-Total Correlation threshold (default = 0.2).}
#'     \item{\code{p_val}}{P-value threshold for the Critical Ratio t-test
#'       (default = 0.05).}
#'   }
#' @param efa_thresholds A named list of thresholds for the internal EFA check.
#'   \describe{
#'     \item{\code{main_load}}{Minimum acceptable main loading (default = 0.4).}
#'     \item{\code{cross_load}}{Maximum acceptable cross-loading (default = 0.35).
#'       Only evaluated in multidimensional models.}
#'   }
#' @param details Logical.
#'   If \code{FALSE} (default), prints only the summary logs of items
#'   recommended for deletion at each stage.
#'   If \code{TRUE}, prints the injected model syntax, factor correlation
#'   syntax, detailed Item Analysis tables (with significance stars), and
#'   detailed EFA results (via \pkg{bruceR} if available).
#'
#' @return A \code{data.frame} containing the simulated Likert scale data.
#'   All columns are numeric. The data frame includes every item specified in
#'   \code{model_syntax} (including bad items), so that the user can reproduce
#'   the full item analysis and deletion workflow.
#'
#' @section Workflow:
#' The function proceeds through the following stages internally:
#' \enumerate{
#'   \item \strong{Model construction.} Parse \code{model_syntax}, inject random
#'     loadings (if \code{loading_range} is specified), and append factor
#'     covariance syntax (if \code{factor_cor} is specified).
#'   \item \strong{Data generation.} Call \code{lavaan::simulateData()} with the
#'     assembled model. If \code{standardized = TRUE} fails (e.g., due to a
#'     non-positive-definite implied matrix), the function falls back to
#'     \code{standardized = FALSE} with post-hoc standardization.
#'   \item \strong{Signal dilution.} Manipulate designated bad items at the
#'     continuous (z-score) level before Likert conversion.
#'   \item \strong{Likert conversion.} Map z-scores to discrete Likert values
#'     using skewed quantile cutpoints.
#'   \item \strong{Item Analysis (Stage 1).} Evaluate CITC, Critical Ratio,
#'     and Cronbach's alpha if deleted for each factor's items.
#'   \item \strong{EFA check (Stage 2).} Run an internal EFA
#'     (\code{psych::fa()}) on the items surviving Stage 1 and flag those with
#'     low main loadings or high cross-loadings.
#' }
#'
#' @importFrom lavaan lavaanify simulateData
#' @importFrom psych alpha fa describe
#' @importFrom dplyr select mutate filter case_when pull left_join all_of across where any_of
#' @importFrom stringr str_detect str_extract str_replace
#' @importFrom glue glue
#' @importFrom sn qsn
#' @importFrom tibble rownames_to_column
#' @importFrom stats t.test cor cor.test quantile sd var runif rnorm
#' @export
#'
#' @examples
#' \dontrun{
#' # --- Example 1: Unidimensional Model ---
#' # 1 factor, 10 items.
#' # x1, x2 are designed to fail Item Analysis (Stage 1).
#' # x10 is designed to show a low loading in EFA (Stage 2).
#' model_uni <- 'F1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10'
#' bad_uni <- list(stage1 = c("x1", "x2"), stage2 = c("x10"))
#'
#' df_uni <- simulate_efa(
#'   model_syntax = model_uni,
#'   items_to_drop = bad_uni,
#'   n_obs = 300,
#'   loading_range = c(0.65, 0.80),
#'   details = TRUE
#' )
#'
#' # --- Example 2: Multidimensional Model with Scalar factor_cor ---
#' # 3 factors, 15 items. All pairwise factor correlations set to 0.4.
#' # x1 fails Stage 1; x6 cross-loads in Stage 2.
#' model_3f <- '
#'   F1 =~ x1 + x2 + x3 + x4 + x5
#'   F2 =~ x6 + x7 + x8 + x9 + x10
#'   F3 =~ x11 + x12 + x13 + x14 + x15
#' '
#' bad_3f <- list(stage1 = c("x1"), stage2 = c("x6"))
#'
#' df_3f <- simulate_efa(
#'   model_syntax = model_3f,
#'   items_to_drop = bad_3f,
#'   n_obs = 500,
#'   factor_cor = 0.4,
#'   details = FALSE
#' )
#'
#' # --- Example 3: Multidimensional Model with Correlation Matrix ---
#' # Specify different correlations for each factor pair.
#' cor_mat <- matrix(c(
#'   1.0, 0.5, 0.2,
#'   0.5, 1.0, 0.6,
#'   0.2, 0.6, 1.0
#' ), nrow = 3, byrow = TRUE,
#'    dimnames = list(c("F1","F2","F3"), c("F1","F2","F3")))
#'
#' df_custom <- simulate_efa(
#'   model_syntax = model_3f,
#'   items_to_drop = bad_3f,
#'   n_obs = 400,
#'   factor_cor = cor_mat,
#'   loading_range = c(0.70, 0.90),
#'   skew_params = list(good = -0.3, bad = 1.5, jitter = 0.1),
#'   details = TRUE
#' )
#' }
simulate_efa <- function(
    model_syntax,
    items_to_drop = list(),
    n_obs = 300,
    likert_points = 5,
    factor_cor = 0.3,
    loading_range = c(0.70, 0.85),
    skew_params = list(good = -0.5, bad = 2.0, jitter = 0.2),
    ia_thresholds = list(citc = 0.2, p_val = 0.05),
    efa_thresholds = list(main_load = 0.4, cross_load = 0.32),
    details = FALSE
) {

  if (details) {
    cat("==========================================================\n")
    cat("          Simulation Data Generation Report               \n")
    cat("==========================================================\n\n")
  }

  # --- 0. Input Validation ---
  if (!is.null(loading_range)) {
    if (!is.numeric(loading_range) || length(loading_range) != 2) {
      stop("loading_range must be a numeric vector of length 2 (e.g., c(0.70, 0.85)).")
    }
    if (loading_range[1] > loading_range[2]) {
      loading_range <- sort(loading_range)
    }
  }

  # --- 1. Model Parsing ---
  ptable <- lavaan::lavaanify(model_syntax, auto = TRUE)
  factors <- unique(ptable$lhs[ptable$op == "=~"])
  items_map <- ptable[ptable$op == "=~", c("lhs", "rhs")]
  all_items <- unique(items_map$rhs)
  is_unidimensional <- length(factors) == 1
  n_factors <- length(factors)

  # --- 1.3 Inject Factor Loadings ---
  if (!is.null(loading_range)) {
    lines <- strsplit(model_syntax, "\n")[[1]]
    new_lines <- c()

    for (line in lines) {
      trimmed <- trimws(line)
      if (!grepl("=~", trimmed)) {
        new_lines <- c(new_lines, line)
        next
      }

      parts <- strsplit(trimmed, "=~")[[1]]
      lhs <- trimws(parts[1])
      rhs_items <- trimws(strsplit(parts[2], "\\+")[[1]])

      new_rhs <- sapply(rhs_items, function(item) {
        clean_item <- trimws(item)
        if (grepl("\\*", clean_item)) return(clean_item)
        lam <- round(stats::runif(1, loading_range[1], loading_range[2]), 3)
        paste0(lam, "*", clean_item)
      })

      new_lines <- c(new_lines, paste0(lhs, " =~ ", paste(new_rhs, collapse = " + ")))
    }

    model_syntax <- paste(new_lines, collapse = "\n")

    if (details) {
      cat("[Info] Model with injected loadings:\n")
      cat(model_syntax, "\n\n")
    }
  }

  # Re-parse after injection
  ptable <- lavaan::lavaanify(model_syntax, auto = TRUE)
  factors <- unique(ptable$lhs[ptable$op == "=~"])
  items_map <- ptable[ptable$op == "=~", c("lhs", "rhs")]
  all_items <- unique(items_map$rhs)

  # --- 1.5 Build Factor Covariance Syntax ---
  factor_cov_syntax <- ""

  if (!is_unidimensional && !is.null(factor_cor) && n_factors >= 2) {

    if (is.matrix(factor_cor)) {
      cor_mat <- factor_cor
      if (nrow(cor_mat) != n_factors || ncol(cor_mat) != n_factors) {
        stop(paste0("factor_cor matrix dimensions (", nrow(cor_mat), "x", ncol(cor_mat),
                    ") do not match the number of factors (", n_factors, ")."))
      }
      if (!is.null(rownames(cor_mat)) && !is.null(colnames(cor_mat))) {
        missing <- setdiff(factors, rownames(cor_mat))
        if (length(missing) > 0) {
          stop(paste0("factor_cor matrix is missing factor name(s): ",
                      paste(missing, collapse = ", ")))
        }
        cor_mat <- cor_mat[factors, factors]
      }
    } else if (is.numeric(factor_cor) && length(factor_cor) == 1) {
      cor_mat <- matrix(factor_cor, nrow = n_factors, ncol = n_factors)
      diag(cor_mat) <- 1.0
    } else {
      stop("factor_cor must be either a single numeric value or a square matrix.")
    }

    cov_lines <- c()
    for (i in 2:n_factors) {
      for (j in 1:(i - 1)) {
        r <- cor_mat[i, j]
        cov_lines <- c(cov_lines,
                       paste0(factors[i], " ~~ ", round(r, 4), " * ", factors[j]))
      }
    }
    factor_cov_syntax <- paste(cov_lines, collapse = "\n")

    if (details) {
      cat("[Info] Factor correlation syntax appended:\n")
      cat(factor_cov_syntax, "\n\n")
    }
  }

  # --- 2. Stage 2 Manipulation (Syntax Level for Multidimensional) ---
  model_for_sim <- model_syntax

  if (nchar(factor_cov_syntax) > 0) {
    model_for_sim <- paste0(model_for_sim, "\n", factor_cov_syntax)
  }

  if (!is.null(items_to_drop$stage2)) {
    for (bad_item in items_to_drop$stage2) {
      if (!bad_item %in% all_items) next
      orig_factor <- items_map$lhs[items_map$rhs == bad_item]

      if (!is_unidimensional) {
        distractor <- setdiff(factors, orig_factor)[1]
        cross_syntax <- paste0("\n", distractor, " =~ 0.5 * ", bad_item)
        model_for_sim <- paste0(model_for_sim, cross_syntax)
      }
    }
  }

  if (details) {
    cat("[Info] Final model for simulation:\n")
    cat(model_for_sim, "\n\n")
  }

  # --- 3. Generate Base Data (Z-scores) ---
  raw_data <- tryCatch({
    lavaan::simulateData(model_for_sim, sample.nobs = n_obs, standardized = TRUE)
  }, error = function(e) {
    cat("  [System] Standardized generation failed (", conditionMessage(e),
        "). Switching to unstandardized + manual scaling.\n")
    d <- lavaan::simulateData(model_for_sim, sample.nobs = n_obs, standardized = FALSE)
    as.data.frame(scale(d))
  })

  # --- 3.5 Stage 1 & 2 Data Manipulation (Signal Dilution) ---

  # A. Stage 1 Bad Items (Dilute Signal -> Low CITC)
  if (!is.null(items_to_drop$stage1)) {
    for (bad_item in items_to_drop$stage1) {
      if (bad_item %in% colnames(raw_data)) {
        original_signal <- raw_data[[bad_item]]
        noise <- stats::rnorm(n_obs, mean = 0, sd = 1)
        mixed_signal <- 0.15 * original_signal + 0.85 * noise
        raw_data[[bad_item]] <- scale(mixed_signal)[, 1]
      }
    }
  }

  # B. Stage 2 Bad Items (Unidimensional ONLY -> Low Loading)
  if (is_unidimensional && !is.null(items_to_drop$stage2)) {
    for (bad_item in items_to_drop$stage2) {
      if (bad_item %in% colnames(raw_data)) {
        if (details) cat("  [Info] Diluting item", bad_item,
                         "to create low loading (unidimensional model).\n")
        original_signal <- raw_data[[bad_item]]
        noise <- stats::rnorm(n_obs)
        mixed_signal <- 0.30 * original_signal + 0.70 * noise
        raw_data[[bad_item]] <- scale(mixed_signal)[, 1]
      }
    }
  }

  # --- 4. Likert Conversion (Applying Skew) ---
  likert_data <- raw_data
  for (col in colnames(likert_data)) {
    if (!col %in% all_items) next
    sk_val <- 0
    if (col %in% items_to_drop$stage1) {
      sk_val <- skew_params$bad * sample(c(1, -1), 1)
    } else {
      sk_val <- skew_params$good + stats::runif(1, -skew_params$jitter, skew_params$jitter)
    }
    likert_data[[col]] <- .sim_likertize(raw_data[[col]], n_levels = likert_points, skew = sk_val)
  }

  # =======================================================
  # Stage 1: Item Analysis
  # =======================================================
  ia_log <- data.frame(Item = character(), Reason = character(), stringsAsFactors = FALSE)
  ia_drop_list <- c()

  if (details) {
    cat("\n----------------------------------------------------------\n")
    cat(" [Stage 1] Detailed Item Analysis Report\n")
    cat("----------------------------------------------------------\n")
  }

  for (f in factors) {
    f_items <- intersect(items_map$rhs[items_map$lhs == f], colnames(likert_data))
    if (length(f_items) < 2) next

    ia_res <- .sim_run_item_analysis(likert_data, pattern = "{i}", indices = f_items)
    if (is.null(ia_res)) next

    if (details) {
      cat(paste0("\nFactor: ", f, "\n"))
      print(ia_res$formatted, row.names = FALSE)
    }

    raw_stats <- ia_res$raw
    for (i in seq_len(nrow(raw_stats))) {
      item_name <- raw_stats$Item[i]
      reasons <- c()

      if (raw_stats$Val_CR_p[i] > ia_thresholds$p_val) reasons <- c(reasons, "CR Not Sig")
      if (raw_stats$Val_CITC[i] < ia_thresholds$citc) reasons <- c(reasons, paste0("CITC < ", ia_thresholds$citc))
      if (raw_stats$Val_Alpha_Del[i] > raw_stats$Val_Total_Alpha[i] + 0.002) reasons <- c(reasons, "Alpha Increases")
      if (raw_stats$Is_Zero_Var[i]) reasons <- c(reasons, "Zero Variance")

      is_target_bad <- item_name %in% items_to_drop$stage1
      if (length(reasons) > 0) {
        if (is_target_bad || (raw_stats$Val_CITC[i] < ia_thresholds$citc) || raw_stats$Is_Zero_Var[i]) {
          ia_drop_list <- c(ia_drop_list, item_name)
          ia_log <- rbind(ia_log, data.frame(Item = item_name, Reason = paste(reasons, collapse = "; ")))
        }
      }
    }
  }

  # [Stage 1 Log] — message() for highlighted display
  message("\n--- [Stage 1] Log: Items recommended for deletion ---")
  if (nrow(ia_log) > 0) {
    msg_lines <- utils::capture.output(print(ia_log, row.names = FALSE))
    message(paste(msg_lines, collapse = "\n"))
  } else {
    message("  None.")
  }

  # =======================================================
  # Stage 2: EFA Check
  # =======================================================
  ia_drop_unique <- unique(ia_drop_list)
  data_for_efa <- likert_data %>% dplyr::select(-dplyr::any_of(ia_drop_unique))

  sds <- sapply(data_for_efa, sd, na.rm = TRUE)
  data_for_efa <- data_for_efa[, sds > 1e-9, drop = FALSE]

  efa_log <- data.frame(Item = character(), Reason = character(), stringsAsFactors = FALSE)
  can_run_efa <- ncol(data_for_efa) >= 3

  # Determine effective number of factors for EFA
  surviving_items <- colnames(data_for_efa)
  factors_with_items <- unique(items_map$lhs[items_map$rhs %in% surviving_items])
  n_efa_factors <- max(1L, length(factors_with_items))

  if (n_efa_factors < length(factors) && can_run_efa) {
    cat("  [Warning] Only ", n_efa_factors, " of ", length(factors),
        " factors have surviving items. EFA will extract ", n_efa_factors, " factor(s).\n")
  }

  # Detailed EFA report
  if (details) {
    cat("\n----------------------------------------------------------\n")
    cat(" [Stage 2] Detailed EFA Report (Cleaned Data)\n")
    cat("----------------------------------------------------------\n")

    if (can_run_efa && requireNamespace("bruceR", quietly = TRUE)) {
      cat("Running EFA on ", ncol(data_for_efa), " items with ",
          n_efa_factors, " factor(s)...\n")
      tryCatch({
        bruceR::EFA(data_for_efa, vars = colnames(data_for_efa),
                    sort.loadings = FALSE, hide.loadings = efa_thresholds["cross_load"]$cross_load)
      }, error = function(e) {
        cat("  [Error] bruceR::EFA failed: ", conditionMessage(e), "\n")
      })
    } else {
      if (!can_run_efa) cat("Skipping detailed EFA (too few items).\n")
      if (!requireNamespace("bruceR", quietly = TRUE)) cat("Package 'bruceR' is not installed.\n")
    }
  }

  # Internal EFA Logic Check (Always runs)
  efa_check_ran <- FALSE

  if (can_run_efa) {
    tryCatch({
      pca_res <- psych::principal(data_for_efa, nfactors = n_efa_factors,
                                  rotate = "varimax")
      loadings <- unclass(pca_res$loadings)

      for (itm in rownames(loadings)) {
        abs_loads <- abs(loadings[itm, ])
        sorted_loads <- sort(abs_loads, decreasing = TRUE)
        main_load <- sorted_loads[1]
        sec_load  <- if (length(sorted_loads) > 1) sorted_loads[2] else 0

        reasons <- c()
        if (main_load < efa_thresholds$main_load) {
          reasons <- c(reasons, paste0("Main Load < ", efa_thresholds$main_load))
        }
        if (n_efa_factors > 1 && sec_load > efa_thresholds$cross_load) {
          reasons <- c(reasons, paste0("Cross Load > ", efa_thresholds$cross_load))
        }

        if (length(reasons) > 0) {
          efa_log <- rbind(efa_log, data.frame(Item = itm,
                                               Reason = paste(reasons, collapse = "; ")))
        }
      }

      efa_check_ran <- TRUE
      if (details) cat("  [Info] Internal PCA check completed.\n")

    }, error = function(e) {
      if (details) cat("  [Info] PCA check failed: ", conditionMessage(e), "\n")
    })

    if (!efa_check_ran) {
      cat("  [Warning] Internal PCA check failed.\n")
    }
  }

  # [Stage 2 Log]
  message("\n--- [Stage 2] Internal Check Log ---")
  if (!can_run_efa) {
    message("  Skipping EFA check (too few items remaining).")
  } else if (!efa_check_ran) {
    message("  PCA extraction failed. Unable to evaluate factor structure.")
  } else if (nrow(efa_log) > 0) {
    msg_lines <- utils::capture.output(print(efa_log, row.names = FALSE))
    message(paste(msg_lines, collapse = "\n"))
  } else {
    message("  None (Structure Clean).")
  }


  if (details) {
    cat("\n==========================================================\n")
    cat("Process complete. Returning data frame (", nrow(likert_data),
        " x ", ncol(likert_data), ").\n")
    cat("==========================================================\n")
  }

  likert_data[] <- lapply(likert_data, as.numeric)

  return(likert_data)
}


# ==============================================================================
# Internal Helper Functions (Not exported, @noRd)
# ==============================================================================

#' Internal: Convert Z-scores to Likert scale with skewness
#' @noRd
.sim_likertize <- function(z_score, n_levels = 5, skew = 0) {
  if (!is.numeric(z_score)) return(z_score)

  # Cap skewness to prevent numerical overflow in sn::qsn
  safe_skew <- sign(skew) * min(abs(skew), 8)
  probs <- c(0.001, 0.999)

  if (abs(safe_skew) < 1e-4) {
    range_lims <- stats::qnorm(probs)
  } else {
    range_lims <- sn::qsn(probs, xi = 0, omega = 1, alpha = safe_skew)
  }

  endpoints <- seq(range_lims[1], range_lims[2], length.out = n_levels + 1)
  breaks <- endpoints
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf

  likert <- cut(z_score, breaks = breaks, labels = 1:n_levels, include.lowest = TRUE)
  res <- as.integer(as.character(likert))
  if(all(is.na(res))) res <- rep(as.integer(n_levels/2), length(z_score))
  return(res)
}

#' Internal: Robust Item Analysis with "Total" Column
#' @noRd
.sim_run_item_analysis <- function(data, pattern, indices, digits = 3) {
  item_names <- as.character(glue::glue(pattern, i = indices))
  valid_names <- intersect(item_names, colnames(data))
  if (length(valid_names) < 2) return(NULL)

  selected_items <- data %>% dplyr::select(dplyr::all_of(valid_names))

  # 1. Descriptives
  desc_df <- tryCatch({
    psych::describe(selected_items) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Item") %>%
      dplyr::select(Item, Mean = mean, SD = sd, Skew = skew, Kurt = kurtosis)
  }, error = function(e) {
    data.frame(Item = valid_names, Mean = NA, SD = NA, Skew = NA, Kurt = NA)
  })

  # 2. Alpha
  alpha_res <- tryCatch({
    suppressWarnings(psych::alpha(selected_items, warnings = FALSE, check.keys = TRUE))
  }, error = function(e) NULL)

  if (is.null(alpha_res)) {
    drop_vals <- rep(0, length(valid_names)); names(drop_vals) <- valid_names
    total_alpha <- 0
  } else {
    ad <- alpha_res$alpha.drop
    if (is.data.frame(ad)) {
      drop_vals <- if ("raw_alpha" %in% colnames(ad)) ad[["raw_alpha"]] else ad[[1]]
    } else {
      drop_vals <- unlist(ad)
    }
    names(drop_vals) <- rownames(ad)
    total_alpha <- alpha_res$total$raw_alpha
  }

  # 3. CR & CITC Calculation
  row_sums <- rowSums(selected_items, na.rm = TRUE)
  tmp <- selected_items %>%
    dplyr::mutate(
      TotalScore = row_sums,
      PerformanceGroup = dplyr::case_when(
        TotalScore <= stats::quantile(TotalScore, 0.27, na.rm = TRUE) ~ "Low",
        TotalScore >= stats::quantile(TotalScore, 0.73, na.rm = TRUE) ~ "High",
        TRUE ~ NA_character_
      )
    )

  stats_list <- lapply(valid_names, function(it) {
    # CR
    high <- tmp %>% dplyr::filter(PerformanceGroup == "High") %>% dplyr::pull(.data[[it]])
    low  <- tmp %>% dplyr::filter(PerformanceGroup == "Low")  %>% dplyr::pull(.data[[it]])
    t_val <- 0; p_val <- 1
    if (stats::var(high, na.rm=TRUE) > 1e-9 || stats::var(low, na.rm=TRUE) > 1e-9) {
      try({
        tt <- stats::t.test(high, low)
        t_val <- tt$statistic; p_val <- tt$p.value
      }, silent = TRUE)
    }

    # CITC
    others <- setdiff(valid_names, it)
    sum_oth <- rowSums(selected_items[, others, drop = FALSE], na.rm = TRUE)
    r_val <- 0; citc_p <- 1
    if (stats::sd(selected_items[[it]], na.rm=TRUE) > 1e-9 && stats::sd(sum_oth, na.rm=TRUE) > 1e-9) {
      try({
        ct <- stats::cor.test(selected_items[[it]], sum_oth)
        r_val <- ct$estimate; citc_p <- ct$p.value
      }, silent = TRUE)
    }

    stars_cr <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else if (p_val < 0.05) "*" else ""
    stars_citc <- if (citc_p < 0.001) "***" else if (citc_p < 0.01) "**" else if (citc_p < 0.05) "*" else ""

    list(
      Item = it,
      Val_CR_p = p_val, Val_CITC = r_val, Val_Alpha_Del = drop_vals[it], Is_Zero_Var = (stats::sd(selected_items[[it]], na.rm=TRUE) < 1e-9),
      Val_Total_Alpha = total_alpha,
      CR_Fmt = paste0(sprintf(paste0("%.", digits, "f"), t_val), stars_cr),
      CITC_Fmt = paste0(sprintf(paste0("%.", digits, "f"), r_val), stars_citc),
      Alpha_Del_Fmt = sprintf(paste0("%.", digits, "f"), drop_vals[it])
    )
  })

  stats_df <- do.call(rbind, lapply(stats_list, as.data.frame))

  # Format Table: Total as last COLUMN
  formatted_df <- desc_df %>%
    dplyr::left_join(stats_df, by = "Item") %>%
    dplyr::mutate(
      Mean = sprintf("%.2f", Mean), SD = sprintf("%.2f", SD),
      Skew = sprintf("%.2f", Skew), Kurt = sprintf("%.2f", Kurt)
    ) %>%
    dplyr::select(Item, Mean, SD, Skew, Kurt, CR = CR_Fmt, CITC = CITC_Fmt, `Alpha if del` = Alpha_Del_Fmt)

  # Add Total Alpha Column (Last Column)
  formatted_df$Total <- ""
  if(nrow(formatted_df) > 0) formatted_df$Total[1] <- sprintf(paste0("%.", digits, "f"), total_alpha)

  return(list(raw = stats_df, formatted = formatted_df))
}
