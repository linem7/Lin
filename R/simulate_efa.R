#' Simulate Likert Scale Data for Scale Development
#'
#' Generates simulated Likert scale data based on a structural model (lavaan syntax).
#' It allows specific items to be designated as "bad items" designed to be dropped
#' during Item Analysis (Stage 1) or Exploratory Factor Analysis (Stage 2).
#'
#' @param model_syntax A string containing the factor structure in lavaan syntax (e.g., 'f1 =~ x1 + x2').
#' @param items_to_drop A named list specifying items to be manipulated to fail at specific stages.
#'   Example: \code{list(stage1 = c("x1"), stage2 = c("x5"))}.
#' @param n_obs Integer. Sample size (default = 300).
#' @param likert_points Integer. Number of points on the Likert scale (default = 5).
#' @param skew_params A named list controlling skewness for realism.
#'   \code{good}: skewness for normal items (e.g., -0.5).
#'   \code{bad}: skewness for Stage 1 bad items (e.g., 2.0).
#'   \code{jitter}: random noise added to skewness.
#' @param ia_thresholds A named list of thresholds for Item Analysis logic checks.
#'   \code{citc}: Corrected Item-Total Correlation threshold (default 0.2).
#'   \code{p_val}: P-value threshold for Critical Ratio (CR) t-test (default 0.05).
#' @param efa_thresholds A named list of thresholds for EFA logic checks.
#'   \code{main_load}: Minimum main loading (default 0.4).
#'   \code{cross_load}: Maximum allowed cross-loading (default 0.35).
#' @param details Logical.
#'   If \code{FALSE} (default), prints only the summary logs of items recommended for deletion.
#'   If \code{TRUE}, prints detailed statistical tables (Item Analysis with stars) and detailed EFA results (using bruceR).
#'
#' @return A \code{data.frame} containing the simulated Likert scale data.
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
#' # Example 1: Unidimensional Model (1 factor, 10 items)
#' # x1, x2: Bad items for Stage 1 (Low CITC)
#' # x10: Bad item for Stage 2 (Low Loading)
#' # Note: We use 10 items to ensure enough remain for EFA after dropping bad ones.
#' model_uni <- 'F1 =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10'
#' bad_uni <- list(stage1 = c("x1", "x2"), stage2 = c("x10"))
#'
#' df_uni <- simulate_efa(
#'   model_syntax = model_uni,
#'   items_to_drop = bad_uni,
#'   n_obs = 300,
#'   details = TRUE
#' )
#'
#' # Example 2: Multidimensional Model (2 factors, 10 items)
#' # x1: Bad for Stage 1
#' # x6: Bad for Stage 2 (Cross-loading)
#' model_multi <- '
#'   F1 =~ x1 + x2 + x3 + x4 + x5
#'   F2 =~ x6 + x7 + x8 + x9 + x10
#' '
#' bad_multi <- list(stage1 = c("x1"), stage2 = c("x6"))
#'
#' df_multi <- simulate_efa(
#'   model_syntax = model_multi,
#'   items_to_drop = bad_multi,
#'   n_obs = 500,
#'   skew_params = list(good = -0.5, bad = 2.0, jitter = 0.2),
#'   details = FALSE # Will still print summary of drops
#' )
#' }
simulate_efa <- function(
    model_syntax,
    items_to_drop = list(),
    n_obs = 300,
    likert_points = 5,
    skew_params = list(good = -0.5, bad = 2.0, jitter = 0.2),
    ia_thresholds = list(citc = 0.2, p_val = 0.05),
    efa_thresholds = list(main_load = 0.4, cross_load = 0.35),
    details = FALSE
) {

  if(details) {
    cat("==========================================================\n")
    cat("          Simulation Data Generation Report               \n")
    cat("==========================================================\n\n")
  }

  # --- 1. Model Parsing ---
  ptable <- lavaan::lavaanify(model_syntax, auto = TRUE)
  factors <- unique(ptable$lhs[ptable$op == "=~"])
  items_map <- ptable[ptable$op == "=~", c("lhs", "rhs")]
  all_items <- unique(items_map$rhs)
  is_unidimensional <- length(factors) == 1

  # --- 2. Stage 2 Manipulation (Syntax Level for Multidimensional) ---
  model_for_sim <- model_syntax

  if (!is.null(items_to_drop$stage2)) {
    for (bad_item in items_to_drop$stage2) {
      if (!bad_item %in% all_items) next
      orig_factor <- items_map$lhs[items_map$rhs == bad_item]

      if (!is_unidimensional) {
        # Multidimensional: Create Cross-loading via Syntax
        # This is safe because it appends new lines
        distractor <- setdiff(factors, orig_factor)[1]
        cross_syntax <- paste0("\n", distractor, " =~ 0.5 * ", bad_item)
        model_for_sim <- paste0(model_for_sim, cross_syntax)
      }
      # Note: Unidimensional manipulation is handled at data level (Step 3.5)
      # to avoid regex issues that could break the model structure.
    }
  }

  # --- 3. Generate Base Data (Z-scores) ---
  # Fallback to standardized=FALSE if manipulation breaks positive definiteness
  raw_data <- tryCatch({
    lavaan::simulateData(model_for_sim, sample.nobs = n_obs, standardized = TRUE)
  }, error = function(e) {
    if(details) message("  [System] Standardized generation failed. Switching to unstandardized + manual scaling.")
    d <- lavaan::simulateData(model_for_sim, sample.nobs = n_obs, standardized = FALSE)
    as.data.frame(scale(d))
  })

  # --- 3.5 Stage 1 & 2 Data Manipulation (Signal Dilution) ---

  # A. Stage 1 Bad Items (Dilute Signal -> Low CITC)
  if (!is.null(items_to_drop$stage1)) {
    for (bad_item in items_to_drop$stage1) {
      if(bad_item %in% colnames(raw_data)) {
        original_signal <- raw_data[[bad_item]]
        noise <- stats::rnorm(n_obs, mean = 0, sd = 1)
        # Mix: 15% Signal + 85% Noise
        mixed_signal <- 0.15 * original_signal + 0.85 * noise
        raw_data[[bad_item]] <- scale(mixed_signal)[,1]
      }
    }
  }

  # B. Stage 2 Bad Items (Unidimensional ONLY -> Low Loading)
  # This avoids modifying the model syntax for single-factor models, which is error-prone.
  if (is_unidimensional && !is.null(items_to_drop$stage2)) {
    for (bad_item in items_to_drop$stage2) {
      if(bad_item %in% colnames(raw_data)) {
        if(details) message(paste("  [Info] Diluting item", bad_item, "to create low loading (Unidimensional)."))
        original_signal <- raw_data[[bad_item]]
        noise <- stats::rnorm(n_obs)
        # Mix: 30% Signal + 70% Noise -> Loading approx 0.3
        mixed_signal <- 0.30 * original_signal + 0.70 * noise
        raw_data[[bad_item]] <- scale(mixed_signal)[,1]
      }
    }
  }

  # --- 4. Likert Conversion (Applying Skew) ---
  likert_data <- raw_data
  for (col in colnames(likert_data)) {
    if (!col %in% all_items) next
    sk_val <- 0
    if (col %in% items_to_drop$stage1) {
      # Bad items: extreme skew
      sk_val <- skew_params$bad * sample(c(1, -1), 1)
    } else {
      # Good items: specified skew + jitter
      sk_val <- skew_params$good + stats::runif(1, -skew_params$jitter, skew_params$jitter)
    }
    likert_data[[col]] <- .sim_likertize(raw_data[[col]], n_levels = likert_points, skew = sk_val)
  }

  # =======================================================
  # Stage 1: Item Analysis
  # =======================================================
  ia_log <- data.frame(Item = character(), Reason = character(), stringsAsFactors = FALSE)
  ia_drop_list <- c()

  if(details) {
    cat("\n----------------------------------------------------------\n")
    cat(" [Stage 1] Detailed Item Analysis Report\n")
    cat("----------------------------------------------------------\n")
  }

  for (f in factors) {
    f_items <- intersect(items_map$rhs[items_map$lhs == f], colnames(likert_data))
    if (length(f_items) < 2) next

    # Run Internal Helper
    ia_res <- .sim_run_item_analysis(likert_data, pattern = "{i}", indices = f_items)
    if (is.null(ia_res)) next

    # 1. Print Detailed Table ONLY if details = TRUE
    if(details) {
      cat(paste0("\nFactor: ", f, "\n"))
      print(ia_res$formatted, row.names = FALSE)
    }

    # 2. Logic Check (Always runs)
    raw_stats <- ia_res$raw
    for (i in seq_len(nrow(raw_stats))) {
      item_name <- raw_stats$Item[i]
      reasons <- c()

      if (raw_stats$Val_CR_p[i] > ia_thresholds$p_val) reasons <- c(reasons, "CR Not Sig")
      if (raw_stats$Val_CITC[i] < ia_thresholds$citc) reasons <- c(reasons, paste0("CITC < ", ia_thresholds$citc))
      if (raw_stats$Val_Alpha_Del[i] > raw_stats$Val_Total_Alpha[i] + 0.002) reasons <- c(reasons, "Alpha Increases")

      # Safety check for zero variance
      if (raw_stats$Is_Zero_Var[i]) reasons <- c(reasons, "Zero Variance")

      is_target_bad <- item_name %in% items_to_drop$stage1
      if (length(reasons) > 0) {
        # Drop if target OR stats are extremely poor OR Zero Variance
        if (is_target_bad || (raw_stats$Val_CITC[i] < ia_thresholds$citc) || raw_stats$Is_Zero_Var[i]) {
          ia_drop_list <- c(ia_drop_list, item_name)
          ia_log <- rbind(ia_log, data.frame(Item = item_name, Reason = paste(reasons, collapse = "; ")))
        }
      }
    }
  }

  # [Output] Always print Summary Log
  cat("\n--- [Stage 1] Log: Items recommended for deletion ---\n")
  if (nrow(ia_log) > 0) print(ia_log, row.names = FALSE) else cat("  None.\n")

  # =======================================================
  # Stage 2: EFA Check
  # =======================================================
  ia_drop_unique <- unique(ia_drop_list)
  # Use any_of() to handle cases where drop list might be empty
  data_for_efa <- likert_data %>% dplyr::select(-dplyr::any_of(ia_drop_unique))

  # Clean zero variance
  sds <- sapply(data_for_efa, sd, na.rm=TRUE)
  data_for_efa <- data_for_efa[, sds > 1e-9, drop = FALSE]

  efa_log <- data.frame(Item = character(), Reason = character(), stringsAsFactors = FALSE)

  # Check if we have enough items for EFA (at least 3 recommended)
  can_run_efa <- ncol(data_for_efa) >= 3

  # 1. Detailed EFA Output (bruceR) - Only if details=TRUE
  if(details) {
    cat("\n----------------------------------------------------------\n")
    cat(" [Stage 2] Detailed EFA Report (Cleaned Data)\n")
    cat("----------------------------------------------------------\n")

    if (can_run_efa && requireNamespace("bruceR", quietly = TRUE)) {
      cat(paste0("Running EFA on ", ncol(data_for_efa), " items...\n"))
      tryCatch({
        bruceR::EFA(data_for_efa, vars = colnames(data_for_efa), sort.loadings = FALSE,
                    hide.loadings = 0.3)
      }, error = function(e) {
        cat("Error running bruceR::EFA: ", e$message, "\n")
      })
    } else {
      if(!can_run_efa) cat("Skipping detailed EFA (Too few items).\n")
      if(!requireNamespace("bruceR", quietly = TRUE)) cat("Package 'bruceR' is not installed.\n")
    }
  }

  # 2. Internal EFA Logic Check (Always runs)
  if (can_run_efa) {
    tryCatch({
      efa_res <- psych::fa(data_for_efa, nfactors = length(factors), rotate = "promax", fm = "pa", warnings = FALSE)
      loadings <- unclass(efa_res$loadings)

      for (itm in rownames(loadings)) {
        abs_loads <- abs(loadings[itm, ])
        sorted_loads <- sort(abs_loads, decreasing = TRUE)
        main_load <- sorted_loads[1]
        sec_load  <- if(length(sorted_loads) > 1) sorted_loads[2] else 0

        reasons <- c()
        if (main_load < efa_thresholds$main_load) {
          reasons <- c(reasons, paste0("Main Load < ", efa_thresholds$main_load))
        }
        if (!is_unidimensional && sec_load > efa_thresholds$cross_load) {
          reasons <- c(reasons, paste0("Cross Load > ", efa_thresholds$cross_load))
        }

        if (length(reasons) > 0) {
          efa_log <- rbind(efa_log, data.frame(Item = itm, Reason = paste(reasons, collapse = "; ")))
        }
      }
    }, error = function(e) {
      # Silent error for internal check
    })
  }

  # [Output] Always print Summary Log
  cat("\n--- [Stage 2] Internal Check Log ---\n")
  if (!can_run_efa) {
    cat("  Skipping EFA check (Too few items remaining).\n")
  } else if (nrow(efa_log) > 0) {
    print(efa_log, row.names = FALSE)
  } else {
    cat("  None (Structure Clean).\n")
  }

  if(details) {
    cat("\n==========================================================\n")
    cat(paste0("Process Complete. Returning dataframe (", nrow(likert_data), " x ", ncol(likert_data), ").\n"))
    cat("==========================================================\n")
  }

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
