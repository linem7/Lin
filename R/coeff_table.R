#' Create formatted tables for Mplus models
#'
#' Extracts and organises results from an Mplus model
#' (read via \code{MplusAutomation::readModels()}) into up to three
#' publication-ready tables: model fit indices, path coefficients (with
#' MODEL CONSTRAINT parameters when present), and a decomposition of total,
#' indirect, and direct effects. Tables can be rendered as \pkg{kableExtra}
#' HTML tables or returned as plain \code{data.frame}s for downstream
#' processing.
#'
#' @param model A single \code{mplus.model} object returned by
#'   \code{MplusAutomation::readModels()}.
#' @param stat Character vector of statistics to display for each path.
#'   Any combination of \code{"est"} (point estimate), \code{"se"} (standard
#'   error), \code{"z"} (z-statistic, i.e.\ est/se), \code{"p"} (p-value),
#'   and \code{"ci"} (bootstrap confidence interval extracted from Mplus
#'   CINTERVAL output). Defaults to \code{c("est", "se", "p", "ci")}.
#' @param indices Character vector of additional fit indices to include beyond
#'   the default set (\eqn{\chi^2}, df, \eqn{\chi^2}/df, CFI, TLI, RMSEA,
#'   SRMR). Names must match columns in
#'   \code{model$summaries} (e.g.\ \code{"AIC"}, \code{"BIC"}).
#' @param digits_coeff Integer. Number of decimal places for coefficient
#'   estimates, standard errors, and confidence interval bounds.
#'   Defaults to \code{2}.
#' @param digits_fit Integer. Number of decimal places for fit indices.
#'   Defaults to \code{2}.
#' @param digits_pct Integer. Number of decimal places for effect proportions
#'   in the indirect-effects table. Defaults to \code{1}.
#' @param standardized Logical. If \code{TRUE} (default), report
#'   STDYX-standardised estimates for regression paths. MODEL CONSTRAINT (NEW)
#'   parameters are always reported as unstandardised because Mplus does not
#'   produce standardised estimates for them. If \code{FALSE}, unstandardised
#'   estimates are used throughout.
#' @param stars Logical. If \code{TRUE}, append significance stars to point
#'   estimates (\code{*} p < .05, \code{**} p < .01, \code{***} p < .001).
#'   Defaults to \code{TRUE}.
#' @param ci_level Numeric. Confidence level for extracted intervals.
#'   Currently only \code{0.95} is supported (Mplus outputs 2.5\% and 97.5\%
#'   percentiles). Reserved for future extension.
#' @param exclude Character vector of PCRE patterns. Predictors whose names
#'   match any pattern are excluded from the path-coefficient table.
#' @param title_coeff,title_fit,title_indirect Character strings used as table
#'   captions when rendering.
#' @param notes_coeff Character string appended as a general footnote to the
#'   path-coefficient table.
#' @param notes_indirect Character string appended as a general footnote to the
#'   indirect-effects table. When suppression (inconsistent mediation) is
#'   detected, an explanatory note is appended automatically.
#' @param include_indirect Logical. If \code{TRUE}, extract and display the
#'   indirect-effects decomposition table from the MODEL INDIRECT section. If
#'   \code{FALSE} (default), only the fit and coefficient tables are produced,
#'   which allows the function to work with any SEM or regression model that
#'   does not include a MODEL INDIRECT command.
#' @param render Logical. If \code{TRUE} (default), tables are printed as
#'   \pkg{kableExtra} HTML tables and the underlying list is returned
#'   invisibly. If \code{FALSE}, no rendering is performed and the list is
#'   returned visibly, which is convenient for passing results to other
#'   programs or AI-based interpretation.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{\code{fit}}{A one-row \code{data.frame} of model fit indices.
#'     Column names are human-readable labels (e.g.\ \eqn{\chi^2}, df, CFI).}
#'   \item{\code{coeff}}{A \code{data.frame} with predictors as rows and
#'     endogenous-variable \eqn{\times} statistic combinations as columns.
#'     When MODEL CONSTRAINT (NEW) parameters exist, they occupy the bottom
#'     rows (below row \code{.meta$n_reg}).}
#'   \item{\code{indirect}}{A \code{data.frame} with columns
#'     \code{Predictor}, \code{Outcome}, \code{Effect} (Total / Total Indirect
#'     / Specific Indirect / Direct), \code{Path} (arrow-delimited mediator
#'     chain), \code{est}, \code{se}, \code{p}, \code{ci}, and \code{pct}
#'     (effect proportion). \code{NULL} when \code{include_indirect = FALSE}
#'     or when no MODEL INDIRECT section exists in the output.}
#'   \item{\code{.meta}}{A list of rendering metadata: \code{endogenous}
#'     (character vector of endogenous variable names), \code{stat},
#'     \code{n_reg} (number of regression rows before NEW parameters),
#'     \code{has_new}, \code{standardized}, and \code{suppression} (logical
#'     flag for inconsistent mediation).}
#' }
#' When \code{render = TRUE} the list is returned \strong{invisibly}.
#'
#' @section Effect proportions:
#' When the direct and all specific indirect effects share the same sign
#' (consistent mediation), proportions are computed as
#' \eqn{\mathrm{effect} / \mathrm{total\;effect}}.
#' When the direct and any indirect effect have opposite signs (inconsistent
#' mediation / suppression), proportions are computed as
#' \eqn{|\mathrm{indirect} / \mathrm{direct}|} instead, and a footnote plus a
#' console message alert the user. The Total row never shows a proportion.
#'
#' @section Confidence intervals:
#' The function reads bootstrap (or Bayesian credibility) intervals directly
#' from the Mplus CINTERVAL output parsed by \pkg{MplusAutomation}. Symmetric
#' Wald-type intervals are \strong{not} computed. If the output file does not
#' contain a CINTERVAL section, the CI column will be blank.
#'
#' @seealso \code{\link[MplusAutomation]{readModels}} for reading Mplus output;
#'   \code{\link[kableExtra]{kable_classic}} for table styling.
#'
#' @examples
#' \dontrun{
#' ## -- Example Mplus syntax: serial mediation X -> M1 -> M2 -> Y --
#' ##
#' ## TITLE:  Serial mediation with bootstrap CI;
#' ## DATA:   FILE = mydata.dat;
#' ## VARIABLE:
#' ##   NAMES = x m1 m2 y c1 c2;
#' ##   USEVARIABLES = x m1 m2 y c1 c2;
#' ##
#' ## ANALYSIS:
#' ##   BOOTSTRAP = 5000;
#' ##
#' ## MODEL:
#' ##   m1 ON x c1 c2;
#' ##   m2 ON m1 x c1 c2;
#' ##   y  ON m2 m1 x c1 c2;
#' ##
#' ## MODEL INDIRECT:
#' ##   y IND m2 m1 x;   ! X -> M1 -> M2 -> Y  (full serial chain)
#' ##   y IND m1 x;       ! X -> M1 -> Y
#' ##   y IND m2 x;       ! X -> M2 -> Y
#' ##
#' ## MODEL CONSTRAINT:
#' ##   NEW(diff);
#' ##   diff = ind_m1 - ind_m2;  ! difference between two indirect paths
#' ##
#' ## OUTPUT:
#' ##   CINTERVAL(BCBOOTSTRAP);
#' ##   STANDARDIZED;
#'
#' ## -- Read the output in R --
#' library(MplusAutomation)
#' mod <- readModels("serial_med.out")
#'
#' ## -- Default: fit + standardised coefficients only --
#' coeff_table(mod)
#'
#' ## -- Full mediation report with indirect effects --
#' coeff_table(mod, include_indirect = TRUE)
#'
#' ## -- Unstandardised, extra fit indices --
#' coeff_table(mod, standardized = FALSE, indices = c("AIC", "BIC"),
#'                 include_indirect = TRUE)
#'
#' ## -- Raw data.frames for AI interpretation --
#' res <- coeff_table(mod, include_indirect = TRUE, render = FALSE)
#' res$fit
#' res$coeff
#' res$indirect
#' }
#'
#' @importFrom knitr kable
#' @importFrom kableExtra kable_classic add_header_above row_spec footnote
#' @export
coeff_table <- function(
    model,
    stat           = c("est", "se", "p", "ci"),
    indices        = character(0),
    digits_coeff   = 2,
    digits_fit     = 2,
    digits_pct     = 1,
    standardized   = TRUE,
    stars          = TRUE,
    ci_level       = 0.95,
    exclude        = character(0),
    title_coeff    = "Table X. Path coefficients.",
    title_fit      = "Table X. Model fit indices.",
    title_indirect = "Table X. Direct, indirect, and total effects.",
    notes_coeff    = NULL,
    notes_indirect = NULL,
    include_indirect = FALSE,
    render         = TRUE
) {

  ## =========================================================================
  ## 0. Validate inputs
  ## =========================================================================

  if (!inherits(model, c("mplus.model", "list")) || is.null(model$summaries))
    stop("model must be a single mplus.model object returned by readModels().")

  allowed_stat <- c("est", "se", "z", "p", "ci")
  stat <- tolower(stat)
  if (!all(stat %in% allowed_stat))
    stop("`stat` must contain only: ", paste(allowed_stat, collapse = ", "))

  ## =========================================================================
  ## Helper: formatting utilities
  ## =========================================================================

  .fmt_num <- function(x, digits = 3)
    ifelse(is.na(x), "", sprintf(paste0("%.", digits, "f"), x))

  .fmt_p <- function(p)
    ifelse(is.na(p), "",
           ifelse(p < .001, "< .001", sub("^0", "", sprintf("%.3f", p))))

  .add_stars <- function(est_str, p_num) {
    s <- ifelse(is.na(p_num), "",
                ifelse(p_num < .001, "***",
                       ifelse(p_num < .01, "**",
                              ifelse(p_num < .05, "*", ""))))
    paste0(est_str, s)
  }

  .fmt_ci <- function(lo, hi, digits = 2) {
    ifelse(is.na(lo) | is.na(hi), "",
           paste0("[", sprintf(paste0("%.", digits, "f"), lo),
                  ", ", sprintf(paste0("%.", digits, "f"), hi), "]"))
  }

  ## =========================================================================
  ## Helper: pick 95% CI columns
  ## =========================================================================

  .pick_ci_cols <- function(ci_df) {
    nms <- names(ci_df)
    lo_col <- grep("low2\\.5", nms, value = TRUE)[1]
    hi_col <- grep("up2\\.5",  nms, value = TRUE)[1]
    if (is.na(lo_col) || is.na(hi_col)) {
      num_cols <- nms[!nms %in% c("paramHeader", "param", "outcome",
                                  "LatentClass", "BetweenWithin",
                                  "Group", "pred", "intervening",
                                  "summary", "ReferenceClass")]
      lo_col <- num_cols[1]
      hi_col <- num_cols[length(num_cols)]
    }
    list(lo = lo_col, hi = hi_col)
  }

  ## =========================================================================
  ## Helper: parse paramHeader -> lhs / op / rhs
  ## =========================================================================

  .create_formula_parts <- function(header, param) {
    if (grepl("\\.BY$", header)) {
      c(lhs = sub("\\.BY$", "", header), op = "=~", rhs = param)
    } else if (grepl("\\.ON$", header)) {
      c(lhs = sub("\\.ON$", "", header), op = "~", rhs = param)
    } else if (grepl("\\.WITH$", header)) {
      c(lhs = sub("\\.WITH$", "", header), op = "~~", rhs = param)
    } else if (header %in% c("Intercepts", "Means")) {
      c(lhs = param, op = "~ 1", rhs = "")
    } else if (header %in% c("Variances", "Residual.Variances")) {
      c(lhs = param, op = "~~", rhs = param)
    } else if (header == "New.Additional.Parameters") {
      c(lhs = "new", op = ":=", rhs = param)
    } else {
      c(lhs = header, op = "", rhs = param)
    }
  }

  ## =========================================================================
  ## Helper: attach lhs/op/rhs to a parameter data.frame
  ## =========================================================================

  .attach_parts <- function(df) {
    parts <- t(mapply(.create_formula_parts, df$paramHeader, df$param))
    cbind(as.data.frame(parts, stringsAsFactors = FALSE), df)
  }

  ## =========================================================================
  ## PART 1:  Model Fit Table
  ## =========================================================================

  .build_fit_table <- function(model, indices, digits_fit) {

    core <- c("ChiSqM_Value", "ChiSqM_DF", "chisq_df",
              "CFI", "TLI", "RMSEA_Estimate", "SRMR")
    want <- unique(c(core, indices))

    rename_map <- c(ChiSqM_Value   = "\u03C7\u00B2",
                    ChiSqM_DF      = "df",
                    chisq_df       = "\u03C7\u00B2/df",
                    RMSEA_Estimate = "RMSEA")

    s <- model$summaries

    pull_one <- function(idx) {
      if (idx == "chisq_df" &&
          all(c("ChiSqM_Value", "ChiSqM_DF") %in% names(s)) &&
          !is.na(s$ChiSqM_DF) && s$ChiSqM_DF != 0)
        return(s$ChiSqM_Value / s$ChiSqM_DF)
      if (idx %in% names(s)) return(s[[idx]])
      NA_real_
    }

    raw_vals <- vapply(want, pull_one, numeric(1))
    labels   <- ifelse(want %in% names(rename_map), rename_map[want], want)

    fmt <- vapply(seq_along(labels), function(i) {
      v <- raw_vals[i]
      if (is.na(v)) return("")
      pat <- paste0("%.", digits_fit, "f")
      if (labels[i] == "df") return(as.character(as.integer(v)))
      if (labels[i] %in% c("CFI", "TLI", "RMSEA", "SRMR") && v < 1)
        return(sub("^0", "", sprintf(pat, v)))
      sprintf(pat, v)
    }, character(1))

    df_fit <- data.frame(t(fmt), stringsAsFactors = FALSE)
    names(df_fit) <- labels
    df_fit
  }

  fit_table <- .build_fit_table(model, indices, digits_fit)

  ## =========================================================================
  ## PART 2:  Path Coefficient Table (regression + NEW params merged)
  ## =========================================================================

  .build_coeff_table <- function(model, stat, standardized, digits_coeff,
                                 stars, exclude) {

    ## 2a. Get regression parameter table --------------------------------------
    tbl <- if (standardized) {
      model$parameters$stdyx.standardized
    } else {
      model$parameters$unstandardized
    }
    if (is.null(tbl))
      stop("Requested parameter table not found in model output.")

    tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)
    names(tbl)[names(tbl) == "pval"] <- "p"
    tbl <- .attach_parts(tbl)

    ## 2b. Get CI for regression (respects standardized toggle) ----------------
    ci_reg <- NULL
    if ("ci" %in% stat) {
      ci_key <- if (standardized) "ci.stdyx.standardized" else "ci.unstandardized"
      ci_raw <- model$parameters[[ci_key]]
      if (!is.null(ci_raw)) {
        ci_raw <- as.data.frame(ci_raw, stringsAsFactors = FALSE)
        ci_raw <- .attach_parts(ci_raw)
        cc <- .pick_ci_cols(ci_raw)
        ci_raw$ci_lo <- as.numeric(ci_raw[[cc$lo]])
        ci_raw$ci_hi <- as.numeric(ci_raw[[cc$hi]])
        ci_reg <- ci_raw[, c("lhs", "op", "rhs", "ci_lo", "ci_hi"), drop = FALSE]
      }
    }

    ## 2c. Subset regression rows, apply exclude --------------------------------
    reg_tbl <- tbl[tbl$op == "~", , drop = FALSE]
    if (length(exclude) > 0) {
      pat <- paste(exclude, collapse = "|")
      reg_tbl <- reg_tbl[!grepl(pat, reg_tbl$rhs, ignore.case = TRUE, perl = TRUE), ]
    }

    ## merge CI
    if (!is.null(ci_reg) && "ci" %in% stat) {
      reg_tbl <- merge(reg_tbl, ci_reg, by = c("lhs", "op", "rhs"), all.x = TRUE)
    }

    ## 2d. Identify endogenous variables and predictors -------------------------
    endogenous     <- unique(reg_tbl$lhs)
    all_predictors <- unique(reg_tbl$rhs)

    stat_display <- stat

    ## 2e. Build wide-format regression block -----------------------------------
    wide_list <- list()
    for (dv in endogenous) {
      sub <- reg_tbl[reg_tbl$lhs == dv, , drop = FALSE]

      p_numeric <- sub$p
      est_fmt   <- .fmt_num(sub$est, digits_coeff)
      if (stars) est_fmt <- .add_stars(est_fmt, p_numeric)

      block <- data.frame(predictor = sub$rhs, stringsAsFactors = FALSE)
      for (s in stat_display) {
        cn <- paste0(dv, "_", s)
        if (s == "est")     block[[cn]] <- est_fmt
        else if (s == "se") block[[cn]] <- .fmt_num(sub$se, digits_coeff)
        else if (s == "z")  block[[cn]] <- .fmt_num(sub$est_se, digits_coeff)
        else if (s == "p")  block[[cn]] <- .fmt_p(p_numeric)
        else if (s == "ci") {
          block[[cn]] <- if ("ci_lo" %in% names(sub))
            .fmt_ci(sub$ci_lo, sub$ci_hi, digits_coeff) else ""
        }
      }
      wide_list[[dv]] <- block
    }

    merged <- Reduce(function(a, b) merge(a, b, by = "predictor", all = TRUE), wide_list)
    merged[is.na(merged)] <- ""
    merged$predictor <- factor(merged$predictor, levels = all_predictors)
    merged <- merged[order(merged$predictor), , drop = FALSE]
    merged$predictor <- as.character(merged$predictor)

    ## 2f. NEW parameters (always unstandardized) ------------------------------
    n_reg   <- nrow(merged)
    has_new <- FALSE

    new_src <- as.data.frame(model$parameters$unstandardized,
                             stringsAsFactors = FALSE)
    names(new_src)[names(new_src) == "pval"] <- "p"
    new_src <- .attach_parts(new_src)
    new_src <- new_src[new_src$op == ":=", , drop = FALSE]

    if (nrow(new_src) > 0) {
      has_new <- TRUE

      ## CI for NEW: always ci.unstandardized
      ci_new_lo <- rep(NA_real_, nrow(new_src))
      ci_new_hi <- rep(NA_real_, nrow(new_src))
      if ("ci" %in% stat) {
        ci_un <- model$parameters[["ci.unstandardized"]]
        if (!is.null(ci_un)) {
          ci_un <- as.data.frame(ci_un, stringsAsFactors = FALSE)
          ci_un <- .attach_parts(ci_un)
          ci_un <- ci_un[ci_un$op == ":=", , drop = FALSE]
          if (nrow(ci_un) > 0) {
            cc_n <- .pick_ci_cols(ci_un)
            ci_un$ci_lo <- as.numeric(ci_un[[cc_n$lo]])
            ci_un$ci_hi <- as.numeric(ci_un[[cc_n$hi]])
            idx <- match(new_src$rhs, ci_un$rhs)
            ci_new_lo <- ci_un$ci_lo[idx]
            ci_new_hi <- ci_un$ci_hi[idx]
          }
        }
      }

      p_num_new <- new_src$p
      est_new   <- .fmt_num(new_src$est, digits_coeff)
      if (stars) est_new <- .add_stars(est_new, p_num_new)

      new_block <- data.frame(predictor = new_src$rhs, stringsAsFactors = FALSE)

      first_dv <- endogenous[1]
      for (s in stat_display) {
        cn <- paste0(first_dv, "_", s)
        if (s == "est")     new_block[[cn]] <- est_new
        else if (s == "se") new_block[[cn]] <- .fmt_num(new_src$se, digits_coeff)
        else if (s == "z")  new_block[[cn]] <- .fmt_num(new_src$est_se, digits_coeff)
        else if (s == "p")  new_block[[cn]] <- .fmt_p(p_num_new)
        else if (s == "ci") new_block[[cn]] <- .fmt_ci(ci_new_lo, ci_new_hi, digits_coeff)
      }

      for (cn in setdiff(names(merged), names(new_block)))
        new_block[[cn]] <- ""
      new_block <- new_block[, names(merged), drop = FALSE]

      merged <- rbind(merged, new_block)
      row.names(merged) <- NULL
    }

    list(coeff      = merged,
         endogenous = endogenous,
         stat       = stat_display,
         n_reg      = n_reg,
         has_new    = has_new)
  }

  coeff_result <- .build_coeff_table(model, stat, standardized, digits_coeff,
                                     stars, exclude)

  ## =========================================================================
  ## PART 3:  Indirect Effects Table (conditional on include_indirect)
  ## =========================================================================

  .build_indirect_table <- function(model, standardized, digits_coeff,
                                    stars, digits_pct) {

    ind_key <- if (standardized) "stdyx.standardized" else "unstandardized"
    ind_section <- model$indirect[[ind_key]]

    if (is.null(ind_section)) {
      warning("No indirect effects found for type: ", ind_key,
              ". Check that MODEL INDIRECT is specified and the output ",
              "contains the corresponding section.")
      return(NULL)
    }

    overall  <- ind_section$overall
    specific <- ind_section$specific

    ## CI
    ci_key <- if (standardized) "ci.stdyx.standardized" else "ci.unstandardized"
    ci_section <- model$indirect[[ci_key]]

    ci_overall  <- NULL
    ci_specific <- NULL
    if (!is.null(ci_section)) {
      if (!is.null(ci_section$overall)) {
        ci_overall <- ci_section$overall
        cc <- .pick_ci_cols(ci_overall)
        ci_overall$ci_lo <- as.numeric(ci_overall[[cc$lo]])
        ci_overall$ci_hi <- as.numeric(ci_overall[[cc$hi]])
      }
      if (!is.null(ci_section$specific)) {
        ci_specific <- ci_section$specific
        cc_sp <- .pick_ci_cols(ci_specific)
        ci_specific$ci_lo <- as.numeric(ci_specific[[cc_sp$lo]])
        ci_specific$ci_hi <- as.numeric(ci_specific[[cc_sp$hi]])
      }
    }

    .norm_p <- function(df) {
      if (!is.null(df) && "pval" %in% names(df) && !"p" %in% names(df))
        names(df)[names(df) == "pval"] <- "p"
      df
    }
    overall  <- .norm_p(overall)
    specific <- .norm_p(specific)

    rows_list <- list()
    suppression_flag <- FALSE
    pairs <- unique(overall[, c("pred", "outcome")])

    for (r in seq_len(nrow(pairs))) {
      pr <- pairs$pred[r]
      oc <- pairs$outcome[r]

      ov_sub <- overall[overall$pred == pr & overall$outcome == oc, , drop = FALSE]

      total_row <- ov_sub[ov_sub$summary == "Total", , drop = FALSE]
      tind_row  <- ov_sub[grepl("indirect|Indirect", ov_sub$summary), , drop = FALSE]
      dir_row   <- ov_sub[ov_sub$summary == "Direct", , drop = FALSE]

      total_est <- if (nrow(total_row) > 0) as.numeric(total_row$est[1]) else NA_real_
      dir_est   <- if (nrow(dir_row) > 0) as.numeric(dir_row$est[1]) else NA_real_

      spec_sub <- if (!is.null(specific)) {
        specific[specific$pred == pr & specific$outcome == oc, , drop = FALSE]
      } else { data.frame() }

      is_suppression <- FALSE
      if (nrow(spec_sub) > 0 && !is.na(dir_est)) {
        spec_ests <- as.numeric(spec_sub$est)
        if (any(!is.na(spec_ests) & spec_ests != 0 & dir_est != 0 &
                sign(spec_ests) != sign(dir_est))) {
          is_suppression <- TRUE
          suppression_flag <- TRUE
        }
      }

      .make_row <- function(effect_type, path_label, est_val, se_val,
                            p_val, ci_lo, ci_hi, show_pct = TRUE) {
        est_n <- as.numeric(est_val)
        se_n  <- as.numeric(se_val)
        p_n   <- as.numeric(p_val)

        pct_str <- ""
        if (show_pct && effect_type != "Total") {
          if (!is.na(est_n)) {
            if (is_suppression) {
              if (!is.na(dir_est) && dir_est != 0)
                pct_str <- paste0(sprintf(paste0("%.", digits_pct, "f"),
                                          abs(est_n / dir_est) * 100), "%")
            } else {
              if (!is.na(total_est) && abs(total_est) > 1e-10)
                pct_str <- paste0(sprintf(paste0("%.", digits_pct, "f"),
                                          (est_n / total_est) * 100), "%")
            }
          }
        }

        est_fmt <- .fmt_num(est_n, digits_coeff)
        if (stars) est_fmt <- .add_stars(est_fmt, p_n)

        data.frame(
          Predictor = pr, Outcome = oc,
          Effect = effect_type, Path = path_label,
          est = est_fmt, se = .fmt_num(se_n, digits_coeff),
          p = .fmt_p(p_n), ci = .fmt_ci(ci_lo, ci_hi, digits_coeff),
          pct = pct_str,
          stringsAsFactors = FALSE
        )
      }

      .ci_overall_lookup <- function(summary_pattern) {
        lo <- NA_real_; hi <- NA_real_
        if (!is.null(ci_overall)) {
          m <- ci_overall[ci_overall$pred == pr & ci_overall$outcome == oc &
                            grepl(summary_pattern, ci_overall$summary,
                                  ignore.case = TRUE), , drop = FALSE]
          if (nrow(m) > 0) { lo <- m$ci_lo[1]; hi <- m$ci_hi[1] }
        }
        list(lo = lo, hi = hi)
      }

      ## Total
      if (nrow(total_row) > 0) {
        ci_t <- .ci_overall_lookup("^Total$")
        rows_list[[length(rows_list) + 1]] <-
          .make_row("Total", "", total_row$est[1], total_row$se[1],
                    total_row$p[1], ci_t$lo, ci_t$hi, show_pct = FALSE)
      }
      ## Total Indirect
      if (nrow(tind_row) > 0) {
        ci_ti <- .ci_overall_lookup("indirect|Indirect")
        rows_list[[length(rows_list) + 1]] <-
          .make_row("Total Indirect", "", tind_row$est[1], tind_row$se[1],
                    tind_row$p[1], ci_ti$lo, ci_ti$hi)
      }
      ## Specific Indirect(s)
      if (nrow(spec_sub) > 0) {
        for (sp in seq_len(nrow(spec_sub))) {
          path_label <- paste(pr, "->",
                              gsub("\\.", " -> ", spec_sub$intervening[sp]),
                              "->", oc)
          ci_s_lo <- NA_real_; ci_s_hi <- NA_real_
          if (!is.null(ci_specific)) {
            ci_m <- ci_specific[ci_specific$pred == pr &
                                  ci_specific$outcome == oc &
                                  ci_specific$intervening == spec_sub$intervening[sp],
                                , drop = FALSE]
            if (nrow(ci_m) > 0) { ci_s_lo <- ci_m$ci_lo[1]; ci_s_hi <- ci_m$ci_hi[1] }
          }
          rows_list[[length(rows_list) + 1]] <-
            .make_row("Specific Indirect", path_label,
                      spec_sub$est[sp], spec_sub$se[sp],
                      spec_sub$p[sp], ci_s_lo, ci_s_hi)
        }
      }
      ## Direct
      if (nrow(dir_row) > 0) {
        ci_d <- .ci_overall_lookup("^Direct$")
        rows_list[[length(rows_list) + 1]] <-
          .make_row("Direct", paste(pr, "->", oc),
                    dir_row$est[1], dir_row$se[1],
                    dir_row$p[1], ci_d$lo, ci_d$hi)
      }
    }

    indirect_df <- do.call(rbind, rows_list)
    row.names(indirect_df) <- NULL
    attr(indirect_df, "suppression") <- suppression_flag
    indirect_df
  }

  ## Only build indirect table when requested
  indirect_table <- NULL
  if (include_indirect) {
    indirect_table <- .build_indirect_table(model, standardized, digits_coeff,
                                            stars, digits_pct)
  }

  ## =========================================================================
  ## PART 4:  Assemble output
  ## =========================================================================

  result <- list(
    fit      = fit_table,
    coeff    = coeff_result$coeff,
    indirect = indirect_table,
    .meta    = list(
      endogenous   = coeff_result$endogenous,
      stat         = coeff_result$stat,
      n_reg        = coeff_result$n_reg,
      has_new      = coeff_result$has_new,
      standardized = standardized,
      suppression  = if (!is.null(indirect_table)) attr(indirect_table, "suppression") else FALSE
    )
  )

  ## =========================================================================
  ## PART 5:  Render
  ## =========================================================================

  if (!render) return(result)

  ## --- 5a. Fit table --------------------------------------------------------
  fit_kbl <- knitr::kable(fit_table, format = "html",
                          caption = title_fit, align = "c") %>%
    kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria",
                              position = "left")

  ## --- 5b. Coefficient table ------------------------------------------------
  coeff_df   <- coeff_result$coeff
  endogenous <- coeff_result$endogenous
  stat_disp  <- coeff_result$stat
  n_reg      <- coeff_result$n_reg
  has_new    <- coeff_result$has_new

  mapping <- if (standardized) {
    c(est = "\u03B2", se = "SE", z = "z", p = "p", ci = "95% CI")
  } else {
    c(est = "b", se = "SE", z = "z", p = "p", ci = "95% CI")
  }
  stat_labels <- mapping[stat_disp]

  coeff_colnames <- c("Predictor", rep(stat_labels, times = length(endogenous)))

  coeff_header <- c(" " = 1)
  for (dv in endogenous) coeff_header[dv] <- length(stat_disp)

  ## auto-generate notes_coeff if NULL
  if (is.null(notes_coeff)) {
    notes_coeff <- if (standardized) {
      "Standardized coefficients are reported."
    } else {
      "Unstandardized coefficients are reported."
    }
  }
  coeff_notes_full <- notes_coeff
  if (has_new)
    coeff_notes_full <- paste0(coeff_notes_full,
                               " Parameters below the line are from MODEL CONSTRAINT (unstandardized).")

  coeff_kbl <- knitr::kable(coeff_df, format = "html",
                            col.names = coeff_colnames,
                            caption = title_coeff,
                            align = c("l", rep("c", ncol(coeff_df) - 1))) %>%
    kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria",
                              position = "left") %>%
    kableExtra::add_header_above(coeff_header)

  if (has_new) {
    coeff_kbl <- coeff_kbl %>%
      kableExtra::row_spec(n_reg, extra_css = "border-bottom: 1px solid;")
  }

  coeff_kbl <- coeff_kbl %>%
    kableExtra::footnote(general = coeff_notes_full, footnote_as_chunk = TRUE)

  ## --- 5c. Indirect effects table (only when include_indirect = TRUE) -------
  ind_kbl <- NULL
  if (!is.null(indirect_table)) {
    ind_colnames <- c("Predictor", "Outcome", "Effect", "Path",
                      mapping["est"], "SE", "p", "95% CI", "Proportion")

    ind_kbl <- knitr::kable(indirect_table, format = "html",
                            col.names = ind_colnames,
                            caption = title_indirect,
                            align = c("l","l","l","l","c","c","c","c","c")) %>%
      kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria",
                                position = "left")

    ind_notes <- if (is.null(notes_indirect)) character(0) else notes_indirect
    if (isTRUE(attr(indirect_table, "suppression"))) {
      supp_note <- paste0(
        "Suppression (inconsistent mediation) detected: ",
        "the direct and indirect effects have opposite signs. ",
        "Proportions are computed as |indirect / direct| instead of effect / total.")
      ind_notes <- c(ind_notes, supp_note)
      message("Note: Suppression effect detected. ",
              "Effect proportions use |indirect / direct| instead of effect / total.")
    }

    if (length(ind_notes) > 0) {
      ind_kbl <- ind_kbl %>%
        kableExtra::footnote(general = paste(ind_notes, collapse = " "),
                             footnote_as_chunk = TRUE)
    }
  }

  ## --- 5d. Print -------------------------------------------------------------
  cat("\n")
  print(fit_kbl)
  cat("\n\n")
  print(coeff_kbl)
  if (!is.null(ind_kbl)) { cat("\n\n"); print(ind_kbl) }

  invisible(result)
}
