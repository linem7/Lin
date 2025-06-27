#' Fit Table of Model Fit Indices and Difference Metrics
#'
#' This function creates a table of fit indices for one or more models and computes
#' difference metrics (e.g., ΔCFI, ΔRMSEA, Δχ²) between comparative models. If the
#' primary chi-square values or degrees of freedom are unavailable, the function will
#' automatically switch to a log-likelihood (LL) based approach. When available, the
#' LLCorrectionFactor is used; otherwise, the LL difference is computed using LL alone.
#'
#' The final output always displays the default fit indices:
#' \code{"Models", "χ²/df", "χ²/df", "CFI", "TLI", "RMSEA", "SRMR"}
#' and the default difference metrics:
#' \code{"Δχ²", "p_Δχ²", "ΔCFI", "ΔRMSEA"}. If the user specifies additional indices
#' (e.g., "AIC", "BIC") via the \code{indices} argument, then those indices will be included
#' in the fit chunk, and their corresponding difference metrics (ΔAIC and/or ΔBIC) will be added
#' to the diff chunk.
#'
#' @param ... One or more model objects (or lists of models) that contain a \code{summaries} element.
#' @param indices Character vector of additional indices to report. These will be added after the default
#'   fit indices. If "AIC" or "BIC" are included, the corresponding difference metrics (ΔAIC, ΔBIC) are
#'   also reported.
#' @param ref Character string indicating the reference model to compare to. Options are "last" (default)
#'   or "first".
#' @param diffTest Logical. If \code{TRUE} (default), difference metrics will be computed.
#' @param digits Integer. Number of digits to round numeric outputs.
#' @param model_names Optional character vector of custom model names. If supplied, its length must equal the number of models passed via `...`; these names will override the default names (extracted from each model’s `Filename` field). Defaults to `NULL`.
#'
#' @return A data frame containing the fit indices and difference metrics.
#'
#' @examples
#' \dontrun{
#'   # Assuming mod1 and mod2 are valid model objects with a 'summaries' element:
#'   fit_table(mod1, mod2, indices = c("AIC", "BIC"), ref = "first", digits = 2)
#' }
#'
#' @export
fit_table <- function(...,
                      indices      = character(0),     # extra indices to show
                      ref          = c("last", "first"),
                      diffTest     = TRUE,
                      digits       = 2,
                      model_names  = NULL) {

  ref <- match.arg(ref)
  args <- list(...)

  ## ---------------- helpers ---------------------------------------------
  .get_num <- function(df, nm) {
    if (nm %in% names(df)) as.numeric(df[[nm]][1]) else NA_real_
  }
  .strip0 <- function(x) sub("^(-?)0\\.", "\\1.", x)
  .fnum   <- function(x, d) sprintf(paste0("%.", d, "f"), x)

  .fmt_fit <- function(x, d) {
    vapply(x, function(xx) {
      if (is.na(xx)) return(NA_character_)
      out <- .fnum(xx, d)
      if (xx < 1) out <- .strip0(out)
      out
    }, character(1))
  }

  .fmt_diff <- function(x, d) {
    vapply(x, function(xx) {
      if (is.na(xx)) return(NA_character_)
      if (abs(xx) < 10^(-d))
        return(paste0(".", paste(rep("0", d), collapse = "")))
      .strip0(.fnum(xx, d))
    }, character(1))
  }

  .fmt_p <- function(p, d = 3) {
    vapply(p, function(pp) {
      if (is.na(pp)) return(NA_character_)
      if (pp < .001) return("< .001")
      .strip0(.fnum(pp, d))
    }, character(1))
  }

  ## --------------- gather model objects --------------------------------
  modelList <- list()
  for (obj in args) {
    if (is.list(obj) && !("summaries" %in% names(obj))) {
      modelList <- c(modelList, obj)
    } else if (is.list(obj) && ("summaries" %in% names(obj))) {
      modelList[[length(modelList) + 1]] <- obj
    } else {
      stop("Input is not a valid model or list of models.")
    }
  }
  n <- length(modelList)
  if (n == 0) stop("No models provided.")

  ## --------------- decide which defaults to show -----------------------
  sm0       <- as.data.frame(modelList[[1]]$summaries)
  analysis0 <- modelList[[1]]$input$analysis %||% ""
  has_classic <- all(c("ChiSqM_Value","ChiSqM_DF",
                       "CFI","TLI","RMSEA_Estimate","SRMR") %in% names(sm0)) &&
    !grepl("random", analysis0, ignore.case = TRUE)

  if (has_classic) {
    base_idx <- c("Models", "χ²(df)", "χ²/df", "CFI", "TLI", "RMSEA", "SRMR")
  } else {
    base_idx <- c("Models", "LL", "Parameters", "AIC", "BIC", "aBIC")
  }
  add_idx <- indices

  # build the diff‐columns list
  if (diffTest && n > 1) {
    diff_idx <- c("Δχ²", "p_Δχ²", "ΔCFI", "ΔRMSEA")
    if ("AIC" %in% add_idx) diff_idx <- union(diff_idx, "ΔAIC")
    if ("BIC" %in% add_idx) diff_idx <- union(diff_idx, "ΔBIC")
  } else {
    diff_idx <- character(0)
  }

  all_cols <- unique(c(base_idx, add_idx, diff_idx))
  holder   <- setNames(
    data.frame(matrix(NA_character_, n, length(all_cols)),
               stringsAsFactors = FALSE),
    all_cols
  )

  ## --------------- fill each model’s row ------------------------------
  for (i in seq_len(n)) {
    sm <- as.data.frame(modelList[[i]]$summaries)
    holder$Models[i] <- if (!is.null(model_names)) model_names[i] else sm$Filename[1]

    if (has_classic) {
      cs_val <- .get_num(sm, "ChiSqM_Value")
      df_val <- as.integer(.get_num(sm, "ChiSqM_DF"))
      if (!is.na(cs_val) && !is.na(df_val) && df_val != 0) {
        holder$`χ²(df)`[i] <- paste0(.fnum(cs_val, digits), " (", df_val, ")")
        holder$`χ²/df` [i] <- .fnum(cs_val / df_val, digits)
      }
      holder$CFI  [i] <- .fmt_fit(.get_num(sm, "CFI"),            digits)
      holder$TLI  [i] <- .fmt_fit(.get_num(sm, "TLI"),            digits)
      holder$RMSEA[i] <- .fmt_fit(.get_num(sm, "RMSEA_Estimate"), digits)
      holder$SRMR [i] <- .fmt_fit(.get_num(sm, "SRMR"),           digits)
    } else {
      holder$LL         [i] <- .fnum(.get_num(sm, "LL"),   digits)
      holder$Parameters [i] <- as.character(as.integer(.get_num(sm, "Parameters")))
      holder$AIC        [i] <- .fnum(.get_num(sm, "AIC"),  digits)
      holder$BIC        [i] <- .fnum(.get_num(sm, "BIC"),  digits)
      holder$aBIC       [i] <- .fnum(.get_num(sm, "aBIC"), digits)
    }

    # extra indices, with integer‐only for Parameters/df
    for (idx in add_idx) {
      if (!(idx %in% names(sm))) next
      val <- .get_num(sm, idx)
      if (idx %in% c("Parameters", "ChiSqM_DF")) {
        holder[[idx]][i] <- as.character(as.integer(val))
      } else {
        holder[[idx]][i] <- .fnum(val, digits)
      }
    }
  }

  ## --------------- difference testing ----------------------------------
  if (diffTest && n > 1) {
    # compute deltas
    delta_CFI   <- delta_RMSEA <- delta_AIC <- delta_BIC <-
      delta_chisq <- delta_df <- p_delta_chisq <- rep(NA_real_, n)

    for (i in 2:n) {
      ref_i <- if (ref == "first") 1 else i - 1
      s1 <- as.data.frame(modelList[[i]]$summaries)
      s0 <- as.data.frame(modelList[[ref_i]]$summaries)

      delta_CFI   [i] <- .get_num(s1, "CFI")            - .get_num(s0, "CFI")
      delta_RMSEA [i] <- .get_num(s1, "RMSEA_Estimate") - .get_num(s0, "RMSEA_Estimate")
      delta_AIC   [i] <- .get_num(s1, "AIC")            - .get_num(s0, "AIC")
      delta_BIC   [i] <- .get_num(s1, "BIC")            - .get_num(s0, "BIC")

      if (!is.na(.get_num(s1, "ChiSqM_Value")) &&
          !is.na(.get_num(s0, "ChiSqM_Value"))) {
        delta_chisq[i] <- .get_num(s1, "ChiSqM_Value") - .get_num(s0, "ChiSqM_Value")
        delta_df   [i] <- as.integer(
          .get_num(s1, "ChiSqM_DF") - .get_num(s0, "ChiSqM_DF")
        )
        if (!is.na(delta_df[i]) && delta_df[i] > 0) {
          p_delta_chisq[i] <- 1 - pchisq(delta_chisq[i], delta_df[i])
        }
      }
    }

    # format into holder
    if ("ΔCFI"   %in% all_cols) holder$`ΔCFI`   <- .fmt_diff(delta_CFI,   digits)
    if ("ΔRMSEA"%in% all_cols) holder$`ΔRMSEA` <- .fmt_diff(delta_RMSEA, digits)
    if ("ΔAIC"  %in% all_cols) holder$`ΔAIC`   <- .fmt_diff(delta_AIC,   digits)
    if ("ΔBIC"  %in% all_cols) holder$`ΔBIC`   <- .fmt_diff(delta_BIC,   digits)
    if ("Δχ²"   %in% all_cols) holder$`Δχ²`    <- ifelse(
      !is.na(delta_chisq) & !is.na(delta_df),
      paste0(.fnum(delta_chisq, digits), " (", delta_df, ")"),
      NA_character_
    )
    if ("p_Δχ²" %in% all_cols) holder$`p_Δχ²`  <- .fmt_p(p_delta_chisq)

    # first‐row placeholders
    to_fill <- intersect(c("Δχ²","p_Δχ²","ΔCFI","ΔRMSEA","ΔAIC","ΔBIC"), all_cols)
    holder[1, to_fill] <- "—"
  }

  ## --------------- return requested columns ----------------------------
  holder[, all_cols, drop = FALSE]
}


