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
                      indices      = character(0),   # extra indices to show
                      ref          = c("last", "first"),
                      diffTest     = TRUE,
                      digits       = 2,
                      model_names  = NULL) {

  ref <- match.arg(ref)
  args <- list(...)

  ## ---------------- helpers ---------------------------------------------
  # remove leading zero for both positive (0.83 → .83) and negative (-0.83 → -.83)
  .strip0 <- function(x) sub("^(-?)0\\.", "\\1.", x)
  .fnum   <- function(x, d) sprintf(paste0("%.", d, "f"), x)

  ## APA-style for CFI / TLI / RMSEA / SRMR (vectorised)
  .fmt_fit <- function(x, d) {
    vapply(x, function(xx) {
      if (is.na(xx)) return(NA_character_)
      out <- .fnum(xx, d)
      if (xx < 1) out <- .strip0(out)
      out
    }, character(1))
  }

  ## APA-style for Δ columns (vectorised, never "< .001")
  .fmt_diff <- function(x, d) {
    vapply(x, function(xx) {
      if (is.na(xx)) return(NA_character_)
      if (abs(xx) < 10^(-d)) return(paste0(".", paste(rep("0", d), collapse = "")))
      .strip0(.fnum(xx, d))
    }, character(1))
  }

  ## p-value formatting (always 3‑digit when printed)
  .fmt_p <- function(p, d = 3) {
    vapply(p, function(pp) {
      if (is.na(pp)) return(NA_character_)
      if (pp < .001) return("< .001")
      .strip0(.fnum(pp, d))
    }, character(1))
  }

  ## --------------- gather model objects ----------------------------------
  modelList <- list()
  for (obj in args) {
    if (is.list(obj) && !("summaries" %in% names(obj))) {
      modelList <- c(modelList, obj)
    } else if (is.list(obj) && ("summaries" %in% names(obj))) {
      modelList[[length(modelList) + 1]] <- obj
    } else stop("Input is not a valid model or list of models.")
  }
  n <- length(modelList); if (n == 0) stop("No models provided.")

  ## --------------- holder ------------------------------------------------
  def_cols <- c("Models", "χ²(df)", "χ²/df",
                "CFI", "TLI", "RMSEA", "SRMR",
                "LL", "Parameters", "AIC", "BIC")
  holder   <- as.data.frame(matrix(NA, n, length(def_cols)),
                            stringsAsFactors = FALSE)
  names(holder) <- def_cols

  ## --------------- fill fit indices --------------------------------------
  for (i in seq_len(n)) {
    sm <- as.data.frame(modelList[[i]]$summaries)
    holder$Models[i] <- if (!is.null(model_names)) model_names[i] else sm$Filename[1]

    ## χ²(df) + χ²/df
    cs <- as.numeric(sm$ChiSqM_Value[1]); df <- as.numeric(sm$ChiSqM_DF[1])
    if (!is.na(cs) && !is.na(df))
      holder$`χ²(df)`[i] <- paste0(.fnum(cs, digits), " (", df, ")")
    if (!is.na(cs) && !is.na(df) && df != 0)
      holder$`χ²/df`[i] <- .fnum(cs / df, digits)

    ## 4 main fit indices
    holder$CFI[i]   <- .fmt_fit(as.numeric(sm$CFI[1]),            digits)
    holder$TLI[i]   <- .fmt_fit(as.numeric(sm$TLI[1]),            digits)
    holder$RMSEA[i] <- .fmt_fit(as.numeric(sm$RMSEA_Estimate[1]), digits)
    holder$SRMR[i]  <- .fmt_fit(as.numeric(sm$SRMR[1]),           digits)

    ## ICs & parameters
    holder$LL[i]         <- .fnum(as.numeric(sm$LL[1]), digits)
    holder$Parameters[i] <- as.integer(sm$Parameters[1])
    holder$AIC[i]        <- .fnum(as.numeric(sm$AIC[1]), digits)
    holder$BIC[i]        <- .fnum(as.numeric(sm$BIC[1]), digits)
  }

  ## --------------- differences -------------------------------------------
  if (diffTest && n > 1) {
    delta_CFI   <- delta_RMSEA <- delta_AIC <- delta_BIC <-
      delta_chisq <- delta_df <- p_delta_chisq <- rep(NA_real_, n)

    for (i in 2:n) {
      ref_i <- if (ref == "first") 1 else i - 1
      s1    <- as.data.frame(modelList[[i]]$summaries)
      s0    <- as.data.frame(modelList[[ref_i]]$summaries)

      delta_CFI[i]   <- as.numeric(s1$CFI[1])            - as.numeric(s0$CFI[1])
      delta_RMSEA[i] <- as.numeric(s1$RMSEA_Estimate[1]) - as.numeric(s0$RMSEA_Estimate[1])
      delta_AIC[i]   <- as.numeric(s1$AIC[1])            - as.numeric(s0$AIC[1])
      delta_BIC[i]   <- as.numeric(s1$BIC[1])            - as.numeric(s0$BIC[1])

      if (!is.na(s1$ChiSqM_Value[1]) && !is.na(s0$ChiSqM_Value[1])) {
        delta_chisq[i] <- as.numeric(s1$ChiSqM_Value[1]) - as.numeric(s0$ChiSqM_Value[1])
        delta_df[i]    <- as.numeric(s1$ChiSqM_DF[1])    - as.numeric(s0$ChiSqM_DF[1])
        if (!is.na(delta_df[i]) && delta_df[i] > 0)
          p_delta_chisq[i] <- 1 - pchisq(delta_chisq[i], delta_df[i])
      }
    }

    holder$`ΔCFI`   <- .fmt_diff(delta_CFI,   digits)
    holder$`ΔRMSEA` <- .fmt_diff(delta_RMSEA, digits)
    holder$`ΔAIC`   <- .fmt_diff(delta_AIC,   digits)
    holder$`ΔBIC`   <- .fmt_diff(delta_BIC,   digits)

    stars <- vapply(p_delta_chisq, function(p)
      if (is.na(p)) "" else if (p < .001) "***"
      else if (p < .01) "**" else if (p < .05) "*" else if (p < .1) "." else "",
      character(1))

    holder$`Δχ²`  <- ifelse(!is.na(delta_chisq) & !is.na(delta_df),
                            paste0(.fnum(delta_chisq, digits),
                                   " (", delta_df, ") ", stars), NA)
    holder$`p_Δχ²` <- .fmt_p(p_delta_chisq)  # always 3‑digit output

    holder[1, c("Δχ²","p_Δχ²","ΔCFI","ΔRMSEA","ΔAIC","ΔBIC")] <- "—"
  }

  ## --------------- column selection --------------------------------------
  base_idx <- c("Models","χ²(df)","χ²/df","CFI","TLI","RMSEA","SRMR")
  diff_idx <- if (diffTest && n > 1)
    c("Δχ²","p_Δχ²","ΔCFI","ΔRMSEA") else character(0)

  add_idx  <- indices
  if ("AIC" %in% add_idx) diff_idx <- union(diff_idx, "ΔAIC")
  if ("BIC" %in% add_idx) diff_idx <- union(diff_idx, "ΔBIC")

  holder[, unique(c(base_idx, add_idx, diff_idx)), drop = FALSE]
}


