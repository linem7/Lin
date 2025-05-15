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
                      indices = character(0),  # additional indices to report
                      ref = c("last", "first"),
                      diffTest = TRUE,
                      digits = 3,
                      model_names = NULL) {

  ref <- match.arg(ref)

  ## --- STEP 0: Flatten input objects into a single list of models ---
  args <- list(...)
  modelList <- list()
  for (obj in args) {
    if (is.list(obj) && !("summaries" %in% names(obj))) {
      for (m in obj) {
        if (!("summaries" %in% names(m))) {
          stop("One of the models in the provided list does not contain a 'summaries' element.")
        }
        modelList[[length(modelList) + 1]] <- m
      }
    } else if (is.list(obj) && ("summaries" %in% names(obj))) {
      modelList[[length(modelList) + 1]] <- obj
    } else {
      stop("Input is not a valid model object or list of models.")
    }
  }

  ## --- STEP1: Validation for each arguments
  n <- length(modelList)
  if(n == 0) stop("No models provided.")
  if (!is.character(indices))
    stop("`indices` must be a character vector, extra index including LL, Parameters, AIC, and BIC.")
  if (!is.logical(diffTest) || length(diffTest) != 1)
    stop("`diffTest` must be a single TRUE or FALSE.")
  if (!is.numeric(digits) || length(digits) != 1 || digits < 0)
    stop("`digits` must be a single non‐negative number.")

  if (!is.null(model_names)) {
    if (!is.character(model_names))
      stop("`model_names` must be a character vector.")
    if (length(model_names) != n)
      stop("`model_names` must have length equal to the number of models (", n, ").")
  }

  #––– Vectorized helper –––
  .format_value <- function(x, digits) {
    x_num    <- as.numeric(x)
    thresh   <- 10^(-digits)
    thresh_s <- sprintf(paste0("%.", digits, "f"), thresh)
    vapply(x_num, function(xi) {
      if (is.na(xi)) return(NA_character_)
      if (abs(xi) < thresh) paste0("<", thresh_s)
      else                 sprintf(paste0("%.", digits, "f"), xi)
    }, FUN.VALUE = character(1), USE.NAMES = FALSE)
  }

  ## --- STEP 2: Create an empty data.frame holder with default columns ---
  # Note: "Parameters" follows "LL".
  default_cols <- c("Models", "χ²(df)", "χ²/df", "CFI", "TLI", "RMSEA", "SRMR",
                    "LL", "Parameters", "AIC", "BIC")
  holder <- data.frame(matrix(NA, nrow = n, ncol = length(default_cols)),
                       stringsAsFactors = FALSE)
  colnames(holder) <- default_cols

  ## --- STEP 3: Populate the holder for each model ---

  for(i in seq_len(n)) {
    summ <- as.data.frame(modelList[[i]][["summaries"]])

    # Model name: prefer model_names over Filename
    holder$Models[i] <- if (!is.null(model_names)) {
      model_names[i]
    } else if ("Filename" %in% names(summ)) {
      summ$Filename[1]
    } else {
      NA
    }

    if("ChiSqM_Value" %in% names(summ)) {
      chisq_val <- sprintf(paste0("%.", digits, "f"), as.numeric(summ$ChiSqM_Value[1]))
    } else {
      chisq_val <- NA
    }

    df_val <- if("ChiSqM_DF" %in% names(summ)) as.numeric(summ$ChiSqM_DF[1]) else NA

    # Combine chi-square and df into a single column (χ²(df))
    holder$`χ²(df)`[i] <- if (!is.na(chisq_val) && !is.na(df_val)) {
      paste0(chisq_val, " (", df_val, ")")
    } else {
      NA
    }

    # Compute χ²/df if both chi-square and df are available.
    {
      chisq_raw <- if("ChiSqM_Value" %in% names(summ)) summ$ChiSqM_Value[1] else NA
      df_raw    <- if("ChiSqM_DF" %in% names(summ)) summ$ChiSqM_DF[1] else NA

      chisq_val <- suppressWarnings(as.numeric(as.character(chisq_raw)))
      df_val    <- suppressWarnings(as.numeric(as.character(df_raw)))

      if (!is.na(chisq_val) && !is.na(df_val)) {
        # Handle zero degrees of freedom
        if (df_val == 0) {
          # If both chisq and df are zero, ratio is 0; otherwise, Inf.
          ratio <- Inf
        } else {
          ratio <- chisq_val / df_val
        }
        holder[["χ²/df"]][i] <- .format_value(as.numeric(ratio), digits)
      } else {
        holder[["χ²/df"]][i] <- NA
      }
    }

    # Other fit indices
    holder$CFI[i]   <- if("CFI" %in% names(summ)) sprintf(paste0("%.", digits, "f"), as.numeric(summ$CFI[1])) else NA
    holder$TLI[i]   <- if("TLI" %in% names(summ)) sprintf(paste0("%.", digits, "f"), as.numeric(summ$TLI[1])) else NA
    holder$RMSEA[i] <- if("RMSEA_Estimate" %in% names(summ)) sprintf(paste0("%.", digits, "f"), as.numeric(summ$RMSEA_Estimate[1])) else NA
    holder$SRMR[i]  <- if("SRMR" %in% names(summ)) sprintf(paste0("%.", digits, "f"), as.numeric(summ$SRMR[1])) else NA

    # Information criteria and parameters (Parameters follows LL)
    holder$LL[i]         <- if("LL" %in% names(summ)) sprintf(paste0("%.", digits, "f"), as.numeric(summ$LL[1])) else NA
    holder$Parameters[i] <- if("Parameters" %in% names(summ)) as.numeric(summ$Parameters[1]) else NA
    holder$AIC[i]        <- if("AIC" %in% names(summ)) sprintf(paste0("%.", digits, "f"), as.numeric(summ$AIC[1])) else NA
    holder$BIC[i]        <- if("BIC" %in% names(summ)) sprintf(paste0("%.", digits, "f"), as.numeric(summ$BIC[1])) else NA
  }


  ## --- STEP 4: Compute difference metrics (if diffTest is TRUE) ---
  if(diffTest && n > 1) {
    # Preallocate difference vectors.
    ΔCFI      <- rep(NA, n)
    ΔRMSEA    <- rep(NA, n)
    ΔAIC      <- rep(NA, n)
    ΔBIC      <- rep(NA, n)
    `Δχ²`   <- rep(NA, n)
    Δdf      <- rep(NA, n)
    `p_Δχ²` <- rep(NA, n)

    for(i in 2:n) {
      # Choose reference row based on 'ref'
      ref_row <- if(ref == "first") 1 else i - 1

      # Retrieve raw summaries for current and reference models.
      summ_current <- as.data.frame(modelList[[i]][["summaries"]])
      summ_ref     <- as.data.frame(modelList[[ref_row]][["summaries"]])

      # Compute differences for CFI, RMSEA, AIC, and BIC using raw numeric values.
      if("CFI" %in% names(summ_current) && "CFI" %in% names(summ_ref))
        ΔCFI[i]   <- as.numeric(summ_current$CFI[1]) - as.numeric(summ_ref$CFI[1])
      if("RMSEA_Estimate" %in% names(summ_current) && "RMSEA_Estimate" %in% names(summ_ref))
        ΔRMSEA[i] <- as.numeric(summ_current$RMSEA_Estimate[1]) - as.numeric(summ_ref$RMSEA_Estimate[1])
      if("AIC" %in% names(summ_current) && "AIC" %in% names(summ_ref))
        ΔAIC[i]   <- as.numeric(summ_current$AIC[1]) - as.numeric(summ_ref$AIC[1])
      if("BIC" %in% names(summ_current) && "BIC" %in% names(summ_ref))
        ΔBIC[i]   <- as.numeric(summ_current$BIC[1]) - as.numeric(summ_ref$BIC[1])

      ## --- Chi-square difference ---
      # Use chi-square and df if available; otherwise, use LL-based approach.
      if("ChiSqM_Value" %in% names(summ_current) && "ChiSqM_Value" %in% names(summ_ref) &&
         "ChiSqM_DF" %in% names(summ_current) && "ChiSqM_DF" %in% names(summ_ref)) {
        delta_chisq  <- as.numeric(summ_current$ChiSqM_Value[1]) - as.numeric(summ_ref$ChiSqM_Value[1])
        delta_df_val <- as.numeric(summ_current$ChiSqM_DF[1]) - as.numeric(summ_ref$ChiSqM_DF[1])
      } else if("LL" %in% names(summ_current) && "LL" %in% names(summ_ref) &&
                "Parameters" %in% names(summ_current) && "Parameters" %in% names(summ_ref)) {
        if(!is.null(summ_current[["LLCorrectionFactor"]]) &&
           !is.null(summ_ref[["LLCorrectionFactor"]])) {
          cf0 <- as.numeric(summ_ref$LLCorrectionFactor[1])
          cf1 <- as.numeric(summ_current$LLCorrectionFactor[1])
          cd <- (as.numeric(summ_ref$Parameters[1]) * cf0 - as.numeric(summ_current$Parameters[1]) * cf1) /
            (as.numeric(summ_ref$Parameters[1]) - as.numeric(summ_current$Parameters[1]))
          delta_chisq <- -2 * (as.numeric(summ_ref$LL[1]) - as.numeric(summ_current$LL[1])) / cd
        } else {
          delta_chisq <- -2 * (as.numeric(summ_ref$LL[1]) - as.numeric(summ_current$LL[1]))
        }
        delta_df_val <- as.numeric(summ_current$Parameters[1]) - as.numeric(summ_ref$Parameters[1])
      } else {
        delta_chisq  <- NA
        delta_df_val <- NA
      }

      `Δχ²`[i] <- delta_chisq
      Δdf[i]    <- delta_df_val
      p_val <- if(!is.na(delta_df_val) && delta_df_val > 0 && !is.na(delta_chisq))
        1 - pchisq(delta_chisq, df = delta_df_val) else NA
      `p_Δχ²`[i] <- p_val
    }

    # Append formatted difference metrics to the holder.
    holder[["ΔCFI"]]   <- .format_value(ΔCFI,   digits)
    holder[["ΔRMSEA"]] <- .format_value(ΔRMSEA, digits)
    if("AIC" %in% colnames(holder)) {
      holder[["ΔAIC"]] <- .format_value(ΔAIC, digits)
    }
    if("BIC" %in% colnames(holder)) {
      holder[["ΔBIC"]] <- .format_value(ΔBIC, digits)
    }

    stars <- vapply(`p_Δχ²`, function(p) {
      if      (is.na(p))      ""
      else if (p < 0.001)     "***"
      else if (p < 0.01)      "**"
      else if (p < 0.05)      "*"
      else if (p < 0.1)       "."
      else                    ""
    }, FUN.VALUE = character(1))

    holder[["Δχ²"]] <- ifelse(
      !is.na(`Δχ²`) & !is.na(Δdf),
      paste0(.format_value(`Δχ²`, digits), " (", Δdf, ") ", stars),
      NA_character_
    )
    holder[["p_Δχ²"]] <- ifelse(!is.na(`p_Δχ²`), .format_value(as.numeric(`p_Δχ²`), digits), NA)
  }

  if (diffTest && n > 1) {
    diff_cols <- c("Δχ²", "p_Δχ²", "ΔCFI", "ΔRMSEA",
                   intersect(colnames(holder), c("ΔAIC", "ΔBIC")))
    holder[1, diff_cols] <- "\u2013"
  }

  ## --- STEP 5: Select and rearrange columns based on the 'indices' argument ---
  # Default indices that are always reported.
  default_indices <- c("Models", "χ²(df)", "χ²/df", "CFI", "TLI", "RMSEA", "SRMR")
  # Default difference metrics to include.
  default_diff    <- if (diffTest && n > 1) {
    c("Δχ²","p_Δχ²","ΔCFI","ΔRMSEA")
  } else {
    character(0)
  }

  if(length(indices) > 0) {
    # If "AIC" is mentioned, add both AIC and its difference; similarly for "BIC".
    extra_diff <- c()
    if("AIC" %in% indices) extra_diff <- c(extra_diff, "ΔAIC")
    if("BIC" %in% indices) extra_diff <- c(extra_diff, "ΔBIC")
    fit_cols <- union(default_indices, indices)
    change_cols <- union(default_diff, extra_diff)
    final_cols <- union(fit_cols, change_cols)
    final_holder <- holder[, final_cols, drop = FALSE]
  } else {
    # If no indices provided, show default fit indices and default difference metrics.
    final_holder <- holder[, unique(c(default_indices, default_diff)), drop = FALSE]
  }

  rownames(final_holder) <- NULL
  return(final_holder)
}

