#' Fit Table of Model Fit Indices and Difference Metrics
#'
#' This function creates a table of fit indices for one or more models and computes
#' difference metrics (e.g., ΔCFI, ΔRMSEA, ΔChisq) between comparative models. If the
#' primary chi-square values or degrees of freedom are unavailable, the function will
#' automatically switch to a log-likelihood (LL) based approach. When available, the
#' LLCorrectionFactor is used; otherwise, the LL difference is computed using LL alone.
#'
#' The final output always displays the default fit indices:
#' \code{"Models", "Chisq", "df", "chisq/df", "CFI", "TLI", "RMSEA", "SRMR"}
#' and the default difference metrics:
#' \code{"ΔChisq", "p_ΔChisq", "ΔCFI", "ΔRMSEA"}. If the user specifies additional indices
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
                      digits = 3) {

  ref <- match.arg(ref)

  ## --- STEP 1: Flatten input objects into a single list of models ---
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

  n <- length(modelList)
  if(n == 0) stop("No models provided.")

  ## --- STEP 2: Create an empty data.frame holder with default columns ---
  # Note: "Parameters" follows "LL".
  default_cols <- c("Models", "Chisq", "df", "chisq/df", "CFI", "TLI", "RMSEA", "SRMR",
                    "LL", "Parameters", "AIC", "BIC")
  holder <- data.frame(matrix(NA, nrow = n, ncol = length(default_cols)),
                       stringsAsFactors = FALSE)
  colnames(holder) <- default_cols

  ## --- STEP 3: Populate the holder for each model ---
  for(i in seq_len(n)) {
    summ <- as.data.frame(modelList[[i]][["summaries"]])

    # Model name: if available, from 'Filename'
    holder$Models[i] <- if("Filename" %in% names(summ)) summ$Filename[1] else NA

    # Primary chi-square and df values (if available)
    holder$Chisq[i] <- if("ChiSqM_Value" %in% names(summ)) round(as.numeric(summ$ChiSqM_Value[1]), digits) else NA
    holder$df[i]    <- if("ChiSqM_DF"    %in% names(summ)) as.numeric(summ$ChiSqM_DF[1]) else NA

    # Compute chisq/df if both chi-square and df are available and df != 0.
    if(!is.na(holder$Chisq[i]) && !is.na(holder$df[i]) && holder$df[i] != 0) {
      holder[["chisq/df"]][i] <- holder$Chisq[i] / holder$df[i]
    } else {
      holder[["chisq/df"]][i] <- NA
    }

    # Other fit indices
    holder$CFI[i]   <- if("CFI"            %in% names(summ)) round(as.numeric(summ$CFI[1]), digits) else NA
    holder$TLI[i]   <- if("TLI"            %in% names(summ)) round(as.numeric(summ$TLI[1]), digits) else NA
    holder$RMSEA[i] <- if("RMSEA_Estimate" %in% names(summ)) round(as.numeric(summ$RMSEA_Estimate[1]), digits) else NA
    holder$SRMR[i]  <- if("SRMR"           %in% names(summ)) round(as.numeric(summ$SRMR[1]), digits) else NA

    # Information criteria and parameters (Parameters follows LL)
    holder$LL[i]         <- if("LL"         %in% names(summ)) round(as.numeric(summ$LL[1]), digits) else NA
    holder$Parameters[i] <- if("Parameters" %in% names(summ)) as.numeric(summ$Parameters[1]) else NA
    holder$AIC[i]        <- if("AIC"        %in% names(summ)) round(as.numeric(summ$AIC[1]), digits) else NA
    holder$BIC[i]        <- if("BIC"        %in% names(summ)) round(as.numeric(summ$BIC[1]), digits) else NA
  }

  ## --- STEP 4: Compute difference metrics (if diffTest is TRUE) ---
  if(diffTest && n > 1) {
    # Preallocate difference vectors.
    ΔCFI      <- rep(NA, n)
    ΔRMSEA    <- rep(NA, n)
    ΔAIC      <- rep(NA, n)
    ΔBIC      <- rep(NA, n)
    ΔChisq   <- rep(NA, n)
    Δdf      <- rep(NA, n)
    p_ΔChisq <- rep(NA, n)

    for(i in 2:n) {
      # Choose reference row based on 'ref'
      ref_row <- if(ref == "first") 1 else i - 1

      # Compute differences for CFI, RMSEA, AIC, and BIC.
      ΔCFI[i]   <- holder$CFI[i]   - holder$CFI[ref_row]
      ΔRMSEA[i] <- holder$RMSEA[i] - holder$RMSEA[ref_row]
      ΔAIC[i]   <- holder$AIC[i]   - holder$AIC[ref_row]
      ΔBIC[i]   <- holder$BIC[i]   - holder$BIC[ref_row]

      ## --- Chi-square difference ---
      # Use chi-square and df if available; otherwise, use LL-based approach.
      if(!is.na(holder$Chisq[i]) && !is.na(holder$Chisq[ref_row]) &&
         !is.na(holder$df[i]) && !is.na(holder$df[ref_row])) {
        delta_chisq  <- holder$Chisq[i] - holder$Chisq[ref_row]
        delta_df_val <- holder$df[i] - holder$df[ref_row]
      } else if(!is.na(holder$LL[i]) && !is.na(holder$LL[ref_row]) &&
                !is.na(holder$Parameters[i]) && !is.na(holder$Parameters[ref_row])) {
        # Use LLCorrectionFactor if available; otherwise, use LL alone.
        if(!is.null(modelList[[i]][["summaries"]][["LLCorrectionFactor"]]) &&
           !is.null(modelList[[ref_row]][["summaries"]][["LLCorrectionFactor"]])) {
          cf0 <- as.numeric(modelList[[ref_row]][["summaries"]][["LLCorrectionFactor"]][1])
          cf1 <- as.numeric(modelList[[i]][["summaries"]][["LLCorrectionFactor"]][1])
          cd <- (holder$Parameters[ref_row] * cf0 - holder$Parameters[i] * cf1) /
            (holder$Parameters[ref_row] - holder$Parameters[i])
          delta_chisq <- -2 * (holder$LL[ref_row] - holder$LL[i]) / cd
        } else {
          delta_chisq <- -2 * (holder$LL[ref_row] - holder$LL[i])
        }
        delta_df_val <- holder$Parameters[i] - holder$Parameters[ref_row]
      } else {
        delta_chisq  <- NA
        delta_df_val <- NA
      }

      ΔChisq[i] <- delta_chisq
      Δdf[i]    <- delta_df_val
      p_ΔChisq[i] <- if(!is.na(delta_df_val) && delta_df_val > 0 && !is.na(delta_chisq)) {
        1 - pchisq(delta_chisq, df = delta_df_val)
      } else NA
    }

    ## Append the difference metrics to the holder.
    holder[["ΔCFI"]]   <- ifelse(!is.na(ΔCFI), round(ΔCFI, digits), NA)
    holder[["ΔRMSEA"]] <- ifelse(!is.na(ΔRMSEA), round(ΔRMSEA, digits), NA)
    holder[["ΔChisq"]] <- ifelse(!is.na(ΔChisq) & !is.na(Δdf),
                                 paste0(round(ΔChisq, digits), " (", round(Δdf, 0), ") ",
                                        sapply(p_ΔChisq, function(p) {
                                          if(is.na(p)) return("")
                                          else if(p < 0.001) return("***")
                                          else if(p < 0.01)  return("**")
                                          else if(p < 0.05)  return("*")
                                          else if(p < 0.1)   return(".")
                                          else return("")
                                        })),
                                 NA)
    holder[["p_ΔChisq"]] <- ifelse(!is.na(p_ΔChisq), round(p_ΔChisq, digits), NA)
    if("AIC" %in% colnames(holder)) {
      holder[["ΔAIC"]] <- ifelse(!is.na(ΔAIC), round(ΔAIC, digits), NA)
    }
    if("BIC" %in% colnames(holder)) {
      holder[["ΔBIC"]] <- ifelse(!is.na(ΔBIC), round(ΔBIC, digits), NA)
    }
  }

  ## --- STEP 5: Select and rearrange columns based on the 'indices' argument ---
  # Default indices that are always reported.
  default_indices <- c("Models", "Chisq", "df", "chisq/df", "CFI", "TLI", "RMSEA", "SRMR")
  # Default difference metrics to include.
  default_diff <- c("ΔChisq", "p_ΔChisq", "ΔCFI", "ΔRMSEA")

  if(length(indices) > 0) {
    # If "AIC" is mentioned, add both AIC and its difference; similarly for "BIC".
    extra_diff <- c()
    if("AIC" %in% indices) extra_diff <- c(extra_diff, "ΔAIC")
    if("BIC" %in% indices) extra_diff <- c(extra_diff, "ΔBIC")
    # Also include any other indices provided that exist in to holder.
    extra_change <- union(intersect(indices, colnames(holder)), extra_diff)
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
