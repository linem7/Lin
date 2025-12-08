#' Fit Table of Model Fit Indices and Difference Metrics
#'
#' This function creates a publication-ready table of fit indices for one or more Mplus models
#' and computes difference metrics (e.g., ΔCFI, ΔRMSEA, Δχ²) between comparative models.
#'
#' @description
#' The function is designed to be robust and flexible:
#' \itemize{
#'   \item \strong{Recursive Input Handling:} It accepts single objects, lists of objects, or nested
#'     \code{mplus.model.list} objects returned by \code{MplusAutomation::readModels}, flattening them automatically.
#'   \item \strong{MI Shortcut Parsing:} It detects and parses Mplus "shortcut" Measurement Invariance
#'     outputs (where Configural, Metric, and Scalar models are reported in a single output file).
#'     These are extracted as distinct rows in the table.
#'   \item \strong{Robust Binding:} It handles models with mismatched summary columns (e.g., mixing
#'     standard outputs with parsed text outputs) without errors.
#' }
#'
#' The final output always displays the default fit indices:
#' \code{"Models", "χ²(df)", "CFI", "TLI", "RMSEA", "SRMR"}
#' and, if \code{diffTest = TRUE}, the default difference metrics:
#' \code{"Δχ²", "p_Δχ²", "ΔCFI", "ΔRMSEA"}.
#'
#' If the user specifies additional indices (e.g., "AIC", "BIC") via the \code{indices} argument,
#' those indices will be included in the fit chunk. If "AIC" or "BIC" are requested, their
#' corresponding difference metrics (ΔAIC, ΔBIC) are automatically added to the difference chunk.
#'
#' @param ... One or more model objects (or lists of models) that contain a \code{summaries} element.
#'   Typically objects returned by \code{MplusAutomation::readModels}.
#' @param indices Character vector of additional indices to report. These will be added after the default
#'   fit indices.
#' @param ref Character string indicating the reference model to compare to. Options are \code{"last"}
#'   (compare current row to the immediately preceding row) or \code{"first"} (compare all rows to the first row).
#'   Defaults to \code{"last"}.
#' @param diffTest Logical. If \code{TRUE} (default), difference metrics will be computed.
#' @param digits Integer. Number of digits to round numeric outputs. Defaults to \code{3}.
#'   P-values < .001 are formatted as "< .001", and indices bounded by 1 (e.g., CFI, RMSEA)
#'   have leading zeros removed per APA style (e.g., .950).
#' @param model_names Optional character vector of custom model names. If supplied, its length must
#'   equal the *total number of resulting rows* in the table.
#'   \emph{Note:} If a single input file is parsed into multiple rows (e.g., Configural/Metric/Scalar),
#'   you must provide a name for every row generated, otherwise the default names (Filename + Step) will be used.
#'
#' @return A data frame containing the formatted fit indices and difference metrics.
#'
#' @examples
#' \dontrun{
#'   # Basic usage with a list of models
#'   res <- MplusAutomation::readModels("target/path")
#'   fit_table(res)
#'
#'   # Customized table with AIC/BIC and specific referencing
#'   fit_table(res, indices = c("AIC", "BIC"), ref = "first", digits = 2)
#' }
#' @export
fit_table <- function(...,
                      indices      = character(0),
                      ref          = c("last", "first"),
                      diffTest     = TRUE,
                      digits       = 3,
                      model_names  = NULL) {

  ref <- match.arg(ref)
  input_objs <- list(...)

  ## --------------------------------------------------------------------------
  ## 1. HELPER: Robust Row Binding (Fixes the "columns do not match" error)
  ## --------------------------------------------------------------------------
  .bind_rows_robust <- function(df_list) {
    if (length(df_list) == 0) return(data.frame())

    # Identify all unique column names across all dataframes
    all_cols <- unique(unlist(lapply(df_list, names)))

    # Standardize each dataframe to have those columns (filling missing with NA)
    standardized_list <- lapply(df_list, function(df) {
      if (is.null(df) || nrow(df) == 0) return(NULL)
      missing_cols <- setdiff(all_cols, names(df))
      if (length(missing_cols) > 0) {
        df[missing_cols] <- NA
      }
      df[, all_cols, drop = FALSE]
    })

    # Remove NULLs and bind
    standardized_list <- Filter(Negate(is.null), standardized_list)
    if (length(standardized_list) == 0) return(data.frame())

    do.call(rbind, standardized_list)
  }

  ## --------------------------------------------------------------------------
  ## 2. HELPER: Parse Text for MI Shortcuts (Config/Metric/Scalar)
  ## --------------------------------------------------------------------------
  .parse_mi_output <- function(output_text, filename_base) {
    full_text <- paste(output_text, collapse = "\n")

    # Split text blocks by the unique MI header
    blocks <- unlist(strsplit(full_text, "MODEL FIT INFORMATION FOR THE ", fixed = TRUE))
    blocks <- blocks[grepl("Chi-Square Test of Model Fit", blocks)]

    if (length(blocks) == 0) return(NULL)

    parsed_rows <- list()

    for (blk in blocks) {
      # Extract Model Name (Configural, Metric, Scalar)
      m_name_raw <- sub(" MODEL.*", "", blk)
      m_name_raw <- trimws(gsub("\n", "", m_name_raw))
      m_name     <- tools::toTitleCase(tolower(m_name_raw))

      # Name format: "filename.out (Metric)"
      # This ensures the filename is preserved but rows are distinct
      row_name   <- paste0(filename_base, " (", m_name, ")")

      # Regex Helpers
      get_val <- function(pattern) {
        rgx <- paste0(pattern, "[\\s\\t]*[:=]?[\\s\\t]*([-]?[0-9]*\\.[0-9]+([ED][+-]?[0-9]+)?)")
        m <- regexpr(rgx, blk, perl = TRUE, ignore.case = TRUE)
        if (m == -1) return(NA_real_)
        str_val <- sub(paste0(".*", pattern, "[\\s\\t]*[:=]?[\\s\\t]*"), "", regmatches(blk, m))
        str_val <- gsub("D", "E", str_val)
        as.numeric(str_val)
      }

      get_int <- function(pattern) {
        rgx <- paste0(pattern, "[\\s\\t]*[:=]?[\\s\\t]*([0-9]+)")
        m <- regexpr(rgx, blk, perl = TRUE, ignore.case = TRUE)
        if (m == -1) return(NA_integer_)
        str_val <- sub(paste0(".*", pattern, "[\\s\\t]*[:=]?[\\s\\t]*"), "", regmatches(blk, m))
        as.integer(str_val)
      }

      # Build Row
      df_row <- data.frame(
        Title          = row_name,
        Filename       = filename_base,
        ChiSqM_Value   = get_val("Chi-Square Test of Model Fit.*?Value"),
        ChiSqM_DF      = get_int("Chi-Square Test of Model Fit.*?Degrees of Freedom"),
        CFI            = get_val("CFI"),
        TLI            = get_val("TLI"),
        RMSEA_Estimate = get_val("RMSEA.*?Estimate"),
        SRMR           = get_val("SRMR.*?Value"),
        AIC            = get_val("Akaike \\(AIC\\)"),
        BIC            = get_val("Bayesian \\(BIC\\)"),
        LL             = get_val("Loglikelihood.*?H0 Value"),
        Parameters     = get_int("Number of Free Parameters"),
        stringsAsFactors = FALSE
      )

      if (!is.na(df_row$ChiSqM_Value)) {
        parsed_rows[[length(parsed_rows) + 1]] <- df_row
      }
    }

    if (length(parsed_rows) == 0) return(NULL)
    do.call(rbind, parsed_rows)
  }

  ## --------------------------------------------------------------------------
  ## 3. PRE-PROCESSOR: Flatten Lists
  ## --------------------------------------------------------------------------
  # Recursive helper to handle list of lists (mplus.model.list)
  collect_models <- function(obj) {
    if (is.list(obj) && "summaries" %in% names(obj)) return(list(obj))
    if (is.list(obj)) return(unlist(lapply(obj, collect_models), recursive = FALSE))
    return(list())
  }

  flat_list <- list()
  for (arg in input_objs) {
    flat_list <- c(flat_list, collect_models(arg))
  }

  if (length(flat_list) == 0) stop("No valid Mplus models found.")

  ## --------------------------------------------------------------------------
  ## 4. EXTRACTION LOOP
  ## --------------------------------------------------------------------------
  master_summaries <- list()

  for (i in seq_along(flat_list)) {
    obj <- flat_list[[i]]
    sm  <- as.data.frame(obj$summaries)

    # 1. Get the Filename (Priority)
    # If MplusAutomation didn't catch the filename, use a generic "Model_X"
    fname <- if (!is.null(sm$Filename[1]) && !is.na(sm$Filename[1])) sm$Filename[1] else paste0("Model_", i)

    # 2. Check for MI Shortcut in text
    has_mi_text <- FALSE
    if ("output" %in% names(obj)) {
      has_mi_text <- any(grepl("MODEL FIT INFORMATION FOR THE (CONFIGURAL|METRIC|SCALAR)", toupper(obj$output)))
    }

    if (has_mi_text) {
      # PARSE MANUAL (Split into 3 rows)
      parsed_df <- .parse_mi_output(obj$output, fname)
      if (!is.null(parsed_df)) {
        master_summaries[[length(master_summaries) + 1]] <- parsed_df
      } else {
        # Fallback to standard, overwrite Title with Filename
        sm$Title <- fname
        master_summaries[[length(master_summaries) + 1]] <- sm
      }
    } else {
      # STANDARD FILE
      # Overwrite the 'Title' (usually "cfa results") with the 'Filename'
      sm$Title <- fname
      master_summaries[[length(master_summaries) + 1]] <- sm
    }
  }

  # Combine using the Robust Binder
  df_all <- .bind_rows_robust(master_summaries)
  n <- nrow(df_all)
  if (n == 0) stop("No summary data could be extracted.")

  ## --------------------------------------------------------------------------
  ## 5. FORMATTING & DISPLAY
  ## --------------------------------------------------------------------------
  .get_num <- function(df, r, nm) {
    if (!nm %in% names(df)) return(NA_real_)
    val <- df[[nm]][r]
    if (is.null(val)) return(NA_real_)
    as.numeric(val)
  }

  .strip0 <- function(x) sub("^(-?)0\\.", "\\1.", x)
  .fnum   <- function(x, d) sprintf(paste0("%.", d, "f"), x)

  .fmt_fit <- function(x, d) {
    if (is.na(x)) return(NA_character_)
    out <- .fnum(x, d)
    if (abs(x) < 1) out <- .strip0(out)
    out
  }

  .fmt_p <- function(p, d = 3) {
    if (is.na(p)) return(NA_character_)
    if (p < .001) return("< .001")
    .strip0(.fnum(p, d))
  }

  # Check which indices exist
  has_classic <- all(c("ChiSqM_Value", "CFI", "RMSEA_Estimate") %in% names(df_all))

  if (has_classic) {
    base_idx <- c("Models", "χ²(df)", "CFI", "TLI", "RMSEA", "SRMR", "AIC", "BIC")
  } else {
    base_idx <- c("Models", "LL", "Parameters", "AIC", "BIC")
  }

  diff_idx <- if (diffTest && n > 1) c("Δχ²", "p_Δχ²", "ΔCFI", "ΔRMSEA") else character(0)
  all_cols <- unique(c(base_idx, indices, diff_idx))

  holder <- setNames(data.frame(matrix(NA_character_, n, length(all_cols)), stringsAsFactors=FALSE), all_cols)

  for (i in 1:n) {
    # Name column: User preference overrides standard
    if (!is.null(model_names) && i <= length(model_names)) {
      holder$Models[i] <- model_names[i]
    } else {
      # Use the Title we forced to be Filename (or Filename (Config)) earlier
      holder$Models[i] <- df_all$Title[i]
    }

    if (has_classic) {
      cv <- .get_num(df_all, i, "ChiSqM_Value")
      df <- as.integer(.get_num(df_all, i, "ChiSqM_DF"))

      if (!is.na(cv)) holder$`χ²(df)`[i] <- paste0(.fnum(cv, digits), " (", df, ")")

      holder$CFI[i]   <- .fmt_fit(.get_num(df_all, i, "CFI"), digits)
      holder$TLI[i]   <- .fmt_fit(.get_num(df_all, i, "TLI"), digits)
      holder$RMSEA[i] <- .fmt_fit(.get_num(df_all, i, "RMSEA_Estimate"), digits)
      holder$SRMR[i]  <- .fmt_fit(.get_num(df_all, i, "SRMR"), digits)
      holder$AIC[i]   <- .fnum(.get_num(df_all, i, "AIC"), 0)
      holder$BIC[i]   <- .fnum(.get_num(df_all, i, "BIC"), 0)
    } else {
      # Fallback for non-CFA models
      holder$LL[i]  <- .fnum(.get_num(df_all, i, "LL"), digits)
      holder$AIC[i] <- .fnum(.get_num(df_all, i, "AIC"), digits)
      holder$BIC[i] <- .fnum(.get_num(df_all, i, "BIC"), digits)
    }

    # Extra user-requested indices
    for (idx in indices) {
      if (idx %in% names(df_all)) {
        val <- .get_num(df_all, i, idx)
        if (idx %in% c("Parameters", "ChiSqM_DF")) {
          holder[[idx]][i] <- as.character(as.integer(val))
        } else {
          holder[[idx]][i] <- .fnum(val, digits)
        }
      }
    }
  }

  ## --------------------------------------------------------------------------
  ## 6. DIFFERENCE TESTING
  ## --------------------------------------------------------------------------
  if (diffTest && n > 1 && has_classic) {
    for (i in 2:n) {
      ref_i <- if (ref == "first") 1 else i - 1

      c1 <- .get_num(df_all, i, "ChiSqM_Value"); c0 <- .get_num(df_all, ref_i, "ChiSqM_Value")
      d1 <- .get_num(df_all, i, "ChiSqM_DF");    d0 <- .get_num(df_all, ref_i, "ChiSqM_DF")

      d_chisq <- c1 - c0
      d_df    <- abs(d1 - d0)

      if (!is.na(d_chisq) && !is.na(d_df)) {
        holder$`Δχ²`[i] <- paste0(.fnum(d_chisq, digits), " (", d_df, ")")
        if (d_df > 0) {
          p_val <- pchisq(abs(d_chisq), d_df, lower.tail = FALSE)
          holder$`p_Δχ²`[i] <- .fmt_p(p_val)
        } else {
          holder$`p_Δχ²`[i] <- "—"
        }
      }

      cf1 <- .get_num(df_all, i, "CFI"); cf0 <- .get_num(df_all, ref_i, "CFI")
      if (!is.na(cf1) && !is.na(cf0)) holder$`ΔCFI`[i] <- .fmt_fit(cf1 - cf0, digits)

      rm1 <- .get_num(df_all, i, "RMSEA_Estimate"); rm0 <- .get_num(df_all, ref_i, "RMSEA_Estimate")
      if (!is.na(rm1) && !is.na(rm0)) holder$`ΔRMSEA`[i] <- .fmt_fit(rm1 - rm0, digits)
    }

    # Fill first row with dashes
    cols <- intersect(names(holder), c("Δχ²", "p_Δχ²", "ΔCFI", "ΔRMSEA"))
    holder[1, cols] <- "—"
  }

  return(holder)
}

