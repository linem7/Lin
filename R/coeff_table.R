#' Generate a Regression‑Coefficient Table with Model‑Fit Indices
#'
#' @title Regression Coefficient Table (Mplus ⇢ HTML)
#'
#' @description
#' `coeff_table()` ingests one or more model objects produced by
#' **\pkg{MplusAutomation}** (e.g., via [readModels()] or
#' [extractModelSummaries()]) and outputs a publication‑ready HTML table that
#' combines regression coefficients and model‑fit indices.  The result is a
#' **\pkg{kableExtra}** object that can be used directly in R Markdown, Quarto,
#' or Shiny.
#'
#' @param ...          One or more Mplus model objects *or* nested lists of such
#'                     objects.
#' @param stat         Character vector of coefficient statistics to display;
#'                     defaults to `c("est", "se", "p")`.  Use `"all"` for all
#'                     four statistics (est, se, z, p).
#' @param indices      Extra fit indices to add (e.g., `"AIC"`, `"BIC"`).
#' @param digits_coeff Decimal precision for coefficients (default `3`).
#' @param digits_fit   Decimal precision for fit indices (default `3`).
#' @param standardized Logical; if `TRUE`, use `stdyx.standardized` values.
#' @param params       Which parameter classes to keep (e.g., `"regression"`).
#' @param model_names  Optional vector of display names for models.
#' @param stars        Logical; append significance stars when `TRUE`.
#' @param title,notes  Caption and foot‑note strings.
#' @param exclude      Regex patterns: rows whose *rhs* matches are dropped.
#' @param order_by     Either `"lhs"` (default) or `"rhs"` for row ordering.
#' @param annotations  Either `NULL`, a `data.frame`, or a named `list` that
#'                     defines extra annotation columns (see Details).
#' @param full_width   Passed to [kableExtra::kable_classic()].
#'
#' @return A **kableExtra** HTML table.
#'
#' @examples
#' \dontrun{
#' ## Single model
#' m1 <- MplusAutomation::readModels("cfa.out")
#' coeff_table(m1)
#'
#' ## Multiple models with labels & extra indices
#' mods <- MplusAutomation::readModels("sem_dir/")
#' coeff_table(mods, model_names = c("Base", "Mediated"), indices = c("AIC", "BIC"))
#'
#' ## With annotation columns
#' coeff_table(mods$model1, mods$model2,
#'             annotations = list(Type = c("Father→Father", rep("", 7))),
#'             stars = TRUE)
#' }
#'
#' @import knitr
#' @importFrom dplyr bind_rows full_join
#' @importFrom magrittr %>%
#' @export
coeff_table <- function(...,
                        stat = c("est", "se", "p"),
                        indices = character(0),
                        digits_coeff = 3,
                        digits_fit   = 3,
                        standardized = FALSE,
                        params       = c("regression"),
                        model_names  = NULL,
                        stars        = FALSE,
                        title        = "Table X. Unstandardized Regression Coefficients and Model Fit Indices.",
                        notes        = "Unstandardized coefficients are reported.",
                        exclude      = character(0),
                        order_by     = c("lhs", "rhs"),
                        annotations  = NULL,
                        full_width   = FALSE) {
  if (length(stat) == 1 && stat == "all") stat <- c("est","se","z","p")

  add_stars <- function(est, p) {
    stars_chr <- ifelse(is.na(p), "",
                        ifelse(p < 0.001, "***",
                               ifelse(p < 0.01,  "**",
                                      ifelse(p < 0.05,  "*", ""))))
    paste0(est, stars_chr)
  }

  flatten_models <- function(x, parent_name = NULL) {
    flat <- list()
    if (is.list(x) && !is.null(x$summaries)) {
      nm <- parent_name %||% "model"
      setNames(list(x), nm)
    } else if (is.list(x)) {
      for (i in seq_along(x)) {
        nm <- names(x)[i] %||% paste0("Model", i)
        flat <- c(flat, flatten_models(x[[i]], nm))
      }
      flat
    } else {
      list()
    }
  }

  args      <- list(...)
  modelList <- unlist(lapply(args, flatten_models), recursive = FALSE)
  if (length(modelList) == 0) stop("No valid model objects were found in the input.")

  ## ── NEW ──
  ann_cols <- if (!is.null(annotations)) names(annotations) else character(0)

  num_models <- length(modelList)

  # determine model_names
  if (is.null(model_names)) {
    model_names <- vapply(modelList, function(m) {
      fn <- m$summaries$Filename[1] %||% ""
      if (nzchar(fn)) fn else NA_character_
    }, character(1))
    miss <- which(is.na(model_names))
    if (length(miss)>0) model_names[miss] <- paste0("Model ", miss)
  } else if (length(model_names) != num_models) {
    warning(sprintf("Length of model_names (%d) does not match number of models (%d); using default names.",
                    length(model_names), num_models))
    model_names <- paste0("Model ", seq_len(num_models))
  }

  model_tables <- vector("list", num_models)
  fit_list     <- vector("list", num_models)

  for (i in seq_len(num_models)) {
    mod   <- modelList[[i]]
    coeffs <- if (standardized) {
      mod$parameters$stdyx.standardized %||%
        stop(sprintf("No standardized results for model %d", i))
    } else {
      mod$parameters$unstandardized
    }
    req_cols <- c("paramHeader","param","est","se","pval","est_se")
    if (!all(req_cols %in% names(coeffs)))
      stop("Coefficient data missing columns: ", paste(setdiff(req_cols,names(coeffs)), collapse=","))
    if (!"p" %in% names(coeffs)) names(coeffs)[names(coeffs)=="pval"] <- "p"

    coeffs$ptype <- mapply(classify_param, coeffs$paramHeader, coeffs$param)
    coeffs       <- coeffs[coeffs$ptype %in% params, , drop=FALSE]

    parts <- t(mapply(create_formula_parts, coeffs$paramHeader, coeffs$param,
                      SIMPLIFY=TRUE, USE.NAMES=FALSE))
    parts <- as.data.frame(parts, stringsAsFactors=FALSE)
    coeffs <- cbind(parts, coeffs)

    # exclude rows
    if (length(exclude)>0) {
      pat <- paste(exclude, collapse="|")
      coeffs <- coeffs[!grepl(pat, coeffs$rhs, ignore.case=TRUE), ]
    }

    pvals <- as.numeric(coeffs$p)
    if ("est" %in% stat) coeffs$est <- fmt_num(as.numeric(coeffs$est), digits_coeff)
    if ("se"  %in% stat) coeffs$se  <- fmt_num(as.numeric(coeffs$se),  digits_coeff)
    if ("z"   %in% stat) coeffs$z   <- fmt_num(as.numeric(coeffs$est_se), digits_coeff)
    if ("p"   %in% stat) coeffs$p   <- ifelse(pvals<0.001, "<0.001", fmt_num(pvals, digits_coeff))
    if (stars && "est"%in%stat) coeffs$est <- add_stars(coeffs$est, pvals)

    # rename stat columns
    for (s in stat) names(coeffs)[names(coeffs)==s] <- paste0("M",i,"_",s)

    # keep only needed
    coef_cols <- paste0("M",i,"_",stat)
    tab_cols  <- c("lhs","op","rhs", coef_cols)
    if (nrow(coeffs)==0) coeffs <- setNames(data.frame(matrix(ncol=length(tab_cols), nrow=0)), tab_cols)
    else coeffs <- coeffs[, tab_cols, drop=FALSE]
    model_tables[[i]] <- coeffs

    ## --- Fit indices ----------------------------------------------------
    summ        <- mod$summaries
    default_ix  <- c("ChiSqM_Value","ChiSqM_DF","χ²/df","CFI","TLI",
                     "RMSEA_Estimate","SRMR")
    all_ix      <- unique(c(default_ix, indices))

    fit_vals <- vapply(all_ix, function(idx) {
      if (idx == "ChiSqM_DF") {
        if (!is.null(summ$ChiSqM_DF)) {
          as.character(as.integer(summ$ChiSqM_DF))
        } else {
          NA_character_
        }
      } else if (idx == "χ²/df") {
        if (!is.null(summ$ChiSqM_Value) &&
            !is.null(summ$ChiSqM_DF)   &&
            summ$ChiSqM_DF != 0) {
          fmt_num(summ$ChiSqM_Value / summ$ChiSqM_DF, digits_fit)
        } else {
          NA_character_
        }
      } else if (!is.null(summ[[idx]])) {
        fmt_num(summ[[idx]], digits_fit)
      } else {
        NA_character_
      }
    }, character(1))

    # pretty labels
    rename_map <- c(ChiSqM_Value = "χ²", ChiSqM_DF = "df",
                    RMSEA_Estimate = "RMSEA")
    labels <- ifelse(all_ix %in% names(rename_map),
                     rename_map[all_ix], all_ix)

    ## ---- NEW: build row with annotation stub ---------------------------
    n_lab  <- length(labels)
    # empty columns for every part of the table --------------------------
    fit_row <- data.frame(matrix("", nrow = n_lab,
                                 ncol = length(c(ann_cols, "lhs","op","rhs"))),
                          stringsAsFactors = FALSE)
    names(fit_row) <- c(ann_cols, "lhs", "op", "rhs")

    if (length(ann_cols) > 0) {
      # put the label in the FIRST annotation column (e.g., Type)
      fit_row[[ann_cols[1]]] <- labels
    } else {
      # fall back to lhs if no annotation cols present
      fit_row$lhs <- labels
    }

    # drop operator/rhs for fit rows
    fit_row$op  <- ""
    fit_row$rhs <- ""

    # attach values (only into the first stat column for this model) ------
    first_stat  <- paste0("M", i, "_", stat[1])
    for (s in stat) fit_row[[paste0("M", i, "_", s)]] <- ""
    fit_row[[first_stat]] <- fit_vals

    fit_list[[i]] <- fit_row
  }

  # merge coefficients
  merged_coeff <- Reduce(function(a,b) full_join(a,b,by=c("lhs","op","rhs")), model_tables)



  # order coefficients
  op_order <- c("=~"=1, "~"=2, "~ 1"=3, "~~"=4, ":="=5)
  merged_coeff$order <- op_order[merged_coeff$op]
  by <- match.arg(order_by)
  merged_coeff <- if (by=="lhs") merged_coeff[order(merged_coeff$order, merged_coeff$lhs),]
  else           merged_coeff[order(merged_coeff$order, merged_coeff$rhs),]
  merged_coeff$order <- NULL

  # --- Annotation block (updated to handle manual vectors) ---
  if (!is.null(annotations)) {
    if (is.data.frame(annotations)) {
      if (nrow(annotations) != nrow(merged_coeff))
        stop("annotation data.frame must have ", nrow(merged_coeff), " rows")
      merged_coeff <- cbind(annotations, merged_coeff)
    } else if (is.list(annotations)) {
      for (col in names(annotations)) {
        vals <- annotations[[col]]
        if (is.character(vals) && length(vals) == nrow(merged_coeff)) {
          # manual vector
          merged_coeff[[col]] <- vals
        } else if (is.character(vals) && !is.null(names(vals))) {
          # regex mapping: named vector of patterns
          merged_coeff[[col]] <- ""
          for (val in names(vals)) {
            hits <- grepl(vals[[val]], merged_coeff$rhs)
            merged_coeff[[col]][hits] <- val
          }
          # only first occurrence of each label
          merged_coeff[[col]][duplicated(merged_coeff[[col]])] <- ""
        } else {
          stop("Invalid annotations[['", col, "']]: must be a length-n character vector or named regex vector.")
        }
      }
    } else {
      stop("annotations must be either a data.frame or a named list of character vectors.")
    }
  }

  # merge fit
  merged_fit <- Reduce(function(a,b) full_join(a,b,by=c("lhs","op","rhs")), fit_list)

  # bind and clean
  final_table <- bind_rows(merged_coeff, merged_fit)
  rownames(final_table) <- NULL
  final_table[is.na(final_table)] <- ""

  # ensure all expected cols
  ann_cols <- if (!is.null(annotations)) names(annotations) else character(0)
  coef_cols <- unlist(lapply(seq_len(num_models), function(i) paste0("M",i,"_",stat)))
  expected <- c(ann_cols, "lhs","op","rhs", coef_cols)
  for (col in expected) if (!col %in% names(final_table)) final_table[[col]] <- ""
  final_table <- final_table[, expected, drop=FALSE]

  # set display names
  mapping <- if (standardized) c(est="β", se="SE", z="z", p="p") else c(est="b", se="SE", z="z", p="p")
  disp_stat <- sapply(stat, function(s) mapping[s] %||% s)
  header_labels <- c(ann_cols, "LHS","Operator","RHS", rep(disp_stat, times=num_models))

  # header grouping
  n_static <- 3 + length(ann_cols)
  header_groups <- c(" " = n_static)
  for (nm in model_names) header_groups[nm] <- length(stat)

  knitr::kable(final_table, caption=title, col.names=header_labels) %>%
    kableExtra::kable_classic(html_font="Cambria", full_width=full_width, position="left") %>%
    kableExtra::add_header_above(header_groups) %>%
    kableExtra::row_spec(nrow(merged_coeff), extra_css="border-bottom: 1px solid;") %>%
    kableExtra::footnote(general=notes, footnote_as_chunk=TRUE)
}
