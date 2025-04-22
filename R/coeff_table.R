#' Generate a Regression Coefficient Table with Fit Indices
#'
#' @description
#' `coeff_table()` takes one or more MplusAutomation model objects (or nested lists thereof),
#' extracts their regression (and other) parameter estimates plus fit‐index statistics,
#' and returns a publication‐ready table via **knitr**/**kableExtra**.  You can customize:
#'
#' - Which statistics to display (estimates, SE, z, p; “all” for every one).
#' - Whether to include standardized coefficients.
#' - Significance stars on estimates.
#' - Which parameter types to keep (regression paths, loadings, intercepts, etc.).
#' - Regex filters to drop unwanted rows.
#' - Row order by LHS or RHS of the formula.
#' - User‐supplied annotation columns (must be a named list of character vectors).
#' - Model header labels (overrides filenames).
#' - Caption and footnote text.
#'
#' @param ...            One or more MplusAutomation model objects (or lists thereof),
#'                       each containing a `$summaries` slot.
#' @param stat           Character vector of coefficient statistics to include.
#'                       Choices: `"est"`, `"se"`, `"z"`, `"p"`.  Use `"all"` to select all four.
#' @param indices        Character vector of additional fit‐index names to append
#'                       (e.g. `"AIC"`, `"BIC"`).
#' @param digits_coeff   Integer: number of decimal places for coefficient values.
#' @param digits_fit     Integer: number of decimal places for fit‐index values.
#' @param standardized   Logical: if `TRUE`, use standardized coefficients (β) instead of b.
#' @param stars          Logical: if `TRUE`, append significance stars based on p‐value.
#' @param params         Character vector of parameter classes to include.
#'                       Options: `"loading"`, `"regression"`, `"expectation"`,
#'                       `"undirected"`, `"variability"`, `"new"`, `"other"`.
#' @param exclude        Character vector of regular expressions; any row whose RHS
#'                       matches will be dropped.
#' @param order_by       Character: `"lhs"` or `"rhs"`, specifying the secondary sort key.
#' @param annotations    Named list of character vectors for manual row annotations.
#'                       Length must equal the number of retained coefficient rows.
#' @param model_names    Optional character vector to override the auto‐detected
#'                       model labels (filenames).  Must match the number of models.
#' @param title          Character: caption for the output table.
#' @param notes          Character: footnote text to display under the table.
#'
#' @return A **knitr**/**kableExtra** table object, with:
#'   - Annotation columns (if any) on the left
#'   - LHS | Operator | RHS structural columns
#'   - One block of coefficient rows, ordered by operator and `order_by`
#'   - One block of fit‐index rows appended below
#'   - Header spanners grouping each model’s statistic columns
#'
#' @examples
#' \dontrun{
#' # Single model, unstandardized regression paths
#' m1 <- readModels("path/to/model1.out")
#' coeff_table(m1, stat = c("est","se","p"))
#'
#' # Two models, all stats, custom AIC/BIC, and manual annotations
#' m2 <- readModels("path/to/model2.out")
#' coeff_table(m1, m2,
#'             stat        = "all",
#'             indices     = c("AIC","BIC"),
#'             annotations = list(
#'               Group = rep(c("A","B"), each = 5)
#'             ))
#' }
#'
#' @seealso \code{\link[knitr]{kable}}, \code{\link[kableExtra]{add_header_above}}
#'
#' @importFrom tibble tibble
#' @importFrom dplyr full_join bind_rows
#' @importFrom knitr kable
#' @importFrom kableExtra kable_classic add_header_above row_spec footnote
#' @export

coeff_table <- function(
    ...,
    ## basic appearance ---------------------------------------------------------
    stat          = c("est", "se", "p"),     # statistics to print
    indices       = character(0),            # extra fit indices
    digits_coeff  = 3,
    digits_fit    = 3,
    standardized  = FALSE,
    stars         = TRUE,
    ## rows / columns to keep ---------------------------------------------------
    params        = "regression",
    exclude       = character(0),
    order_by      = c("lhs", "rhs"),
    ## adornments ---------------------------------------------------------------
    annotations   = NULL,
    model_names   = NULL,
    title         = "Table X. Unstandardized regression coefficients and model‑fit indices.",
    notes         = "Unstandardized coefficients are reported."
) {

  # 0.  Handle the “all” shortcut and capture raw inputs ------------------------
  # If stat == "all", expand to every available statistic.
  # Then bundle ... into a list for downstream helpers.
  if (length(stat)==1 && stat=="all") {
    stat <- c("est","se","z","p")
  }
  raw_inputs <- list(...)


  ##############################################################################
  ## 1.  Validate & normalize inputs — internal helper `.prep_args()`         ##
  ##############################################################################

  # Processing funtion
  .prep_args <- function(
    dots = raw_inputs,
    annotations   = NULL,
    stat          = c("est","se","p"),
    params        = c("regression"),
    indices       = character(0),
    order_by      = c("lhs","rhs"),
    standardized  = FALSE
  ) {

    # Checks in turn:
    #  • stat: allowed values est, se, z, p
    allowed_stat <- c("est","se","z","p")

    stat <- tolower(stat)
    if (!all(stat %in% allowed_stat))
      stop("`stat` must contain only ",
           paste(allowed_stat, collapse = ", "),
           ' or the single string "all".')

    #  • params: must be one of the seven parameter classes
    allowed_params <- c("loading","regression","expectation",
                        "undirected","variability","new","other")
    if (!all(params %in% allowed_params))
      stop("`params` must be a subset of ",
           paste(allowed_params, collapse = ", "), ".")

    #  • order_by: match to "lhs" or "rhs"
    order_by <- match.arg(order_by)

    #  • annotations: if non‑NULL, a named list of character vectors
    if (!is.null(annotations)) {
      if (!is.list(annotations)) {
        stop("`annotations` must be a *named list* of character vectors.")
      }

      if (is.null(names(annotations)) || any(names(annotations) == "")) {
        stop("`annotations` list must be *named*; names become column headers.")
      }

      for (nm in names(annotations)) {
        vals <- annotations[[nm]]
        if (!is.character(vals)) {
          stop("annotations[['", nm, "']] must be a character vector.")
        }
      }
    }

    # Returns a list of cleaned options (stat, params, indices, etc.).
    list(
      dots = dots,
      annotations  = annotations,
      stat         = stat,
      params       = params,
      indices      = indices,
      order_by     = order_by,
      standardized = standardized
    )
  }

  # Processing results
  v <- .prep_args(raw_inputs, annotations,
                  stat, params, indices, order_by, standardized)

  ##############################################################################
  ## 2.  Flatten nested model inputs → tibble of (model, label)               ##
  ##############################################################################

  # Processing function
  .flatten_models <- function(model_inputs, override_names = NULL) {

    # Walks arbitrarily deep lists, finds each object with a `$summaries` slot,
    # pulls out its summaries$Filename as the label, and returns a tibble
    # with columns:
    #   • model — list‑column of model objects
    #   • label — character vector of unique model names
    collect <- function(obj) {
      out <- list()
      if (is.list(obj) && "summaries" %in% names(obj)) {     # a model object
        out[[ length(out) + 1 ]] <- list(
          model = obj,
          label = obj$summaries$Filename           # use filename directly
        )
      } else if (is.list(obj)) {                             # nested list
        for (elt in obj) out <- c(out, collect(elt))
      }
      out
    }
    flat <- collect(model_inputs)

    models_vec <- lapply(flat, `[[`, "model")
    name_vec   <- vapply(flat, `[[`, character(1), "label")

    tib <- tibble::tibble(
      model = models_vec,
      label = name_vec
    )

    ## --- user override ---------------------------------------------------------
    if (!is.null(override_names)) {
      if (length(override_names) != nrow(tib))
        stop("`model_names` length (", length(override_names),
             ") does not match the number of detected models (", nrow(tib), ").")
      tib$label <- override_names
    }

    ## --- guarantee unique header labels ---------------------------------------
    tib$label <- make.unique(tib$label)

    tib
  }

  # Processing results
  flatten_tib <- .flatten_models(raw_inputs, model_names)
  model_labels <- flatten_tib$label

  ##############################################################################
  ## 3.  Build per‑model coefficient block (`.coef_block()`)                  ##
  ##############################################################################
  # For one model, this:
  #  1. Chooses unstandardized or standardized table
  #  2. Parses each row into lhs | op | rhs
  #  3. Filters by params (operator type) and exclude (RHS regex)
  #  4. Formats numbers (digits, "<0.001", stars)
  #  5. Renames columns to ModelLabel_stat and returns a data.frame

  # Processing functions
  # Helper functions
  .create_formula_parts <- function(header, param) {
    if (grepl("\\.BY$", header)) {
      c(lhs = sub("\\.BY$", "", header), op = "=~", rhs = param)
    } else if (grepl("\\.ON$", header)) {
      c(lhs = sub("\\.ON$", "", header), op = "~", rhs = param)
    } else if (grepl("\\.WITH$", header)) {
      c(lhs = sub("\\.WITH$", "", header), op = "~~", rhs = param)
    } else if (header %in% c("Intercepts","Means")) {
      c(lhs = param, op = "~ 1", rhs = "")
    } else if (header %in% c("Variances","Residual.Variances")) {
      c(lhs = param, op = "~~", rhs = param)
    } else if (header == "New.Additional.Parameters") {
      c(lhs = "new", op = ":=", rhs = param)
    } else {
      c(lhs = header, op = "", rhs = param)
    }
  }

  .fmt_num <- function(x, digits = 3)
    ifelse(is.na(x), NA_character_, sprintf(paste0("%.", digits, "f"), x))

  .add_stars <- function(est, p) {
    pieces <- ifelse(is.na(p), "",
                     ifelse(p < .001, "***",
                            ifelse(p < .01, "**",
                                   ifelse(p < .05, "*", ""))))
    paste0(est, pieces)
  }

  # Main function to retrieve coefficients
  .coef_block <- function(model_obj, model_label,
                          stat, params, exclude,
                          standardized, digits_coeff, stars) {

    tbl <- if (standardized){
      model_obj$parameters$stdyx.standardized
    } else {
      model_obj$parameters$unstandardized
    }
    if (is.null(tbl)) {
      stop("Model ", model_label, ": requested ",
           if (standardized) "standardized" else "unstandardized",
           " parameters not found.")
    }

    tbl <- as.data.frame(tbl, stringsAsFactors = FALSE)

    ## unify p column ------------------------------------------------------
    names(tbl)[names(tbl) == "pval"] <- "p"


    ## build lhs/op/rhs --------------------------------------------------------
    parts <- t(mapply(.create_formula_parts, tbl$paramHeader, tbl$param))
    tbl   <- cbind(as.data.frame(parts, stringsAsFactors = FALSE), tbl)

    ## keep only rows whose operator matches `params` ---------------------------
    op_map <- c(loading     = "=~",
                regression   = "~",
                expectation  = "~ 1",
                undirected   = "~~",
                variability  = "~~",   # variance/covariance share '~~'
                new          = ":=")

    keep_ops <- unname(op_map[ params ])
    tbl      <- tbl[ tbl$op %in% keep_ops , ]

    ## drop by RHS pattern -----------------------------------------------------
    if (length(exclude)) {
      bad <- grepl(paste(exclude, collapse = "|"),
                   tbl$rhs,
                   ignore.case = TRUE,
                   perl = TRUE)
      tbl <- tbl[ !bad , ]
    }

    ## numeric formatting ------------------------------------------------------
    if ("est" %in% stat) tbl$est <- .fmt_num(tbl$est, digits_coeff)
    if ("se"  %in% stat) tbl$se  <- .fmt_num(tbl$se,  digits_coeff)
    if ("z"   %in% stat) tbl$z <- .fmt_num(tbl$est_se, digits_coeff)
    if ("p" %in% stat)
      tbl$p <- ifelse(!is.na(tbl$p) & tbl$p < .001,
                      "<0.001",
                      .fmt_num(tbl$p, digits_coeff))

    if (stars && "est" %in% stat)
      tbl$est <- .add_stars(tbl$est, tbl$p)

    ## keep only needed columns + rename per‑model -----------------------------
    keep <- c("lhs","op","rhs", stat)
    tbl  <- tbl[ , keep , drop = FALSE]
    names(tbl)[names(tbl) %in% stat] <- paste0(model_label, "_", stat)
    tbl
  }

  # Processing results
  coeff_list <- Map(.coef_block, flatten_tib$model, flatten_tib$label,
                    MoreArgs = list(stat          = v$stat,
                                    params        = v$params,
                                    exclude       = exclude,
                                    standardized  = v$standardized,
                                    digits_coeff  = digits_coeff,
                                    stars         = stars))


  ##############################################################################
  ## 4.  Merge & order coefficient rows                                       ##
  ##############################################################################
  #  • full_join all per‑model blocks by lhs/op/rhs
  #  • assign numeric precedence to each operator (=~, ~, ~ 1, ~~ , :=)
  #  • order by precedence then by user’s choice (lhs or rhs)

  ## ---- merge all coefficient blocks -----------------------------------------
  merged_coeff <- Reduce(function(a, b)
    dplyr::full_join(a, b, by = c("lhs","op","rhs")),
    coeff_list)

  ## ---- simple numeric precedence --------------------------------------------
  precedence <- match(merged_coeff$op,
                      c("=~", "~", "~ 1", "~~", ":="),         # desired order
                      nomatch = 99)                            # unknown ops last

  ## ---- final ordering --------------------------------------------------------
  o <- order(precedence, merged_coeff[[v$order_by]], na.last = TRUE) # using the arg from .prep_arg
  ordered_coeff <- merged_coeff[o, , drop = FALSE]             # done


  ##############################################################################
  ## 5.  Prepend user annotations                                              ##
  ##############################################################################
  # `.apply_annotations()` expects a named list of length‑n character vectors,
  # validates length, then moves annotation columns to the front.

  # Processing function
  .apply_annotations <- function(coeff_df, annotations) {

    if (is.null(annotations)) return(coeff_df)

    ## must be a *named* list ---------------------------------------------------
    if (!is.list(annotations) ||
        is.null(names(annotations)) || any(names(annotations) == ""))
      stop("`annotations` must be a *named list* of character vectors.")

    ## validate each column -----------------------------------------------------
    for (col in names(annotations)) {
      vals <- annotations[[col]]

      if (!is.character(vals) || length(vals) != nrow(coeff_df))
        stop("annotations[['", col, "']] must be a character vector with ",
             nrow(coeff_df), " elements.")

      coeff_df[[col]] <- vals     # add column
    }

    ## move annotation columns to the very front -------------------------------
    ann_cols <- names(annotations)
    coeff_df[ , c(ann_cols, setdiff(names(coeff_df), ann_cols)), drop = FALSE ]
  }

  # Processing results
  annotated_coeff <- .apply_annotations(ordered_coeff, annotations)


  ##############################################################################
  ## 6.  Build and merge fit‑index block (`.fit_block()`)                     ##
  ##############################################################################
  # For each model:
  #  1. Determine indices (core + user)
  #  2. Pull/compute χ²/df
  #  3. Create a skeleton with annotation + lhs/op/rhs
  #  4. Fill first stat column with values, others blank
  # Then full_join all models into `merged_fit`.

  # Processing function
  .fit_block <- function(model_obj, model_label,
                         stat, indices, digits_fit,
                         annotations = NULL) {

    ## work out which annotation columns exist (if any) ------------------------
    ann_cols <- if (is.null(annotations)) character(0) else names(annotations)

    ## 1. indices to pull ------------------------------------------------------
    core  <- c("ChiSqM_Value","ChiSqM_DF","χ²/df",
               "CFI","TLI","RMSEA_Estimate","SRMR")
    want  <- unique(c(core, indices))

    rename_map <- c(ChiSqM_Value = "χ²",
                    ChiSqM_DF    = "df",
                    RMSEA_Estimate = "RMSEA")

    s <- model_obj$summaries

    pull_one <- function(idx) {
      if (idx == "χ²/df" &&
          all(c("ChiSqM_Value","ChiSqM_DF") %in% names(s)) &&
          s$ChiSqM_DF != 0)
        return(s$ChiSqM_Value / s$ChiSqM_DF)
      if (idx %in% names(s)) return(s[[idx]])
      NA_real_
    }
    raw <- vapply(want, pull_one, numeric(1))

    lbl <- ifelse(want %in% names(rename_map), rename_map[want], want)

    show <- .fmt_num(raw, digits_fit)
    show[lbl=="df"] <- ifelse(is.na(raw[lbl=="df"]), "",
                              as.integer(raw[lbl=="df"]))
    if ("Parameters" %in% want)
      show[lbl=="Parameters"] <- as.integer(raw[lbl=="Parameters"])

    ## 2. build skeleton (ann_cols + lhs/op/rhs) -------------------------------
    skel <- data.frame(matrix("", nrow = length(lbl),
                              ncol = length(c(ann_cols,"lhs","op","rhs")),
                              dimnames = list(NULL, c(ann_cols,"lhs","op","rhs"))),
                       stringsAsFactors = FALSE)

    if (length(ann_cols))
      skel[[ann_cols[1]]] <- lbl else skel$lhs <- lbl

    skel$op  <- ""
    skel$rhs <- ""

    ## 3. put numbers only in the first stat column ----------------------------
    for (sname in stat)
      skel[[paste0(model_label,"_",sname)]] <- if (sname==stat[1]) show else ""

    skel
  }

  # Processing result
  ## Build list of per‑model fit blocks ----------------------------------------
  fit_list <- Map(.fit_block,
                  flatten_tib$model,
                  flatten_tib$label,
                  MoreArgs = list(stat        = stat,
                                  indices     = indices,
                                  digits_fit  = digits_fit,
                                  annotations = annotations))

  ## First merge all model fit tables ------------------------------------------
  merged_fit <- Reduce(function(a,b)
    dplyr::full_join(a, b,
                     by = c(names(annotations),
                            "lhs","op","rhs")),
    fit_list)


  ##############################################################################
  ## 7.  Bind coefficients + fit rows into `final_df`                         ##
  ##############################################################################
  #  • Ensure both pieces have identical columns (fill missing with "")
  #  • bind_rows(coeff, fit) and replace NA → ""

  ## 7‑A Ensure both frames have identical columns ----------------------------
  all_cols <- union(names(annotated_coeff), names(merged_fit))

  fill_missing <- function(df, cols) {
    missing <- setdiff(cols, names(df))
    if (length(missing))
      df[missing] <- ""
    df[ , cols, drop = FALSE ]            # identical column ordering
  }

  annotated_coeff <- fill_missing(annotated_coeff, all_cols)
  merged_fit      <- fill_missing(merged_fit,      all_cols)

  ## 7‑B Bind and clean --------------------------------------------------------
  final_df <- dplyr::bind_rows(annotated_coeff, merged_fit)
  final_df[is.na(final_df)] <- ""         # replace any lingering NA’s

  ## 7‑C  Drop row names so kable() doesn’t invent a first column -------------
  row.names(final_df) <- NULL

  ## 7‑D (For debugging) return, or pass to Step 8 (header & rendering) -------
  final_df

  ##############################################################################
  ## 8.  Header names & spanners for kableExtra                               ##
  ##############################################################################
  #  8A Map stat keys → b/β, SE, z, p
  #  8B Build col.names: annotations, LHS, Operator, RHS, then each model’s stats
  #  8C Compute header_groups:
  #     • stub: length(ann_cols) + 3
  #     • one spanner per model spanning length(stat)
  #  8D kable() + add_header_above() + row_spec() + footnote()

  # 8A. Map statistic keys to display labels (b vs β, SE, z, p) ----------------
  mapping <- if (standardized) {
    c(est = "β", se = "SE", z = "z", p = "p")
  } else {
    c(est = "b", se = "SE", z = "z", p = "p")
  }
  stat_labels    <- mapping[stat]        # e.g. c("b","SE","p")

  # 8B. Build full column‐name vector -----------------------------------------
  ann_cols       <- if (is.null(annotations)) character(0) else names(annotations)
  n_models       <- length(model_labels)
  model_colnames <- rep(stat_labels, times = n_models)

  new_colnames <- c(
    ann_cols,      # user’s annotation columns, if any
    "LHS", "Operator", "RHS",  # structural cols
    model_colnames             # each model’s stat columns
  )

  # 8C. Compute header spanner widths -----------------------------------------
  # Stub spans annotation + 3 structural columns
  stub_width    <- length(ann_cols) + 3
  header_groups <- c(" " = stub_width) # for add_header_above(" "=3, "Model A" = 3, "Model B" = 3)

  # One spanner per model, each covering length(stat) columns
  for (lbl in model_labels) {
    header_groups[lbl] <- length(stat)
  }

  # 8D. Render table with kableExtra -------------------------------------------
  knitr::kable(final_df,
               col.names = new_colnames,
               caption   = title) %>%
    kableExtra::kable_classic(full_width = FALSE, html_font = "Cambria", position="left") %>%
    kableExtra::add_header_above(header_groups) %>%
    kableExtra::row_spec(nrow(merged_coeff), extra_css="border-bottom: 1px solid;") %>%
    kableExtra::footnote(general = notes, footnote_as_chunk = TRUE)
}
