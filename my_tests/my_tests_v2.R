# Creating flow ----
# Creating the within function first, then using this function to retrieve key
# Information



# My test 2.R for coeff_table
library(bruceR)
library(MplusAutomation)

set.wd()
#  Read models
models <- readModels("./database")
## Arguments
m1 <- models$ex5.11.out
m2 <- models$ex5.13.out

# Examples for validation
coeff_table(m1)
coeff_table(m2)
coeff_table(m1, m2)
coeff_table(m1, m2,
            annotations = list(Type = c("good", rep("",2)),
                               Region = c(rep("", 2), "Bad")), model_names = cc("good, bad"))



# Main function ----

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
  
  ## 0a.  Expand `stat = "all"` into the full vector -------------------------
  if (length(stat)==1 && stat=="all") {
    stat <- c("est","se","z","p")
  }
  
  ## 0b.  From here on, just use `stat` everywhere, no further checks needed. --
  
  raw_inputs <- list(...)
  
  
  ##############################################################################
  ##  Step 1: prep_args(): Validate & normalise inputs                        ##
  ##############################################################################
  #  This helper is *not* exported; coeff_table() will call it internally.
  #  It returns a clean, named list of all the options coeff_table() needs later.
  #  ==> If any test fails, it stops with a clear, informative error message.
  
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
    
    ## 1. Expand / validate `stat` --------------------------------------------
    allowed_stat <- c("est","se","z","p")

    stat <- tolower(stat)
    if (!all(stat %in% allowed_stat))
      stop("`stat` must contain only ",
           paste(allowed_stat, collapse = ", "),
           ' or the single string "all".')
    
    ## 2. Validate `params` ----------------------------------------------------
    allowed_params <- c("loading","regression","expectation",
                        "undirected","variability","new","other")
    if (!all(params %in% allowed_params))
      stop("`params` must be a subset of ",
           paste(allowed_params, collapse = ", "), ".")
    
    ## 3. Normalise `order_by` -------------------------------------------------
    order_by <- match.arg(order_by)
    
    ## 4. Validate `annotations` as named list -------------------------------
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
    
    ## 6. Return cleaned list --------------------------------------------------
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
  ##  Step 2 helper  –  flatten & label models                                 ##
  ##############################################################################
  
  # Processing function
  .flatten_models <- function(model_inputs, override_names = NULL) {
    
    ## --- recursive collector ---------------------------------------------------
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
    
    ## Convert to tibble -------------------------------------------------------
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
  ##  Step 3 — extract one coefficient block                                  ##
  ##############################################################################
  
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
    
    ## choose parameter table --------------------------------------------------
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
                   tbl$rhs, ignore.case = TRUE)
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
  ##  Step 4 : merge coefficient blocks and order rows                         ##
  ##############################################################################
  
  ## ---- merge all coefficient blocks -----------------------------------------
  merged_coeff <- Reduce(function(a, b)
    dplyr::full_join(a, b, by = c("lhs","op","rhs")),
    coeff_list)
  
  ## ---- simple numeric precedence --------------------------------------------
  precedence <- match(merged_coeff$op,
                      c("=~", "~", "~ 1", "~~", ":="),         # desired order
                      nomatch = 99)                            # unknown ops last
  
  ## ---- final ordering --------------------------------------------------------
  order_side <- match.arg(order_by)      # "lhs" or "rhs"  (already validated)
  o <- order(precedence, merged_coeff[[order_side]], na.last = TRUE)
  
  ordered_coeff <- merged_coeff[o, , drop = FALSE]             # done
  
  
  ##############################################################################
  ## Step 5 — Add annotations, produce master table                           ##
  ##############################################################################
  
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
  ## Step 6 — Build fit index block                                           ##
  ##############################################################################
  
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
  ##  STEP 7 : bind coefficient‑ and fit‑blocks                                ##
  ##############################################################################
  
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
  
  ## 7‑C (For debugging) return, or pass to Step 8 (header & rendering) -------
  final_df
  
  ##############################################################################
  ## Step 8: Header & Spanner Construction for kableExtra                      ##
  ##############################################################################
  
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
  header_groups <- c(" " = stub_width)
  
  # One spanner per model, each covering length(stat) columns
  for (lbl in model_labels) {
    header_groups[lbl] <- length(stat)
  }
  
  # Paranoid repair: adjust stub if counts drift
  extra <- ncol(final_df) - sum(header_groups)
  if (extra != 0) {
    header_groups[1] <- header_groups[1] + extra
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
