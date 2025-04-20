library(bruceR)
library(MplusAutomation)

set.wd()

#  Read models
models <- readModels("./database")
m1 <- models$ex5.11.out
m2 <- models$ex5.13.out

coeff_table(m1, m2)
coeff_table(m1, m2,
            annotations = list(Type = c("good", rep("",2)),
                               Region = c(rep("", 2), "Bad")), model_names = cc("good, bad"))



coeff_table <- function(
    ...,
    ## basic appearance ---------------------------------------------------------
    stat          = c("est", "se", "p"),     # statistics to print
    indices       = character(0),            # extra fit indices
    digits_coeff  = 3,
    digits_fit    = 3,
    standardized  = FALSE,
    stars         = FALSE,
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

  ## 0 ▸ capture raw inputs ---------------------------------------------------
  raw_inputs <- list(...)

  ##############################################################################
  ##  .prep_args(): validate & normalise inputs                                ##
  ##############################################################################
  .prep_args <- function(
    dots,
    annotations   = NULL,
    stat          = c("est","se","p"),
    params        = c("regression"),
    indices       = character(0),
    order_by      = c("lhs","rhs"),
    standardized  = FALSE
  ) {

    ## 1. Expand / validate `stat` --------------------------------------------
    allowed_stat <- c("est","se","z","p")
    if (length(stat) == 1 && identical(stat, "all"))
      stat <- allowed_stat                      # user shortcut

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

    ## 4. Validate `annotations` ----------------------------------------------
    if (!is.null(annotations)) {
      ## 4‑a  shorthand: bare character vector -------------------------------
      if (is.atomic(annotations) && is.character(annotations)) {
        annotations <- list(Annotation = annotations)

        ## 4‑b data.frame: all char columns ------------------------------------
      } else if (is.data.frame(annotations)) {
        if (!all(vapply(annotations, is.character, logical(1))))
          stop("All columns in the `annotations` data.frame must be character.")

        ## 4‑c named list -------------------------------------------------------
      } else if (is.list(annotations)) {
        if (is.null(names(annotations)) || any(names(annotations) == ""))
          stop("`annotations` list must be *named*; names become column headers.")
        for (nm in names(annotations)) {
          x <- annotations[[nm]]
          ok_vec <- is.character(x)
          ok_map <- ok_vec && !is.null(names(x))  # regex map
          if (!(ok_vec || ok_map))
            stop("annotations[['", nm, "']] must be a character vector ",
                 "or a *named* character vector (regex map).")
        }

        ## 4‑d anything else ----------------------------------------------------
      } else {
        stop("`annotations` must be NULL, a character vector, ",
             "a character data.frame, or a named list of character vectors.")
      }
    }

    ## 5. Ensure at least one model object present ----------------------------
    is_model <- function(x) is.list(x) && ("summaries" %in% names(x))
    is_any_model <- function(obj) {
      if (is_model(obj)) TRUE
      else if (is.list(obj)) any(vapply(obj, is_any_model, logical(1)))
      else FALSE
    }
    if (!any(vapply(dots, is_any_model, logical(1))))
      stop("No MplusAutomation model objects detected in `...` ",
           "(each must contain a `$summaries` slot).")

    ## 6. Return cleaned list --------------------------------------------------
    list(
      dots         = dots,
      annotations  = annotations,
      stat         = stat,
      params       = params,
      indices      = indices,
      order_by     = order_by,
      standardized = standardized
    )
  }
  v <- .prep_args(raw_inputs, annotations,
                  stat, params, indices, order_by, standardized)
  ## ---- END of .prep_args() --------------------------------------------------
  ##############################################################################
  ##  Step 2 helper  –  flatten & label models                                 ##
  ##############################################################################
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
  model_tbl <- .flatten_models(raw_inputs, model_names)

  ##############################################################################
  ##  Step 3 helpers – parse & format                                         ##
  ##############################################################################

  .classify_param <- function(header, param) {
    if      (grepl("\\.BY$",   header)) "loading"
    else if (grepl("\\.ON$",   header)) "regression"
    else if (grepl("\\.WITH$", header)) "undirected"
    else if (header %in% c("Intercepts","Means"))  "expectation"
    else if (header %in% c("Variances","Residual.Variances")) "variability"
    else if (header == "New.Additional.Parameters")       "new"
    else                                                  "other"
  }

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

  ##############################################################################
  ##  Step 3 — extract one coefficient block                                  ##
  ##############################################################################
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
  }

  coeff_list <- Map(.coef_block, tib$model, tib$label,
                    MoreArgs = list(stat          = v$stat,
                                    params        = v$params,
                                    exclude       = exclude,
                                    standardized  = v$standardized,
                                    digits_coeff  = digits_coeff,
                                    stars         = stars))

  ##############################################################################
  ##  Step 4 : merge coefficient blocks and order rows                         ##
  ##############################################################################

  # Process function
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
  ##  Step 4 — Map over models to get list of fit blocks                      ##
  ##############################################################################
  .fit_block <- function(model_obj, model_label,
                         stat, indices, digits_fit) {

    ## 1. candidate indices ----------------------------------------------------
    default <- c("ChiSqM_Value","ChiSqM_DF","χ²/df",
                 "CFI","TLI","RMSEA_Estimate","SRMR")
    want    <- unique(c(default, indices))

    ## 2. convenient rename map -------------------------------------------------
    rename_map <- c(ChiSqM_Value   = "χ²",
                    ChiSqM_DF      = "df",
                    RMSEA_Estimate = "RMSEA")

    summ <- model_obj$summaries

    ## 3. pull / compute --------------------------------------------------------
    get_one <- function(idx) {
      if (idx == "χ²/df") {
        if (all(c("ChiSqM_Value","ChiSqM_DF") %in% names(summ)) &&
            summ$ChiSqM_DF != 0)
          return(summ$ChiSqM_Value / summ$ChiSqM_DF)
        return(NA_real_)
      }
      if (idx %in% names(summ)) return(summ[[idx]])
      NA_real_
    }
    raw_vals <- vapply(want, get_one, numeric(1))

    ## 4. rename & label --------------------------------------------------------
    labels <- ifelse(want %in% names(rename_map),
                     rename_map[want],
                     want)

    show   <- .fmt_num(raw_vals, digits_fit)
    ## keep df and Parameters integers un‑formatted
    show[ labels == "df" ] <- ifelse(is.na(raw_vals[labels == "df"]),
                                     NA_real_, as.integer(raw_vals[labels == "df"]))
    if ("Parameters" %in% want) {
      show[ labels == "Parameters" ] <- as.integer(raw_vals[labels == "Parameters"])
    }


    ## 5. build tidy block ------------------------------------------------------
    out <- data.frame(lhs = labels,
                      op  = "",
                      rhs = "",
                      stringsAsFactors = FALSE)

    for (s in stat) {
      out[[ paste0(model_label, "_", s) ]] <-
        if (s == stat[1]) show else ""        # only first stat column holds numbers
    }
  }

  fit_list   <- Map(.fit_block , tib$model, tib$label,
                    MoreArgs = list(stat       = v$stat,
                                    indices    = v$indices,
                                    digits_fit = digits_fit))

  ##############################################################################
  ## Step 5 — Merge model blocks, add annotations, produce master table            ##
  ##############################################################################
  .apply_annotations <- function(coeff_df, annotations) {
    if (is.null(annotations)) return(coeff_df)

    # Must be a named list of character vectors with matching length
    if (!is.list(annotations))
      stop("`annotations` must be a named list of character vectors.")

    for (col in names(annotations)) {
      vals <- annotations[[col]]

      if (!is.character(vals) || length(vals) != nrow(coeff_df)) {
        stop("Each element of `annotations[['", col, "']]` must be a character vector ",
             "with ", nrow(coeff_df), " rows.")
      }

      coeff_df[[col]] <- vals
    }
    coeff_df
  }

  ## 5‑A & 5‑B   consolidate model blocks --------------------------------------
  merged_coeff <- Reduce(function(x, y)
    dplyr::full_join(x, y, by = c("lhs","op","rhs")),
    coeff_list)

  merged_fit   <- Reduce(function(x, y)
    dplyr::full_join(x, y, by = c("lhs","op","rhs")),
    fit_list)

  ## 5‑C – apply user annotations ----------------------------------------------
  annotated_coeff <- .apply_annotations(merged_coeff, v$annotations)

  ## 5‑D – bind coeff and fit rows ---------------------------------------------
  final_df <- dplyr::bind_rows(annotated_coeff, merged_fit)
  final_df[is.na(final_df)] <- ""

  ## For debugging right now, RETURN the assembled data‑frame ------------------
  return(final_df)     # for now

}

