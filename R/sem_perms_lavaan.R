#' sem_perms_lavaan: Run permutation-based SEM/CFA using lavaan
#'
#' This function mirrors the functionality of a permutation SEM pipeline previously built
#' with MplusAutomation, but implements everything in R via lavaan. It systematically
#' permutes the assignment of measurement models (latent factors) to roles and extracts
#' fit indices and specified parameter p-values for each permutation.
#'
#' @param data A data.frame containing the observed manifest variables.
#' @param latent_defs A named list of lavaan measurement strings, e.g.:
#'   list(
#'     x  = "x =~ x1 + x2 + x3",
#'     m1 = "m1 =~ m11 + m12 + m13",
#'     m2 = "m2 =~ m21 + m22 + m23 + m24 + m25",
#'     y  = "y =~ y1 + y2 + y3 + y4"
#'   )
#'   Each name corresponds to a factor label; these labels are permuted across `roles`.
#' @param roles Character vector of factor names to permute, e.g. `c("x","m1","m2","y")`.
#' @param covariate Character vector of observed covariates to include in every model.
#' @param extra_var Character vector of additional observed variables to include.
#' @param exclusion List of character vectors; factors in each inner vector will not co-occur.
#' @param fix Named list: names are positions in `roles`, values are allowed factor names.
#' @param structural_tpl A glue template for the structural part. Use single braces for placeholders matching `roles`, plus `{cov}` and `{extra}`.
#'   e.g.: `{y} ~ a*{x} + b*{m1} + d*{m2}{cov}{extra}\nind := a * b`
#' @param minimal_run Logical; if `TRUE`, only the first `max_run` permutations are run (quick check).
#' @param max_run Integer; number of permutations to run when `minimal_run = TRUE`.
#' @param parameters Character vector of parameter labels (e.g., `c("a","ind")`) for p-value extraction.
#' @param missing lavaan missing-data method (default: "fiml").
#' @param estimator lavaan estimator (default: "ML").
#' @param digits Decimal places for fit indices & p-values (default: 3).
#' @return A data.frame with one row per permutation:
#'   * `Model`: permutation label (e.g. "m1_m2_x_y")
#'   * `χ²(df)`, `χ²/df`, `CFI`, `TLI`, `RMSEA`, `SRMR`, `AIC`, `BIC`
#'   * one column per entry in `parameters` with p-values
#'   * `Error`: factor (Yes/No) if lavaan parsing/fitting failed
#'   * `Warning`: factor (Yes/No) if lavaan produced a warning during fitting
#'
#' @examples
#' # Example 1: common list format (no helper function)
#' latent_defs <- list(
#'   x  = "x =~ x1 + x2 + x3",
#'   m1 = "m1 =~ m11 + m12 + m13",
#'   m2 = "m2 =~ m21 + m22 + m23 + m24 + m25",
#'   y  = "y =~ y1 + y2 + y3 + y4"
#' )
#' result_list <- sem_perms_lavaan(
#'   data           = latent_data,
#'   latent_defs    = latent_defs,
#'   roles          = c("x","m1","m2","y"),
#'   structural_tpl = "{y} ~ a*{x} + b*{m1} + d*{m2}{cov}{extra}",
#'   minimal_run    = TRUE,
#'   max_run        = 3,
#'   parameters     = c("a","b"),
#'   covariate      = c("age"),
#'   extra_var      = c("gender")
#' )
#'
#' # Example 2: using create_cfa_model helper
#' latent_defs2 <- list(
#'   x  = create_cfa_model("x",  "x{i}", 1:3),
#'   m1 = create_cfa_model("m1", "m1{i}", 1:3),
#'   m2 = create_cfa_model("m2", "m2{i}", 1:5),
#'   y  = create_cfa_model("y",  "y{i}", 1:4)
#' )
#' result_helper <- sem_perms_lavaan(
#'   data           = latent_data,
#'   latent_defs    = latent_defs2,
#'   roles          = c("x","m1","m2","y"),
#'   structural_tpl = "{y} ~ a*{x} + b*{m1} + d*{m2}{cov}{extra}
#'   ind := a * b",
#'   minimal_run    = TRUE,
#'   max_run        = 3,
#'   parameters     = c("a","ind"),
#'   covariate      = c("age","gender"),
#'   extra_var      = c("z1","z2","z3")
#' )
#' @export
sem_perms_lavaan <- function(
    data,
    latent_defs,
    roles          = c("x", "m1", "m2", "y"),
    covariate      = NULL,
    extra_var      = NULL,
    exclusion      = NULL,
    fix            = NULL,
    structural_tpl,
    minimal_run    = FALSE,
    max_run        = 10,
    parameters     = NULL,
    missing        = "fiml",
    estimator      = "ML",
    digits         = 3
) {
  requireNamespace("gtools")
  requireNamespace("lavaan")
  requireNamespace("glue")
  requireNamespace("dplyr")
  requireNamespace("cli")

  pool     <- names(latent_defs)
  perm_mat <- gtools::permutations(n = length(pool), r = length(roles), v = pool)
  perm_df  <- as.data.frame(perm_mat, stringsAsFactors = FALSE)
  names(perm_df) <- roles
  if (!is.null(exclusion)) for (v in exclusion) perm_df <- perm_df[rowSums(perm_df %in% v)<=1, ]
  if (!is.null(fix)) for (pos in names(fix)) perm_df <- perm_df[perm_df[[pos]] %in% fix[[pos]], ]
  if (minimal_run) perm_df <- head(perm_df, max_run)

  results_list <- vector("list", nrow(perm_df))
  cli::cli_progress_bar("Fitting lavaan models", total=nrow(perm_df))
  fit_indices <- c("chisq","df","cfi","tli","rmsea","srmr","aic","bic")

  for (i in seq_len(nrow(perm_df))) {
    cli::cli_progress_update()
    row_i <- perm_df[i, , drop = FALSE]
    idxs  <- unlist(row_i)

    # Measurement part
    meas <- paste(latent_defs[idxs], collapse = "\n")
    cov_s <- if (!is.null(covariate)) paste0(" + ", paste(covariate, collapse = " + ")) else ""
    ext_s <- if (!is.null(extra_var)) paste0(" + ", paste(extra_var, collapse = " + ")) else ""

    # Structural part
    struct <- glue::glue_data(row_i, structural_tpl, cov = cov_s, extra = ext_s)
    mstr   <- paste(meas, struct, sep = "\n\n")

    # --- MODIFICATION START ---
    # Initialize error and warning flags for the current iteration
    ef <- FALSE # Error Flag
    wf <- FALSE # Warning Flag

    # Fit and flag errors/warnings
    fit <- tryCatch({
      # Use withCallingHandlers to catch warnings without stopping execution.
      # The main expression (lavaan::sem) is evaluated. If a warning occurs,
      # the handler is called (setting wf to TRUE), and then execution resumes.
      # The return value of the main expression is preserved.
      withCallingHandlers(
        lavaan::sem(mstr, data = data, missing = missing, estimator = estimator),
        warning = function(w) {
          wf <<- TRUE
          # We can "muffle" the warning so it doesn't print to the console
          # during the loop by uncommenting the next line:
          # invokeRestart("muffleWarning")
        }
      )
    },
    # The outer tryCatch still handles terminal errors as before.
    error = function(e) {
      ef <<- TRUE
      NULL # Return NULL for the fit object on error
    })
    # --- MODIFICATION END ---

    # Extract fit measures
    fm <- tryCatch(
      lavaan::fitMeasures(fit, fit_indices),
      error = function(e) setNames(rep(NA_real_, length(fit_indices)), fit_indices) # Handle missing values
    )

    # Calculate the ratio χ²/df
    ratio_chi2_df <- if (!is.na(fm["chisq"]) && !is.na(fm["df"]) && fm["df"] > 0) {
      fm["chisq"] / fm["df"]
    } else {
      NA
    }

    # Format the fit indices
    fm_formatted <- formatC(fm, digits = 3, format = "f")
    fm_formatted["df"] <- as.character(as.integer(round(fm["df"], 0))) # Make df an integer

    # Calculate and add χ²/df to the formatted results
    fm_formatted["χ²/df"] <- formatC(ratio_chi2_df, digits = 3, format = "f")

    # Create the fit container and rename appropriately
    df_row <- data.frame(
      Model     = paste(idxs, collapse = "_"),
      `χ²(df)`  = sprintf("%s (%s)", fm_formatted["chisq"], fm_formatted["df"]),
      `χ²/df`   = fm_formatted["χ²/df"],
      CFI       = fm_formatted["cfi"],
      TLI       = fm_formatted["tli"],
      RMSEA     = fm_formatted["rmsea"],
      SRMR      = fm_formatted["srmr"],
      AIC       = fm_formatted["aic"],
      BIC       = fm_formatted["bic"],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    # Add p-values if applicable
    if (!ef && !is.null(parameters)) {
      pv <- setNames(rep(NA_character_, length(parameters)), parameters)
      # Check if fit object exists before trying to extract parameters
      if (!is.null(fit)) {
        pe <- tryCatch(lavaan::parameterEstimates(fit, ci = FALSE), error = function(e) NULL)
        if (!is.null(pe)) {
          for (par in parameters) {
            ii <- which(pe$label == par | paste(pe$lhs, pe$op, pe$rhs, sep = "") == par)
            if (length(ii)) pv[par] <- formatC(pe$pvalue[ii[1]], digits = 3, format = "f")
          }
        }
      }
      for (par in parameters) df_row[[par]] <- pv[par]
    }

    # Add error and warning flags
    df_row$Error   <- factor(ifelse(ef, "Yes", "No"), levels = c("Yes", "No"))
    df_row$Warning <- factor(ifelse(wf, "Yes", "No"), levels = c("Yes", "No"))

    results_list[[i]] <- df_row
  }

  cli::cli_progress_done()
  final_df <- dplyr::bind_rows(results_list)
  row.names(final_df) <- NULL
  final_df
}
