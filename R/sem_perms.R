#' -----------------------------------------------------------------------------
#' @title      Run all permutations of a mediation/moderation SEM in Mplus
#' @description Automatically generates Mplus “.inp” files for every permutation
#'              of specified latent roles, runs them, and returns a combined
#'              data.frame of fit indices plus user‐requested parameter p‐values.
#'              Supports exclusion/fix rules, custom structural templates,
#'              covariates, and robust error‐handling with progress bars.
#'
#' @param data           A data.frame containing all manifest indicators and covariates.
#' @param latent_defs    A named character vector of measurement model strings.
#'                        E.g.
#'                        ```r
#'                        c(
#'                          x = "x by x1 x2 x3;",
#'                          m1 = "m1 by m11 m12 m13;",
#'                          …
#'                        )
#'                        ```
#' @param roles          Character vector of latent roles to permute (default `c("x","m1","m2","y")`).
#' @param covariate      Character vector of covariate names included in every model
#'                        (default `c("gender","age")`).
#' @param extra_var      Optional character vector of additional manifest variables.
#' @param exclusion      List of vectors: each vector lists roles that should not
#'                        co‐occur in the same model.
#' @param fix            Named list: each element is a vector of roles allowed
#'                        for one slot (e.g. `list(x = c("x","w"))`).
#' @param analysis       String passed to Mplus “ANALYSIS:” (default `"esti = ml;"`).
#' @param structural_tpl User‐supplied glue template for the “MODEL:” block.
#' @param output         String passed to Mplus “OUTPUT:” (default `"stand;"`).
#' @param out_dir        Directory in which to write `.inp` and `.out` files
#'                        (default `"./Mplus/med"`).
#' @param parameters     Character vector of parameter labels whose p‐values
#'                        you wish to extract (i.e., `"med"`).
#' @param minimal_run    Logical; if `TRUE`, only the first `max_run` permutations are generated and executed,
#'                        for a quick preview.
#' @param max_run        Integer; maximum number of models to generate and run when `minimal_run = TRUE`.
#' @param ...            Additional arguments (currently unused).
#'
#' @return A data.frame combining Mplus fit indices (via `fit_table()`) and
#'         p‐values for the specified parameters, with one row per model.
#'
#' @details
#' - Creates one `.inp` per permutation of `roles`.
#' - Honors `exclusion` and `fix` constraints.
#' - Writes data via `prepareMplusData()`, runs models with `runModels()` and
#'   a CLI progress bar (errors are caught and logged).
#' - Reads results with `readModels()`, builds fit‐index table, then appends
#'   requested p‐values.
#'
#' @examples
#' \dontrun{
#' # Example 1: chain mediation
#' data("latent_data")
#'
#' # Define measurement model for each varaible,
#' # must formatted without hyphen, "x1 x2 x3" rather than "x1-x3"
#' latent_defs <- c(
#'   x  = "x by x1 x2 x3;",
#'   m1 = "m1 by m11 m12 m13;",
#'   m2 = "m2 by m21 m22 m23 m24 m25;",
#'   y  = "y by y1 y2 y3 y4;",
#'   z  = "z by z1 z2 z3;",
#'   w  = "w by w1 w2 w3 w4;"
#' )
#'
#' #-- set up sample constraints for the examples ------------------------------
#’ # roles that should not co-occur (e.g., m1, m2 are from same scale)
#’ exclusive_groups <- list(m = c("m1", "m2"))
#’ # restrict the ‘x’ slot to either x or w
#’ fixed_groups     <- list(x = c("x", "w"))
#'
#' structural_tpl_chain <- "
#'   {m1} ON {x}
#'   {cov};
#'   {m2} ON {x}
#'   {m1}
#'   {cov};
#'   {y} ON {x}
#'   {m1}
#'   {m2}
#'   {cov};
#'
#'   model constraint:
#'   new(med1-med3);
#'   med1 = a1 * b1;
#'   med2 = a2 * b2;
#'   med3 = a1 * d * b2;"
#'
#' out1 <- sem_perms(
#'   data          = latent_data,
#'   latent_defs   = latent_defs,
#'   structural_tpl= structural_tpl_chain,
#'   covariate     = c("age", "gender"),
#'   out_dir       = "./Mplus",
#'   parameters = c("med1", "med2", "med3")
#' )
#'
#' # Example 2: moderation
#' structural_tpl_mod <- "
#'   int | {x} xwith {w};
#'   {y} ON {x} {w} int {cov};
#' "
#'
#' out2 <- sem_perms(
#'   data            = latent_data,
#'   latent_defs     = latent_defs,
#'   roles           = c("x", "w", "y"),
#'   structural_tpl  = structural_tpl_mod,
#'   covariate       = c("age", "gender"),
#'   parameters      = c("int"),
#'   analysis = "type = random;",
#'   out_dir         = "./Mplus/moderation"
#' )
#' }
#' @importFrom MplusAutomation prepareMplusData runModels readModels
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom glue glue glue_data
#' @export
sem_perms <- function(
    data,
    latent_defs,
    roles = NULL,          # Roles to permute
    covariate = NULL,             # Covariates for all models
    extra_var = NULL,
    exclusion = NULL,                   # Exclusive-group constraints
    fix = NULL,                             # Fixed-group constraints
    analysis = "esti = ml;",                      # Mplus analysis command
    structural_tpl,                                 # User-provided structural template
    output = "stand;",                            # Mplus output command
    out_dir = ".",
    parameters = NULL,
    minimal_run = TRUE,
    max_run = 10,
    ...
) {
  # ------------------------------
  # Helper functions
  # ------------------------------

  # 1. Helper: Extracts space-separated indicators from a "BY" statement
  parse_indicators <- function(def) {
    rhs <- sub(".*BY\\s*", "", def, ignore.case = TRUE)
    rhs <- sub("\\s*;\\s*$", "", rhs)
    vars <- strsplit(rhs, "\\s+")[[1]]
    paste(vars, collapse = " ")
  }

  # 2. Helper: run all .inp models with a progress bar and error handling
  run_with_progress <- function(target) {
    files <- list.files(target, "\\.inp$", full.names = TRUE)
    cli::cli_progress_bar("Running models", total = length(files))
    for (f in files) {
      cli::cli_progress_update()
      tryCatch(
        MplusAutomation::runModels(f),  # Robust call: logs failures but continues
        error = function(e) message(sprintf("Model '%s' failed: %s", f, e$message))
      )
    }
    cli::cli_progress_done()
  }

  # 3. Helper: extract parameters
  retrieveNewParams <- function(res, parameters) {
    parameters <- toupper(parameters)
    model_names <- names(res)
    result_list <- lapply(seq_along(res), function(i) {
      mod <- res[[i]]
      df_params <- mod$parameters$unstandardized
      subdf <- df_params[df_params$param %in% parameters, ]
      pvals <- subdf$pval[match(parameters, subdf$param)]

      out <- data.frame(Models = model_names[i], stringsAsFactors = FALSE)
      for (j in seq_along(parameters)) {
        col <- parameters[j]
        out[[col]] <- pvals[j]
      }
      out
    })
    do.call(rbind, result_list)
  }

  # ------------------------------
  # 1. Build permutation grid
  # ------------------------------
  pool <- names(latent_defs)
  n_roles <- length(roles)
  perm_mat <- gtools::permutations(n = length(pool), r = n_roles, v = pool)
  perm_df <- as.data.frame(perm_mat, stringsAsFactors = FALSE)
  names(perm_df) <- roles  # dynamic column names based on user-specified roles

  # ------------------------------
  # 2. Apply exclusion rules
  # ------------------------------
  if (!is.null(exclusion)) {
    for (vals in exclusion) {
      perm_df <- perm_df %>%
        dplyr::filter(rowSums(dplyr::across(all_of(roles), ~ . %in% vals)) <= 1)
      # Namespaced dplyr verbs to avoid conflicts
    }
  }

  # ------------------------------
  # 3. Apply fix rules
  # ------------------------------
  if (!is.null(fix)) {
    for (col in names(fix)) {
      allowed <- fix[[col]]
      perm_df <- perm_df %>%
        dplyr::filter(.data[[col]] %in% allowed)  # restrict specific roles
    }
  }

  # 4. Optionally reduce to a minimal run set
  if (isTRUE(minimal_run)) {
    n_keep <- min(nrow(perm_df), max_run)
    perm_df <- perm_df[seq_len(n_keep), , drop = FALSE]
  }

  # ------------------------------
  # 4. Prepare data and Mplus header
  # ------------------------------
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  header <- MplusAutomation::prepareMplusData(data, file.path(out_dir, "data.dat"))
  var_chunk <- paste0(header[2:5], collapse = "")  # assumes 4 lines in VARIABLE block

  # ------------------------------
  # 5. Loop over permutations to create .inp files
  # ------------------------------
  for (i in seq_len(nrow(perm_df))) {
    this_row <- as.list(perm_df[i, , drop = FALSE])
    indices <- unlist(this_row, use.names = FALSE)

    # Build USEVARIABLES list: indicators + covariates + extra variables
    inds_list <- vapply(latent_defs[indices], parse_indicators, FUN.VALUE = character(1))
    covariate = paste0(covariate, collapse = " ")
    usevariable <- paste(c(inds_list, covariate, extra_var), collapse = "\n")

    # measurement and Structural syntax
    measurement_syntax <- paste0(unname(latent_defs[indices]), collapse = "\n")
    structure_syntax <- glue::glue_data(this_row, structural_tpl, cov = covariate)

    # Assemble full .inp text
    syntax_string <- glue::glue(
      "title: mediation model;
            {var_chunk}
            usev = {usevariable};

            analysis:
            {analysis}

            model:
            !measurement
            {measurement_syntax}

            !structural
            {structure_syntax}

            output: {output}"
    )

    # Write the .inp file
    filename <- paste(indices, collapse = "_")
    create_mplus_inp(
      syntax_string = syntax_string,
      filename = filename,
      target_dir = out_dir,
      overwrite = TRUE
    )
  }

  # ------------------------------
  # 6. Run models and collect results
  # ------------------------------
  run_with_progress(out_dir)
  res <- MplusAutomation::readModels(out_dir, what = c("warn_err", "summaries", "parameters"))
  fit_tb <- fit_table(res, diffTest = FALSE)

  # ------------------------------
  # 7. Extract specified parameters
  # ------------------------------

  coeff_tb <- retrieveNewParams(res = res, parameters = parameters)
  final_df <- dplyr::left_join(fit_tb, coeff_tb, by = "Models")

  # Explicitly return the combined results
  return(final_df)
}

