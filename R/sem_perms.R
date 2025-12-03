#' -----------------------------------------------------------------------------
#' @title       Run all permutations of a mediation/moderation SEM in Mplus
#' @description Automatically generates Mplus ".inp" files for every permutation
#'              of specified latent roles, runs them, and returns a combined
#'              data.frame. The output includes model structure (role assignments),
#'              standard fit indices (Chi-sq, df, CFI, TLI, RMSEA, SRMR), user-requested
#'              additional indices (e.g., AIC, BIC), and specific parameter p-values.
#'              Supports exclusion/fix rules, custom structural templates, covariates,
#'              automatic long-line wrapping for Mplus, and robust error-handling.
#'
#' @param data           A data.frame containing all manifest indicators and covariates.
#' @param latent_defs    A named character vector of measurement model strings.
#'                       The function robustly extracts variable names from these strings
#'                       by intersecting words with `names(data)`.
#'                       E.g.
#'                       ```r
#'                       c(
#'                         x = "x by x1 x2 x3;",
#'                         m1 = "m1 by m11 m12 m13;",
#'                         ...
#'                       )
#'                       ```
#' @param roles          Character vector of latent roles to permute (default `c("x","m","y")`).
#'                       These names will appear as column headers in the output data.frame
#'                       to indicate which variable was assigned to which role.
#' @param covariate      Character vector of covariate names included in every model
#'                       (default `NULL`).
#' @param extra_var      Optional character vector of additional manifest variables to
#'                       include in the `USEVARIABLES` list but not in the permutation logic.
#' @param exclusion      List of vectors: each vector lists variable names (from `names(latent_defs)`)
#'                       that should not co-occur in the same model.
#' @param fix            Named list: each element is a vector of variable names allowed
#'                       for a specific role slot (e.g. `list(x = c("var1","var2"))`).
#' @param analysis       String passed to Mplus "ANALYSIS:" block (default `"esti = ml;"`).
#' @param structural_tpl User-supplied glue template for the "MODEL:" block. Use `{role_name}`
#'                       placeholders (e.g., `{x}`, `{m}`) which will be replaced by actual
#'                       variable names.
#' @param output         String passed to Mplus "OUTPUT:" block (default `"stand;"`).
#' @param out_dir        Directory in which to write `.inp` and `.out` files
#'                       (default `"./"`).
#' @param parameters     Character vector of specific parameter labels (e.g., `"ind_eff"`)
#'                       whose p-values you wish to extract into the results table.
#' @param indices        Character vector of additional fit indices to extract from Mplus
#'                       summaries (e.g., `c("AIC", "BIC", "LL")`). Standard indices
#'                       (ChiSq, df, CFI, TLI, RMSEA, SRMR) are extracted by default.
#' @param minimal_run    Logical; if `TRUE`, only the first `max_run` permutations are
#'                       generated and executed, for a quick preview.
#' @param max_run        Integer; maximum number of models to generate and run when
#'                       `minimal_run = TRUE`.
#' @param ...            Additional arguments (currently unused).
#'
#' @return A data.frame with the following columns:
#' \itemize{
#'   \item \code{Models}: Unique model identifier (based on variable combination).
#'   \item \code{<roles>}: Columns named after `roles` showing the variable assigned to each role.
#'   \item \code{Fit Indices}: Standard indices (\eqn{\chi^2}, df, CFI, TLI, RMSEA, SRMR) plus any specified in `indices`.
#'   \item \code{Parameters}: P-values for parameters specified in `parameters`.
#'   \item \code{Warning/Error}: "Yes"/"No" indicators for Mplus execution status.
#' }
#' @details
#' \itemize{
#'   \item \strong{Workflow}: The function generates all valid combinations of variables based on
#'         \code{roles}, filters them using \code{exclusion} and \code{fix} constraints, and
#'         processes them in three steps: Generate Input -> Run Mplus -> Extract Results.
#'
#'   \item \strong{Template Syntax}: The \code{structural_tpl} uses glue syntax. Placeholders
#'         matching \code{roles} (e.g., \code{{x}}, \code{{m}}) are dynamically replaced by
#'         the specific variable names in each iteration.
#'
#'   \item \strong{Input Safety}: It robustly extracts variable names for \code{USEVARIABLES}
#'         by intersecting syntax words with your data columns (ignoring Mplus keywords like
#'         "BY" or "ON"). It also \strong{automatically wraps long lines} (< 85 chars) to
#'         prevent Mplus truncation errors.
#'
#'   \item \strong{Result Integrity}: Results are matched to models using the internal filename
#'         ID found within the Mplus output. This ensures strict data alignment and prevents
#'         \code{NA} rows caused by file sorting differences.
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Simple Mediation
#' data("latent_data")
#'
#' # 1. Define Measurement Models
#' latent_defs <- c(
#'   x  = "x_lat BY x1 x2 x3;",
#'   m  = "m_lat BY m1 m2 m3;",
#'   y  = "y_lat BY y1 y2 y3;",
#'   w  = "w_lat BY w1 w2 w3;"
#' )
#'
#' # 2. Define Constraints
#' # m and w cannot be in the same model
#' exclusive_groups <- list(c("m", "w"))
#' # y slot can only be y_lat
#' fixed_groups     <- list(y = c("y"))
#'
#' # 3. Define Structural Template
#' structural_tpl_med <- "
#'   {m} ON {x} (a);
#'   {y} ON {m} (b);
#'   {y} ON {x} (c);
#'   MODEL CONSTRAINT:
#'   NEW(ind);
#'   ind = a*b;
#' "
#'
#' # 4. Run Permutations
#' results <- sem_perms(
#'   data           = latent_data,
#'   latent_defs    = latent_defs,
#'   roles          = c("x", "m", "y"),
#'   structural_tpl = structural_tpl_med,
#'   exclusion      = exclusive_groups,
#'   fix            = fixed_groups,
#'   covariate      = c("age", "gender"),
#'   out_dir        = "./Mplus_Perms",
#'   parameters     = c("ind"),         # Extract mediation p-value
#'   indices        = c("AIC", "BIC")   # Extract extra fit indices
#' )
#'
#' # View results
#' head(results)
#' }
#' @importFrom MplusAutomation prepareMplusData runModels readModels
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom glue glue glue_data
#' @importFrom dplyr filter mutate select left_join bind_rows across
#' @importFrom tools file_path_sans_ext
#' @export
sem_perms <- function(
    data,
    latent_defs,
    roles = NULL,            # Roles to permute (e.g., c("x", "m", "y"))
    covariate = NULL,        # Covariates for all models
    extra_var = NULL,
    exclusion = NULL,        # Exclusive-group constraints
    fix = NULL,              # Fixed-group constraints
    analysis = "esti = ml;", # Mplus analysis command
    structural_tpl,          # User-provided structural template
    output = "stand;",       # Mplus output command
    out_dir = ".",
    parameters = NULL,       # Specific model parameters to extract
    indices = NULL,          # Additional fit indices
    minimal_run = TRUE,
    max_run = 10,
    ...
) {

  # --- 1. Helper: Run Models ---
  run_with_progress <- function(target) {
    files <- list.files(target, "\\.inp$", full.names = TRUE)
    if(length(files) == 0) return()
    cli::cli_progress_bar("Running models", total = length(files))
    for (f in files) {
      cli::cli_progress_update()
      tryCatch(MplusAutomation::runModels(f), error = function(e) NULL)
    }
    cli::cli_progress_done()
  }

  # --- 2. Helper: Clean Filename for Matching ---
  # Removes directory path and extension to get the raw ID (e.g., "ait_al_ap")
  get_clean_id <- function(filepath) {
    if(is.null(filepath) || is.na(filepath)) return(NA_character_)
    tools::file_path_sans_ext(basename(filepath))
  }

  # --- 3. Helper: Extract Info from a SINGLE model result ---
  extract_model_info <- function(model_res, target_indices, target_params) {
    # Initialize a single row dataframe
    out <- list()

    # A. Extract ID from Filename (CRITICAL STEP)
    # Use the filename recorded inside the summary to ensure matching
    if (!is.null(model_res$summaries$Filename)) {
      out$Models <- get_clean_id(model_res$summaries$Filename)
    } else {
      return(NULL) # Skip if no filename found
    }

    # B. Extract Fit Indices
    sm <- model_res$summaries
    if (is.data.frame(sm)) {
      defaults <- c("ChiSqM_Value", "ChiSqM_DF", "CFI", "TLI", "RMSEA_Estimate", "SRMR")
      cols_to_get <- unique(c(defaults, target_indices))

      for (col in cols_to_get) {
        if (col %in% names(sm)) {
          val <- sm[[col]]
          out[[col]] <- if(length(val)>0) as.numeric(val[1]) else NA_real_
        } else {
          out[[col]] <- NA_real_
        }
      }
    }

    # C. Extract Parameters
    if (!is.null(target_params)) {
      target_params <- toupper(target_params)
      # Initialize with NA
      for(p in target_params) out[[p]] <- NA_real_

      if (!is.null(model_res$parameters$unstandardized)) {
        dfp <- model_res$parameters$unstandardized
        if ("pval" %in% names(dfp)) {
          for (par in target_params) {
            # Find the parameter
            idx <- which(toupper(dfp$param) == par)
            if (length(idx) > 0) {
              out[[par]] <- dfp$pval[idx[1]]
            }
          }
        }
      }
    }

    # D. Warnings and Errors
    n_warn <- length(model_res$warnings)
    n_err  <- length(model_res$errors)
    out$Warning <- ifelse(n_warn > 0, "Yes", "No")
    out$Error   <- ifelse(n_err > 0, "Yes", "No")

    return(as.data.frame(out, stringsAsFactors = FALSE))
  }

  # ------------------------------
  # MAIN LOGIC
  # ------------------------------

  # 1. Prepare Data & Permutations
  pool <- names(latent_defs)
  n_roles <- length(roles)
  perm_mat <- gtools::permutations(n = length(pool), r = n_roles, v = pool)
  perm_df <- as.data.frame(perm_mat, stringsAsFactors = FALSE)
  names(perm_df) <- roles

  # 2. Exclusion Rules
  if (!is.null(exclusion)) {
    for (vals in exclusion) {
      perm_df <- perm_df %>%
        dplyr::filter(rowSums(dplyr::across(dplyr::all_of(roles), ~ . %in% vals)) <= 1)
    }
  }

  # 3. Fix Rules
  if (!is.null(fix)) {
    for (col in names(fix)) {
      if (col %in% names(perm_df)) {
        perm_df <- perm_df %>% dplyr::filter(.data[[col]] %in% fix[[col]])
      }
    }
  }

  # 4. Minimal Run
  if (isTRUE(minimal_run)) {
    n_keep <- min(nrow(perm_df), max_run)
    perm_df <- perm_df[seq_len(n_keep), , drop = FALSE]
  }

  # 5. Generate ID Column
  perm_df$Models <- apply(perm_df[, roles, drop=FALSE], 1, paste, collapse = "_")

  # 6. Prepare Mplus Files
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  header <- MplusAutomation::prepareMplusData(data, file.path(out_dir, "data.dat"))

  # Variable block extraction
  var_chunk <- paste0(grep("NAMES =", header, value = TRUE, ignore.case = TRUE), collapse="\n")
  if(var_chunk == "") var_chunk <- paste0(header[header %in% grep("NAMES =", header, value=T)], collapse=" ")
  if(var_chunk == "") var_chunk <- paste0("NAMES = ", paste(names(data), collapse=" "), ";")
  miss_line <- grep("MISSING", header, value=TRUE, ignore.case=TRUE)
  if(length(miss_line) > 0) var_chunk <- paste(var_chunk, miss_line[1], sep="\n")

  # 7. Write .inp Files
  write_mplus_input <- function(txt, fn, dir) writeLines(txt, file.path(dir, paste0(fn, ".inp")))

  for (i in seq_len(nrow(perm_df))) {
    this_row <- as.list(perm_df[i, roles, drop = FALSE])
    indices_vec <- unlist(this_row, use.names = FALSE)

    # Robust Variable Extraction
    current_latent_syntax <- paste(unname(latent_defs[indices_vec]), collapse = "\n")
    combined_text <- paste(current_latent_syntax, paste(covariate, collapse = " "), paste(extra_var, collapse = " "))
    all_tokens <- unlist(strsplit(gsub("[\n\t;]", " ", combined_text), "\\s+"))
    valid_items <- intersect(all_tokens, names(data))

    # Handle Line Wrapping
    usevariable <- paste(strwrap(paste(valid_items, collapse = " "), width = 85), collapse = "\n      ")

    cov_str = paste0(covariate, collapse = " ")
    structure_syntax <- glue::glue_data(this_row, structural_tpl, cov = cov_str)

    syntax_string <- glue::glue(
      "TITLE: Permutation run {i};
      DATA: FILE = data.dat;
      VARIABLE:
      {var_chunk}
      USEVARIABLES =
      {usevariable};
      ANALYSIS: {analysis}
      MODEL:
      {current_latent_syntax}
      {structure_syntax}
      OUTPUT: {output}"
    )

    filename <- perm_df$Models[i]
    if (exists("create_mplus_inp") && is.function(create_mplus_inp)) {
      create_mplus_inp(syntax_string, filename, out_dir, overwrite = TRUE)
    } else {
      write_mplus_input(syntax_string, filename, out_dir)
    }
  }

  # 8. Run Models
  run_with_progress(out_dir)

  # 9. Read Models
  # IMPORTANT: recursive=FALSE ensures we don't pick up old nested runs
  res_list <- MplusAutomation::readModels(out_dir, what = c("warn_err", "summaries", "parameters"), quiet = TRUE, recursive = FALSE)

  if (is.null(res_list) || length(res_list) == 0) {
    warning("No Mplus output files read. Please check the 'out_dir' for .out files.")
    return(perm_df) # Return just the structure if run failed
  }

  # 10. Extract Results into a DataFrame
  # Loop through list and build a data frame of results
  results_list <- lapply(res_list, extract_model_info, target_indices = indices, target_params = parameters)
  results_df <- dplyr::bind_rows(results_list)

  if (nrow(results_df) == 0) {
    warning("Outputs were read but no valid summaries found.")
    return(perm_df)
  }

  # 11. Merge Results with Permutation Map
  # Use inner_join or left_join based on 'Models' column
  final_df <- dplyr::left_join(perm_df, results_df, by = "Models")

  # 12. Rename Columns (Indices)
  rename_map <- c(
    "ChiSqM_Value" = "\u03c7\u00b2", # Chi-sq symbol
    "ChiSqM_DF"    = "df",
    "RMSEA_Estimate" = "RMSEA"
  )

  cols <- colnames(final_df)
  for (old in names(rename_map)) {
    if (old %in% cols) colnames(final_df)[cols == old] <- rename_map[old]
  }

  # 13. Final Column Reordering (Models -> Roles -> Indices -> Warn/Err)
  # Identify fit columns (all cols that are not Models, Roles, Warning, Error)
  all_cols <- names(final_df)
  special_cols <- c("Models", roles, "Warning", "Error")
  fit_cols <- setdiff(all_cols, special_cols)

  # Reorder: Models, Roles, Fit Indices, Warning, Error
  final_order <- c("Models", roles, fit_cols, "Warning", "Error")

  # Keep only columns that actually exist
  final_order <- intersect(final_order, names(final_df))

  return(final_df[, final_order, drop = FALSE])
}

