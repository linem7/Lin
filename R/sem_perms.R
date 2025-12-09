#' -----------------------------------------------------------------------------
#' @title        Run all permutations of an SEM (Latent, Observed, or Hybrid) in Mplus
#' @description  Automatically generates Mplus ".inp" files for every permutation
#'               of specified variable roles, runs them, and returns a combined
#'               data.frame.
#'
#'               This function supports:
#'               \itemize{
#'                 \item \strong{Latent Variable Models}: Standard SEM with measurement models.
#'                 \item \strong{Observed Variable Models}: Path analysis using raw data columns.
#'                 \item \strong{Hybrid Models}: Mix of latent and observed variables.
#'                 \item \strong{Complex Syntax}: RI-CLPM or models with constraints (e.g., `@1`, `@0`).
#'               }
#'
#'               The output includes model structure, standard fit indices, user-requested
#'               indices, and specific parameter p-values.
#'
#' @param data           A data.frame containing all manifest indicators and covariates.
#' @param latent_defs    A named list or character vector of variable definitions.
#'                       The function automatically detects the type of definition:
#'                       \itemize{
#'                         \item \strong{Observed Variable}: If the value matches a column name in `data`,
#'                               it is automatically formatted as `"var;"` for Mplus.
#'                         \item \strong{Latent/Complex}: If the value is Mplus syntax (e.g., "F1 BY x1-x3;"),
#'                               it is passed as-is.
#'                       }
#'                       Example:
#'                       ```r
#'                       list(
#'                         # Latent definition
#'                         "Stress" = "Stress BY s1 s2 s3;",
#'                         # Observed definition (Auto-detected)
#'                         "Age"    = "age_col"
#'                       )
#'                       ```
#' @param roles          Character vector of roles to permute (default `c("x","m","y")`).
#'                       These names act as placeholders in `structural_tpl` and column headers in the result.
#' @param covariate      Character vector of covariate names included in every model
#'                       (default `NULL`).
#' @param extra_var      Optional character vector of additional manifest variables to
#'                       include in the `USEVARIABLES` list but not in the permutation logic.
#' @param exclusion      List of vectors: each vector lists variable names (keys from `latent_defs`)
#'                       that should not co-occur in the same model.
#' @param fix            Named list: each element is a vector of variable names allowed
#'                       for a specific role slot (e.g. `list(x = c("Stress","Anxiety"))`).
#' @param analysis       String passed to Mplus "ANALYSIS:" block (default `"esti = ml;"`).
#' @param structural_tpl User-supplied glue template for the "MODEL:" block. Use `{role_name}`
#'                       placeholders (e.g., `{x}`, `{m}`) which will be replaced by the keys
#'                       from `latent_defs`.
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
#'   \item \strong{Robust Variable Extraction}: The function uses a regex-based extractor to identify
#'         required columns for `USEVARIABLES`. It handles syntax like `x1@1`, `x1-x5`, or `x*`
#'         correctly by cleaning symbols before matching against `names(data)`.
#'
#'   \item \strong{Template Syntax}: The \code{structural_tpl} uses glue syntax. Placeholders
#'         matching \code{roles} (e.g., \code{{x}}) are replaced by the \strong{names} (keys)
#'         of the `latent_defs` list.
#'
#'   \item \strong{Result Integrity}: Results are matched to models using the internal filename
#'         ID found within the Mplus output to ensure strict data alignment.
#' }
#'
#' @examples
#' \dontrun{
#'
#' # ==============================================================================
#' # Example 1: Latent Variable Mediation (Standard SEM)
#' # ==============================================================================
#' # 1. Define Measurement Models (Mplus Syntax)
#' latent_defs_sem <- list(
#'   "Stress"     = "Stress BY s1 s2 s3;",
#'   "Anxiety"    = "Anxiety BY a1 a2 a3;",
#'   "Depression" = "Depression BY d1 d2 d3;",
#'   "Burnout"    = "Burnout BY b1 b2 b3;"
#' )
#'
#' # 2. Define Structural Template (Using placeholders)
#' sem_structure <- "
#'   {y} ON {m} (b);
#'   {m} ON {x} (a);
#'   {y} ON {x} (cp);
#'   MODEL CONSTRAINT:
#'   NEW(ind);
#'   ind = a*b;
#' "
#'
#' # 3. Run Permutations
#' res_sem <- sem_perms(
#'   data = my_data,
#'   latent_defs = latent_defs_sem,
#'   roles = c("x", "m", "y"),
#'   structural_tpl = sem_structure,
#'   out_dir = "output_sem",
#'   parameters = c("ind") # Extract indirect effect p-value
#' )
#'
#' # ==============================================================================
#' # Example 2: Observed Variable Mediation (Path Analysis)
#' # ==============================================================================
#' # 1. Define Variables (Direct Column Names)
#' # Note: You do NOT need to add semicolons ";" here. The function adds them automatically.
#' obs_defs <- list(
#'   "SES"       = "ses_score",    # Colname in data is 'ses_score'
#'   "Education" = "edu_years",    # Colname in data is 'edu_years'
#'   "Income"    = "inc_log",      # Colname in data is 'inc_log'
#'   "Health"    = "health_idx"    # Colname in data is 'health_idx'
#' )
#'
#' # 2. Define Structural Template
#' # The structure logic is identical to SEM, but Mplus treats them as observed variables.
#' path_structure <- "
#'   {y} ON {m} (b);
#'   {m} ON {x} (a);
#'   {y} ON {x} (cp);
#'
#'   MODEL CONSTRAINT:
#'   NEW(ind total);
#'   ind = a*b;
#'   total = ind + cp;
#' "
#'
#' # 3. Run Permutations with Constraints
#' res_path <- sem_perms(
#'   data = my_data,
#'   latent_defs = obs_defs,
#'   roles = c("x", "m", "y"),
#'   structural_tpl = path_structure,
#'   # Example constraint: X must be SES or Education
#'   fix = list(x = c("SES", "Education")),
#'   out_dir = "output_path",
#'   parameters = c("ind", "total")
#' )
#' }
#' @importFrom MplusAutomation prepareMplusData runModels readModels
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom glue glue glue_data
#' @importFrom dplyr filter mutate select left_join bind_rows across
#' @importFrom tools file_path_sans_ext
#' @export
sem_perms <- function(
    data,
    latent_defs,             # List of definitions. Keys match {placeholders}, Values are Mplus syntax or var names.
    roles = NULL,            # Roles to permute (e.g., c("x", "m", "y"))
    covariate = NULL,        # Covariates for all models (vector of strings)
    extra_var = NULL,        # Any extra variables to include in USEVARIABLES
    exclusion = NULL,        # List of character vectors for exclusive-group constraints
    fix = NULL,              # List of character vectors for fixed-group constraints
    analysis = "esti = ml;", # Mplus ANALYSIS command
    structural_tpl,          # User-provided structural template (glue syntax)
    output = "stand;",       # Mplus OUTPUT command
    out_dir = ".",           # Directory for output files
    parameters = NULL,       # Specific model parameters to extract (e.g., "ind", "ab")
    indices = NULL,          # Additional fit indices to extract
    minimal_run = TRUE,      # If TRUE, runs only a subset (defined by max_run)
    max_run = 10,            # Maximum number of models to run if minimal_run is TRUE
    ...
) {

  # --- 1. Helper: Run Models ---
  # Iterates through .inp files in the target directory and runs them using MplusAutomation.
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
  # Removes directory path and extension to extract the raw Model ID.
  get_clean_id <- function(filepath) {
    if(is.null(filepath) || is.na(filepath)) return(NA_character_)
    tools::file_path_sans_ext(basename(filepath))
  }

  # --- 3. Helper: Extract Info from a SINGLE model result ---
  # Parses Mplus output to extract fit indices, specific parameters, and warnings/errors.
  extract_model_info <- function(model_res, target_indices, target_params) {
    out <- list()

    # Extract Model ID from filename to ensure correct merging later
    if (!is.null(model_res$summaries$Filename)) {
      out$Models <- get_clean_id(model_res$summaries$Filename)
    } else {
      return(NULL)
    }

    # Extract Fit Indices
    sm <- model_res$summaries
    if (is.data.frame(sm)) {
      defaults <- c("ChiSqM_Value", "ChiSqM_DF", "CFI", "TLI", "RMSEA_Estimate", "SRMR")
      cols_to_get <- unique(c(defaults, target_indices))
      for (col in cols_to_get) {
        # Check if column exists and extract the first value
        out[[col]] <- if(col %in% names(sm) && length(sm[[col]])>0) as.numeric(sm[[col]][1]) else NA_real_
      }
    }

    # Extract Specific Parameters (e.g., P-values)
    if (!is.null(target_params)) {
      target_params <- toupper(target_params)
      # Initialize with NA
      for(p in target_params) out[[p]] <- NA_real_

      if (!is.null(model_res$parameters$unstandardized)) {
        dfp <- model_res$parameters$unstandardized
        if ("pval" %in% names(dfp)) {
          for (par in target_params) {
            # Match parameter name (case-insensitive)
            idx <- which(toupper(dfp$param) == par)
            if (length(idx) > 0) out[[par]] <- dfp$pval[idx[1]]
          }
        }
      }
    }

    # Check for Warnings and Errors
    n_warn <- length(model_res$warnings)
    n_err  <- length(model_res$errors)
    out$Warning <- ifelse(n_warn > 0, "Yes", "No")
    out$Error   <- ifelse(n_err > 0, "Yes", "No")

    return(as.data.frame(out, stringsAsFactors = FALSE))
  }

  # --- [NEW] Helper: Robust Variable Extraction ---
  # Extracts pure variable names from complex Mplus syntax.
  # It cleans out symbols like @1, *, -, ; to identify which columns from the dataset
  # are actually required in USEVARIABLES.
  extract_vars_robust <- function(syntax_str, data_names) {
    if(is.null(syntax_str) || length(syntax_str) == 0) return(character(0))

    # 1. Replace non-variable characters (keep alphanumeric and underscore) with spaces
    #    Example: "x1@1" becomes "x1 1", "x1-x5" becomes "x1 x5"
    clean_str <- gsub("[^A-Za-z0-9_]+", " ", syntax_str)

    # 2. Split by whitespace and remove duplicates
    tokens <- unique(unlist(strsplit(clean_str, "\\s+")))

    # 3. Return only tokens that match actual column names in the dataset
    intersect(tokens, data_names)
  }

  # ------------------------------
  # MAIN LOGIC
  # ------------------------------

  # 1. Prepare Data & Permutations
  # Create a permutation matrix based on available definitions and roles
  pool <- names(latent_defs)
  n_roles <- length(roles)
  perm_mat <- gtools::permutations(n = length(pool), r = n_roles, v = pool)
  perm_df <- as.data.frame(perm_mat, stringsAsFactors = FALSE)
  names(perm_df) <- roles

  # 2. Exclusion Rules
  # Remove rows where multiple variables from the same exclusive group appear
  if (!is.null(exclusion)) {
    for (vals in exclusion) {
      perm_df <- perm_df %>%
        dplyr::filter(rowSums(dplyr::across(dplyr::all_of(roles), ~ . %in% vals)) <= 1)
    }
  }

  # 3. Fix Rules
  # Restrict specific roles to specific subsets of variables
  if (!is.null(fix)) {
    for (col in names(fix)) {
      if (col %in% names(perm_df)) {
        perm_df <- perm_df %>% dplyr::filter(.data[[col]] %in% fix[[col]])
      }
    }
  }

  # 4. Minimal Run (For testing purposes)
  if (isTRUE(minimal_run)) {
    n_keep <- min(nrow(perm_df), max_run)
    perm_df <- perm_df[seq_len(n_keep), , drop = FALSE]
  }

  # 5. Generate Unique Model ID
  perm_df$Models <- apply(perm_df[, roles, drop=FALSE], 1, paste, collapse = "_")

  # 6. Prepare Mplus Files
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Convert R data to Mplus .dat format
  header <- MplusAutomation::prepareMplusData(data, file.path(out_dir, "data.dat"))

  # Extract variable names header from prepareMplusData output
  # This ensures the input file matches the .dat file structure
  var_chunk <- paste0(grep("NAMES =", header, value = TRUE, ignore.case = TRUE), collapse="\n")
  if(var_chunk == "") var_chunk <- paste0(header[header %in% grep("NAMES =", header, value=T)], collapse=" ")
  if(var_chunk == "") var_chunk <- paste0("NAMES = ", paste(names(data), collapse=" "), ";")

  # Check for missing value declaration
  miss_line <- grep("MISSING", header, value=TRUE, ignore.case=TRUE)
  if(length(miss_line) > 0) var_chunk <- paste(var_chunk, miss_line[1], sep="\n")

  # 7. Write .inp Files (Loop through permutations)
  write_mplus_input <- function(txt, fn, dir) writeLines(txt, file.path(dir, paste0(fn, ".inp")))

  for (i in seq_len(nrow(perm_df))) {
    this_row <- as.list(perm_df[i, roles, drop = FALSE])
    indices_vec <- unlist(this_row, use.names = FALSE)

    # --- [MODIFIED LOGIC START] ---

    # A. Smart Handling of Latent Definitions
    # Iterates through selected definitions.
    # - If a definition is just a single variable name found in the data -> Treat as Observed Variable.
    #   (Automatically adds ";" to request variance estimation).
    # - If a definition is complex (contains "BY", "@", etc.) -> Treat as Mplus Syntax.
    current_defs_list <- lapply(latent_defs[indices_vec], function(def) {
      if (length(def) == 1 && def %in% names(data)) {
        return(paste0(def, ";")) # Transform "var" to "var;"
      } else {
        return(def) # Keep complex syntax as is
      }
    })
    current_latent_syntax <- paste(unlist(current_defs_list), collapse = "\n")

    # B. Robust Variable Extraction (Generating USEVARIABLES)
    # Combine all syntax parts that might contain variable names
    combined_text <- paste(
      current_latent_syntax,
      paste(covariate, collapse = " "),
      paste(extra_var, collapse = " ")
    )

    # Extract only valid column names from the text to populate USEVARIABLES
    valid_items <- extract_vars_robust(combined_text, names(data))

    # --- [MODIFIED LOGIC END] ---

    # Format USEVARIABLES string with line wrapping (width=85)
    usevariable <- paste(strwrap(paste(valid_items, collapse = " "), width = 85), collapse = "\n      ")

    # Prepare Covariate string
    cov_str = paste0(covariate, collapse = " ")

    # Fill placeholders in the structural template using current permutation roles
    structure_syntax <- glue::glue_data(this_row, structural_tpl, cov = cov_str)

    # Assemble the final Mplus Input String
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

    # Write to file (Compatibility check for external writer functions)
    if (exists("create_mplus_inp") && is.function(create_mplus_inp)) {
      create_mplus_inp(syntax_string, filename, out_dir, overwrite = TRUE)
    } else {
      write_mplus_input(syntax_string, filename, out_dir)
    }
  }

  # 8. Run Models
  run_with_progress(out_dir)

  # 9. Read Models
  # Recursive=FALSE prevents reading subdirectories or old runs
  res_list <- MplusAutomation::readModels(out_dir, what = c("warn_err", "summaries", "parameters"), quiet = TRUE, recursive = FALSE)

  if (is.null(res_list) || length(res_list) == 0) {
    warning("No Mplus output files read. Please check the directory.")
    return(perm_df)
  }

  # 10. Extract Results into DataFrame
  results_list <- lapply(res_list, extract_model_info, target_indices = indices, target_params = parameters)
  results_df <- dplyr::bind_rows(results_list)

  if (nrow(results_df) == 0) {
    warning("Outputs were read but no valid summaries found.")
    return(perm_df)
  }

  # 11. Merge Results with Permutation Map
  final_df <- dplyr::left_join(perm_df, results_df, by = "Models")

  # 12. Rename Standard Columns
  rename_map <- c(
    "ChiSqM_Value"   = "\u03c7\u00b2", # Chi-square symbol
    "ChiSqM_DF"      = "df",
    "RMSEA_Estimate" = "RMSEA"
  )
  cols <- colnames(final_df)
  for (old in names(rename_map)) {
    if (old %in% cols) colnames(final_df)[cols == old] <- rename_map[old]
  }

  # 13. Final Column Reordering
  # Order: Models -> Roles -> Fit Indices -> Warnings/Errors
  all_cols <- names(final_df)
  special_cols <- c("Models", roles, "Warning", "Error")
  fit_cols <- setdiff(all_cols, special_cols)
  final_order <- intersect(c("Models", roles, fit_cols, "Warning", "Error"), names(final_df))

  return(final_df[, final_order, drop = FALSE])
}

