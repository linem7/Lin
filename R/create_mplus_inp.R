#' Create Mplus .inp File
#'
#' Write a vector of Mplus syntax lines or a single string (with newlines) to an .inp file.
#'
#' @param syntax_string Character. A single string containing newlines or a character vector of lines defining the Mplus syntax.
#' @param filename Character. Desired name of the .inp file. Can include or omit the `.inp` extension.
#' @param target_dir Character. Path to the directory where the file should be written. Will be created if it doesn't exist.
#' @param overwrite Logical. Should an existing file be overwritten? Defaults to `FALSE`.
#' @param verbose Logical. Should a confirmation message be printed? Defaults to `TRUE`.
#' @return Invisibly returns the full path to the created .inp file.
#' @examples
#' # Single string example with embedded newlines
#' syntax <- "TITLE: CFA model;
#' DATA: FILE = 'data.dat';
#' VARIABLE: NAMES = x1-x5;
#' MODEL: f1 BY x1-x3;
#' OUTPUT: STANDARDIZED;"
#' # Write to ./Mplus/cfa.inp, overwriting if exists, and print message
#' create_mplus_inp(syntax_string = syntax, filename = "cfa", "./Mplus", overwrite = TRUE)
#' @export
create_mplus_inp <- function(syntax_string,
                             filename,
                             target_dir = ".",
                             overwrite = FALSE,
                             verbose = TRUE) {
  # Ensure target directory exists
  if (!dir.exists(target_dir)) {
    dir.create(target_dir, recursive = TRUE)
  }

  # Normalize filename and extension
  ext <- tools::file_ext(filename)
  if (tolower(ext) != "inp") {
    filename <- paste0(tools::file_path_sans_ext(filename), ".inp")
  }

  # Build full file path
  file_path <- file.path(target_dir, filename)

  # Prevent unintentional overwrite
  if (file.exists(file_path) && !overwrite) {
    stop("File already exists: ", file_path,
         ". To overwrite, set overwrite = TRUE.",
         call. = FALSE)
  }

  # Prepare lines for writing
  if (is.character(syntax_string) && length(syntax_string) == 1 && grepl("\n", syntax_string)) {
    lines <- strsplit(syntax_string, "\r?\n")[[1]]
  } else if (is.character(syntax_string)) {
    lines <- syntax_string
  } else {
    stop("`syntax_string` must be a character vector or single string containing newlines.")
  }

  # Write the file
  writeLines(lines, file_path)

  # Notify user
  if (verbose) {
    message("Mplus input file written to: ", normalizePath(file_path))
  }

  invisible(file_path)
}
