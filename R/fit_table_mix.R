#' Mixture Model Fit Indices Table and Plot
#'
#' Creates a concise summary table of fit indices for one or more latent-class
#' or latent-profile models read via \code{MplusAutomation}, and (optionally)
#' generates a ggplot2 line chart of these indices across class solutions.
#'
#' The reported indices include:
#' \describe{
#'   \item{Information criteria}{AIC, BIC, and sample-size adjusted BIC (aBIC)}
#'   \item{Entropy}{Model entropy}
#'   \item{Likelihood tests}{Vuong–Lo–Mendell–Rubin (LMR) and bootstrap LRT (BLRT) p-values}
#'   \item{Class-size proportions}{Percentages of most-likely class assignment}
#' }
#'
#' @param ... One or more model objects as returned by \code{MplusAutomation::readModels()}, each containing \code{summaries} and \code{class_counts}.
#' @param digits Integer ≥ 0. Number of decimal places for fit indices (\code{AIC}, \code{BIC}, \code{aBIC}, \code{Entropy}). Default: 3.
#' @param digits_prop Integer ≥ 0. Number of decimal places for class-size proportions (in the Proportion column). Default: 1.
#' @param plot Logical. If \code{TRUE}, returns a \code{list} with components \code{table} (a \code{data.frame} of fit indices) and \code{plot} (a \code{ggplot2} line chart). Default: \code{FALSE}.
#' @return A \code{data.frame} of fit indices, or (if \code{plot = TRUE}) a \code{list} with elements \code{table} and \code{plot}.
#' @details
#' The chart tries to parse the number of classes from the `Model` string via
#' regex. If your model titles don’t follow the “X class” pattern, the x-axis
#' may mislabel.
#'
#' @seealso \code{\link[MplusAutomation]{readModels}}
#'
#' @examples
#' \dontrun{
#'   res <- MplusAutomation::readModels("LCA_results")
#'   fit_table_mix(res)
#'   out <- fit_table_mix(res, digits = 2, plot = TRUE)
#'   out$table; out$plot
#' }
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate select across
#' @importFrom stringr str_replace
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_classic labs
#' @export
fit_table_mix <- function(..., digits = 3, digits_prop = 1, plot = FALSE) {

  ## ── 0. Flatten inputs ────────────────────────────────────────────────────
  raw_inputs <- list(...)
  mod_list   <- list()
  flatten <- function(x) {
    if (is.list(x) && !("summaries" %in% names(x))) {
      invisible(lapply(x, flatten))
    } else if (is.list(x) && ("summaries" %in% names(x))) {
      mod_list[[length(mod_list) + 1]] <<- x
    } else {
      stop("Each input must be a model (with 'summaries') or a list of them.",
           call. = FALSE)
    }
  }
  invisible(lapply(raw_inputs, flatten))
  n <- length(mod_list)
  if (n == 0) stop("No valid models supplied.", call. = FALSE)

  ## ── 1. Helper formatters (vectorised) ────────────────────────────────────
  .fmt_num <- function(v, d) vapply(v, function(x) {
    if (is.na(x) || x == "" || x %in% c("NA","<NA>")) "" else
      sprintf(paste0("%.", d, "f"), as.numeric(x))
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)

  .fmt_p <- function(v, d) vapply(v, function(x) {
    if (is.na(x) || x == "" || x %in% c("NA","<NA>")) "" else
      if (grepl("^<", x)) x else {
        val <- suppressWarnings(as.numeric(x))
        if (is.na(val)) "" else if (val < 0.001) "<0.001"
        else sprintf(paste0("%.", d, "f"), val)
      }
  }, FUN.VALUE = character(1), USE.NAMES = FALSE)

  ## ── 2. Placeholder -------------------------------------------------------
  cols <- c("Model","AIC","BIC","aBIC","Entropy","p_LMR","p_BLRT","Proportion (%)")
  holder <- data.frame(matrix(NA_character_, nrow = n, ncol = length(cols)),
                       stringsAsFactors = FALSE); names(holder) <- cols

  ## ── 3. Populate rows -----------------------------------------------------
  for (i in seq_len(n)) {
    smry <- mod_list[[i]][["summaries"]]
    holder$Model[i]   <- if ("Title"          %in% names(smry)) smry$Title[1]           else ""
    holder$AIC[i]     <- if ("AIC"            %in% names(smry)) smry$AIC[1]             else ""
    holder$BIC[i]     <- if ("BIC"            %in% names(smry)) smry$BIC[1]             else ""
    holder$aBIC[i]    <- if ("aBIC"           %in% names(smry)) smry$aBIC[1]            else ""
    holder$Entropy[i] <- if ("Entropy"        %in% names(smry)) smry$Entropy[1]         else ""
    holder$p_LMR[i]   <- if ("T11_LMR_PValue" %in% names(smry)) smry$T11_LMR_PValue[1]  else ""
    holder$p_BLRT[i]  <- if ("BLRT_PValue"    %in% names(smry)) smry$BLRT_PValue[1]     else ""

    prop <- mod_list[[i]][["class_counts"]][["mostLikely"]][["proportion"]]
    holder$`Proportion (%)`[i] <-
      if (is.null(prop)) "" else paste(
        sprintf(paste0("%.", digits_prop, "f"), sort(prop, TRUE) * 100), collapse="/")
  }

  ## ── 4. Format numbers ----------------------------------------------------
  holder[c("AIC","BIC","aBIC","Entropy")] <- lapply(
    holder[c("AIC","BIC","aBIC","Entropy")], .fmt_num, d = digits)
  holder[c("p_LMR","p_BLRT")] <- lapply(
    holder[c("p_LMR","p_BLRT")], .fmt_p, d = digits)

  rownames(holder) <- NULL

  ## ── 5. Optional ggplot ----------------------------------------------------
  if (!plot) return(holder)

  # prepare tidy data from formatted holder (numeric coercion on-the-fly)
  library(tidyr); library(dplyr); library(ggplot2); library(stringr)
  fits_long <- holder %>%
    mutate(Classes = str_replace(Model, " model.*", "")) %>%   # clean title
    select(Classes, AIC, BIC, aBIC) %>%
    mutate(across(-Classes, as.numeric)) %>%
    pivot_longer(AIC:aBIC, names_to = "Metric", values_to = "Value")

  p <- ggplot(fits_long,
              aes(x = Classes, y = Value,
                  group = Metric, lty = Metric, shape = Metric)) +
    geom_line() + geom_point() +
    theme_classic() +
    labs(x = "Number of Classes",
         title = "Model-fit indices by number of classes")

  list(table = holder, plot = p)
}
