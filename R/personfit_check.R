#' Person-fit diagnostic plots
#'
#' @description
#' Creates histogram plots of person-fit statistics with empirical cutoffs using
#' the PerFit package. Plots parametric lz*, nonparametric Ht, or both side-by-side.
#'
#' @details
#' Creates diagnostic histogram plots for person-fit analysis to identify aberrant
#' response patterns. Person-fit statistics quantify how well an individual's
#' response pattern matches expectations from an Item Response Theory (IRT) model
#' or nonparametric assumptions.
#'
#' \strong{Parametric lz* statistic}: This method is based on the standardized log-likelihood.
#' It requires fitting an IRT model (Rasch, 2PL, or 3PL) to obtain two key components:
#' \emph{item parameters} (e.g., difficulty, discrimination) and \emph{person parameters}
#' (i.e., ability estimates, or theta scores). Both of these are passed to the lz*
#' function to compute a fit score for each person. More extreme negative values
#' indicate a worse fit between the person's responses and the model's predictions.
#'
#' \strong{Nonparametric Ht statistic}: This method evaluates response consistency based
#' on Guttman patterns without an IRT model. It checks if examinees answer easier
#' items correctly and harder items incorrectly, flagging unexpected deviations.
#'
#' Both statistics are visualized as histograms. A critical cutoff value, determined
#' via bootstrap resampling, is drawn on the plot. Individuals beyond this cutoff
#' are flagged as potentially aberrant, aiding in data quality control.
#' The function computes person-fit statistics and their cutoffs, then creates
#' diagnostic histograms. For parametric analysis, it fits an IRT model internally
#' using mirt with IRTpars=TRUE. The returned plot object can be saved or replayed.
#'
#'
#' @param data Matrix or data.frame (N x J) of dichotomous responses (0/1, NA allowed)
#' @param method Character. One of "parametric", "nonparametric", or "both" (default). Accepted aliases, "par"/"para" -> "parametric", "non"/"np" -> "nonparametric", "b"/"all" -> "both"
#' @param model Character. IRT model: "Rasch", "2PL" (default), or "3PL". Used only for parametric. Accepted aliases, "rasch"/"1pl" -> "Rasch", "2pl" -> "2PL", "3pl" -> "3PL"
#' @param cutpoint Numeric. Significance level for cutoff: 0.05 (default) or 0.10
#' @param score_method Character. Ability estimation: "EAP" (default), "MAP", or "ML"
#' @param seed Integer. Random seed for cutoff computation. Default 42
#'
#' @return A recorded plot object that can be replayed with `replayPlot()`
#'
#' @examples
#' # These exmaples take times to generate, please wait
#'
#' # Load example data
#' data(InadequacyData, package = "PerFit")
#'
#' # Example 1: Nonparametric only (no mirt needed)
#' ht_plot <- personfit_check(
#'   data = InadequacyData,
#'   method = "nonparametric",
#'   cutpoint = 0.05
#' )
#'
#' \donttest{
#' # Example 2: Parametric with 2PL model (requires mirt)
#' lz_plot <- personfit_check(
#'   data = InadequacyData,
#'   method = "parametric",
#'   model = "2PL",
#'   score_method = "MAP"
#' )
#'
#' # Example 3: Both methods side-by-side
#' both_plot <- personfit_check(
#'   data = InadequacyData,
#'   method = "both",
#'   model = "3PL",
#'   cutpoint = 0.10
#' )
#'
#' # Example 4: Using aliases
#' quick_plot <- personfit_check(
#'   data = InadequacyData,
#'   method = "par",    # -> "parametric"
#'   model = "rasch"    # -> "Rasch"
#' )
#'
#' # Example 5: Save plot to file
#' my_plot <- personfit_check(InadequacyData, method = "both")
#'
#' # Save as PDF
#' pdf(tempfile(fileext = ".pdf"), width = 12, height = 6)
#' replayPlot(my_plot)
#' dev.off()
#'
#' # Save as PNG
#' png(tempfile(fileext = ".png"), width = 1200, height = 600)
#' replayPlot(my_plot)
#' dev.off()
#' }
#'
#' @seealso
#' \code{\link[PerFit]{lzstar}}, \code{\link[PerFit]{Ht}},
#' \code{\link[PerFit]{cutoff}}, \code{\link[mirt]{mirt}}
#'
#' @export
personfit_check <- function(
    data,                   # N x J matrix/data.frame of 0/1 (can include NA)
    method = "both",        # "parametric", "nonparametric", or "both"
    model  = "2PL",         # "Rasch", "2PL", or "3PL"
    cutpoint = 0.05,        # 0.05 (default) or 0.10
    score_method = "EAP",   # "EAP", "MAP", or "ML"
    seed = 42,              # seed for cutoff() resampling
    plot = TRUE             # produce plots
){
  if (!requireNamespace("PerFit", quietly = TRUE)) {
    stop("Package 'PerFit' is required.")
  }

  # --- Helper functions for argument matching --------------------------------
  .match_method <- function(x){
    x_lower <- tolower(trimws(x[1]))
    # Direct alias mapping
    if (x_lower %in% c("parametric", "para", "par", "p", "param")) return("parametric")
    if (x_lower %in% c("nonparametric", "nonpar", "non", "np", "nonparam")) return("nonparametric")
    if (x_lower %in% c("both", "all", "b")) return("both")

    stop(sprintf("Invalid method '%s'. Must be 'parametric', 'nonparametric', or 'both'.", x[1]), call. = FALSE)
  }

  .match_model <- function(x){
    x_lower <- tolower(trimws(x[1]))
    # Direct alias mapping
    if (x_lower %in% c("rasch", "1pl")) return("Rasch")
    if (x_lower == "2pl") return("2PL")
    if (x_lower == "3pl") return("3PL")
    if (x %in% c("Rasch", "2PL", "3PL")) return(x)  # Exact matches

    stop(sprintf("Invalid model '%s'. Must be 'Rasch', '2PL', or '3PL'.", x[1]), call. = FALSE)
  }

  .match_score_method <- function(x){
    x_upper <- toupper(trimws(x[1]))
    # Direct mapping for score methods
    if (x_upper %in% c("EAP", "MAP", "ML")) return(x_upper)
    if (tolower(x_upper) %in% c("eap", "map", "ml")) return(toupper(x_upper))

    stop(sprintf("Invalid score_method '%s'. Must be 'EAP', 'MAP', or 'ML'.", x[1]), call. = FALSE)
  }

  .irt_pmodel <- function(model){
    switch(model, Rasch = "1PL", `2PL` = "2PL", `3PL` = "3PL")
  }

  # --- Process arguments ------------------------------------------------------
  method <- .match_method(method)
  model <- .match_model(model)
  score_method <- .match_score_method(score_method)

  if (!is.numeric(cutpoint) || !(cutpoint %in% c(0.05, 0.1, 0.10))) {
    stop("cutpoint must be 0.05 or 0.10.")
  }
  cutpoint <- ifelse(abs(cutpoint - 0.1) < .Machine$double.eps^0.5, 0.10, cutpoint)

  # data validation
  invalid_values <- validate_data(data, c(0,1))
  if (!is.null(invalid_values)) {
    stop(paste0("Data contains ", nrow(invalid_values), " invalid value(s). Only 0, 1, and NA are allowed.\n",
                "Invalid locations:\n",
                paste(capture.output(print(invalid_values)), collapse = "\n")),
         call. = FALSE)
  }

  # --- Compute person-fit statistics for plotting ----------------------------
  pmodel <- .irt_pmodel(model)
  lzstar_stat <- lzstar_cut <- Ht_stat <- Ht_cut <- NULL

  # Parametric statistics (lz*)
  if (method %in% c("parametric", "both")) {
    if (!requireNamespace("mirt", quietly = TRUE)) {
      stop("Package 'mirt' is required for parametric person-fit analysis.")
    }

    N <- nrow(data); J <- ncol(data)
    idx_perfect <- rowSums(data, na.rm = TRUE) == J
    n_perfect <- sum(idx_perfect, na.rm = TRUE)
    if (n_perfect > 0) {
      message(sprintf("Parametric step: removing %d perfect scorers before lz* computation.", n_perfect))
    }

    data_np <- data[!idx_perfect, , drop = FALSE]

    # Fit IRT model on non-perfect cases only
    itemtype <- switch(model, Rasch = "Rasch", `2PL` = "2PL", `3PL` = "3PL")
    fit <- mirt::mirt(data_np, 1, itemtype = itemtype, verbose = FALSE)
    ip <- coef(fit, simplify = TRUE, IRTpars = TRUE)$item[, c("a","b","g")]
    ability <- as.numeric(mirt::fscores(fit, method = score_method, full.scores.SE = FALSE))

    # Compute lz* statistic
    lzstar_stat <- PerFit::lzstar(matrix = data_np, IP = ip, IRT.PModel = pmodel, Ability = ability)
    set.seed(seed)
    lzstar_cut <- PerFit::cutoff(lzstar_stat, Blvl = cutpoint)
  }

  # Nonparametric statistics (Ht)
  if (method %in% c("nonparametric", "both")) {
    Ht_stat <- PerFit::Ht(data)
    set.seed(seed)
    Ht_cut <- PerFit::cutoff(Ht_stat, Blvl = cutpoint)
  }

  # --- Create plots -----------------------------------------------------------
  if (isTRUE(plot)) {
    # Save and restore graphics parameters to prevent interference
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    if (method == "both") {
      graphics::par(mfrow = c(1, 2))
      plot(lzstar_stat, cutoff.obj = lzstar_cut, Type = "Histogram",
           main = "lz* PF scores + cutoff")
      plot(Ht_stat, cutoff.obj = Ht_cut, Type = "Histogram",
           main = "Ht PF scores + cutoff")
    } else if (method == "parametric") {
      plot(lzstar_stat, cutoff.obj = lzstar_cut, Type = "Histogram",
           main = "lz* PF scores + cutoff")
    } else { # nonparametric
      plot(Ht_stat, cutoff.obj = Ht_cut, Type = "Histogram",
           main = "Ht PF scores + cutoff")
    }
  }

  invisible(NULL)  # Function returns nothing, focuses purely on plotting
}
