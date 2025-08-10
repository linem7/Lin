#' Plot Person Response Functions (PRFs) for the most (in)consistent respondents
#'
#' @description
#' Rank persons by a person-fit statistic (PFS) and draw their smoothed
#' person response functions (PRFs) with `ggplot2` (one facet per person).
#' Supports nonparametric `Ht` and parametric `lz` (requires item parameters
#' and/or abilities).
#'
#' @details
#' **Workflow.**
#' 1) Validate `X` (must be 0/1/`NA`).
#' 2) Compute PFS (`Ht` or `lz`) on all rows.
#' 3) If PFS fails or yields `NA`s, drop rows with *effective perfect* scores
#'    (all observed 1s or all observed 0s) or with < 2 observed items and recompute.
#' 4) Select `n_cases` highest/lowest per `direction`, build PRFs via `PerFit::PRFplot()`,
#'    evaluate on an equally spaced grid in \[0,1\], and facet with `ggplot2`.
#'
#' **X-axis meaning (how to read PRFs).**
#' - The PRF x-axis is a **standardized item difficulty rank** in \[0,1\].
#' - Let \(p_j = \Pr(\text{correct on item } j)\) computed as the empirical proportion
#'   correct (ignoring `NA`). Items are **ordered by increasing \(p_j\)** (hard → easy).
#' - If an item’s rank among \(J\) items is \(r \in \{1,\dots,J\}\), its standardized
#'   position is \(x = (r-1)/(J-1)\). Thus **x=0 ≈ hardest** items; **x=1 ≈ easiest**.
#' - Ties in \(p_j\) receive averaged ranks; the PRF is then **smoothed** over this
#'   rank axis (controlled by `h` and `N.FPts`).
#' - This x-axis construction is used for the PRF regardless of whether you ranked
#'   persons via `Ht` or `lz`; the **method only affects which persons are shown**, not
#'   how the PRF x-axis is defined.
#'
#' **Y-axis meaning.**
#' - The y-axis is the **smoothed probability of a correct response** across the
#'   difficulty-ranked items. A more monotone increasing curve signals response
#'   patterns consistent with item difficulty; dips/oscillations suggest local
#'   inconsistencies.
#'
#' **Interpretation tips.**
#' - Left side (x≈0): hardest items; right side (x≈1): easiest.
#' - **Steeper rise** indicates more consistent improvement from hard to easy items.
#' - **Dips** at intermediate x suggest unexpected misses on relatively easier items.
#' - The **area under the PRF** roughly tracks a person’s overall raw score (higher
#'   area ↔ more correct responses), but PRFs are primarily for visual **consistency**,
#'   not score estimation.
#'
#' **Parametric (`lz`).**
#' - Provide `IP` and/or `Ability`. A convenient route is to fit a 1-factor model with
#'   **mirt** and pass `IP = coef(fit, IRTpars = TRUE, simplify = TRUE)$items[, c("a","b","g")]`
#'   (use `g = 0` for 2PL). `Ability` can be omitted—`PerFit::lz()` will estimate abilities if needed.
#'
#' @param X Respondent × item data (0/1/`NA`) as matrix or data.frame.
#' @param method `"nonpar"` (default) or `"par"`. Sets the default `pfs`.
#' @param pfs `"Ht"` or `"lz"`. If `NULL`, uses `"Ht"` for `nonpar`, `"lz"` for `par`.
#' @param direction `"worse"` (lowest PFS first) or `"good"` (highest first).
#' @param n_cases Number of persons to display (default 10).
#' @param h Bandwidth for PRF smoothing (passed to `PerFit::PRFplot()`; default .09).
#' @param N.FPts Number of evaluation points in \[0,1\] (default 101).
#' @param IP Optional item parameter matrix for `lz` (e.g., columns `a`, `b`, `g`).
#' @param Ability Optional ability vector for `lz` (length `nrow(X)`).
#' @param facet_ncol Columns in `facet_wrap()` (default 5).
#'
#' @returns A list with:
#' * `plot` — a `ggplot` PRF panel.
#' * `data` — tibble with columns `x` (grid), `id` (row index of `X`), `prf`.
#' * `id` — integer vector of selected row indices.
#' * `raw` — raw 0/1 rows for those `id`s (`orig_id` column + items).
#'
#' @examples
#' # --- Nonparametric (Ht): quick demo on PerFit's data -----------------------
#' data(PerFit::InadequacyData)
#' res_ht <- gg_prf(PerFit::InadequacyData, method = "nonpar", pfs = "Ht",
#'                  n_cases = 9, facet_ncol = 3)
#' print(res_ht$plot)
#' head(res_ht$data)
#' res_ht$id
#' head(res_ht$raw)
#'
#' # --- Parametric (lz) via mirt: obtain IP (and optionally Ability) ----------
#' if (requireNamespace("mirt", quietly = TRUE)) {
#'   library(mirt)
#'   X <- PerFit::InadequacyData
#'   # Fit a unidimensional model (2PL by default)
#'   fit <- mirt(X, 1)
#'
#'   # Item parameters in IRT metric; use columns a, b, g (g=0 under 2PL)
#'   ip <- coef(fit, IRTpars = TRUE, simplify = TRUE, digits = Inf)$items[, c("a","b","g")]
#'
#'   # Optional abilities (EAP); 'lz' can estimate abilities if you omit this
#'   th <- fscores(fit, method = "EAP")[, 1]
#'
#'   res_lz <- gg_prf(X, method = "par", pfs = "lz",
#'                    IP = as.matrix(ip), Ability = as.numeric(th),
#'                    n_cases = 6, facet_ncol = 3)
#'   print(res_lz$plot)
#' }
#'
#' @seealso [PerFit::Ht()], [PerFit::lz()], [PerFit::PRFplot()], [mirt::mirt()]
#' @importFrom PerFit Ht lz PRFplot
#' @importFrom fda eval.fd
#' @importFrom ggplot2 ggplot aes geom_line facet_wrap scale_x_continuous scale_y_continuous labs theme_classic
#' @importFrom tibble tibble
#' @importFrom grDevices pdf dev.off
#' @export
gg_prf <- function(
    data,
    method = c("parametric", "nonparametric"),
    model_type = c("2PL", "3PL", "Rasch"),
    direction = c("worse", "good"),
    n_cases = 10,
    facet_ncol = 5,
    ...                          # additional arguments passed to PRFplot
){

  # --- Helper functions for argument matching --------------------------------
  .match_method <- function(x){
    x_lower <- tolower(trimws(x[1]))
    # Direct alias mapping
    if (x_lower %in% c("parametric", "para", "par", "p", "param")) return("parametric")
    if (x_lower %in% c("nonparametric", "nonpar", "non", "np", "nonparam")) return("nonparametric")

    stop(sprintf("Invalid method '%s'. Must be 'parametric' or 'nonparametric'.", x[1]), call. = FALSE)
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

  # --- Argument validation ----------------------------------------------------
  method <- .match_method(method)
  model_type <- .match_model(model_type)
  direction <- match.arg(direction)

  # --- Dependency checks ------------------------------------------------------
  if (!requireNamespace("PerFit", quietly = TRUE)) stop("Need PerFit package.")
  if (!requireNamespace("fda", quietly = TRUE))    stop("Need fda package.")
  if (!requireNamespace("ggplot2", quietly = TRUE))stop("Need ggplot2 package.")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Need tibble package.")

  # --- Validate data: must be numeric/integer/logical and only {0,1,NA} ------
  check_binary_df <- function(df) {
    bad_type <- !vapply(df, function(x) is.numeric(x) || is.integer(x) || is.logical(x), logical(1))
    if (any(bad_type)) {
      offending <- names(df)[bad_type]
      stop("data has non-numeric columns (e.g., factor/character): ", paste(offending, collapse = ", "))
    }
    bad_val <- vapply(df, function(x) any(!is.na(x) & !(x %in% c(0,1))), logical(1))
    if (any(bad_val)) {
      offending <- names(df)[bad_val]
      stop("data has values not in {0,1,NA} in columns: ", paste(offending, collapse = ", "))
    }
    invisible(TRUE)
  }
  check_binary_df(data)

  N <- nrow(data); J <- ncol(data)
  if (N == 0 || J == 0) stop("data has zero rows or columns after preprocessing.")

  # --- Determine PFS statistic based on method -------------------------------
  pfs <- if (method == "parametric") "lz" else "Ht"

  # --- Auto-estimate parameters if needed for parametric method --------------
  IP <- NULL
  Ability <- NULL

  if (method == "parametric") {
    if (!requireNamespace("mirt", quietly = TRUE)) {
      stop("Need mirt package for parametric method. Install it to proceed.")
    }

    # Fit IRT model with SE estimation
    mod <- suppressMessages(
      mirt::mirt(data, model = 1, itemtype = model_type, verbose = FALSE)
    )

    if (inherits(mod, "try-error")) {
      stop("Failed to estimate IRT model. Check your data or try nonparametric method.")
    }

    # Extract item parameters using IRTpars=TRUE for direct IRT parameterization
    item_params <- mirt::coef(mod, simplify = TRUE, IRTpars = TRUE)$items

    # Format IP matrix for PerFit: always needs 3 columns [a, b, c]
    # Take first 3 columns - mirt returns a, b, g for all models (g=0 for Rasch/2PL)
    IP <- as.matrix(item_params[, 1:3])
    colnames(IP) <- c('a', 'b', 'c')

    # Extract person abilities using EAP with full scores
    person_scores <- mirt::fscores(mod,
                                   method = "EAP",
                                   full.scores = TRUE,
                                   full.scores.SE = TRUE)
    Ability <- as.vector(person_scores[, 1])
  }

  # --- Robust PFS extractor ---------------------------------------------------
  .extract_pfs <- function(pf_obj) {
    s <- pf_obj$PFscores
    if (is.null(s)) stop("PFscores missing in PerFit object.")

    if (is.data.frame(s)) {
      s <- if ("PFscores" %in% names(s)) s[["PFscores"]] else s[[1]]
    }
    if (is.matrix(s)) {
      s <- s[, 1, drop = TRUE]
    }
    if (is.list(s)) {
      s <- unlist(s, use.names = FALSE, recursive = TRUE)
    }

    s <- suppressWarnings(as.numeric(s))
    s
  }

  # --- Define fallback filter -------------------------------------------------
  ones  <- rowSums(data == 1, na.rm = TRUE)
  zeros <- rowSums(data == 0, na.rm = TRUE)
  obs   <- ones + zeros

  effective_perfect <- (obs > 0) & ((ones == obs) | (zeros == obs))
  too_few           <- obs < 2
  drop_mask         <- effective_perfect | too_few

  # --- Compute PFS scores -----------------------------------------------------
  compute_pfs <- function(M_df, pfs, IP = NULL, Ability = NULL) {
    if (pfs == "Ht") {
      PerFit::Ht(M_df)
    } else {
      PerFit::lz(M_df, IP = IP, Ability = Ability)
    }
  }

  # Try computing PFS on full data
  kept_index <- seq_len(N)
  pf_obj <- try(compute_pfs(data, pfs = pfs, IP = IP, Ability = Ability), silent = TRUE)
  scores <- if (!inherits(pf_obj, "try-error")) .extract_pfs(pf_obj) else NULL

  # If failed, retry with filtered data
  if (inherits(pf_obj, "try-error") || is.null(scores) || anyNA(scores)) {
    kept_index <- which(!drop_mask)
    if (length(kept_index) == 0L) stop("No eligible rows after filtering.")
    data_wo <- data[kept_index, , drop = FALSE]
    Ability_wo <- if (!is.null(Ability)) Ability[kept_index] else NULL
    pf_obj <- compute_pfs(data_wo, pfs = pfs, IP = IP, Ability = Ability_wo)
    scores <- .extract_pfs(pf_obj)
  }

  if (!is.numeric(scores) || length(scores) != length(kept_index) || anyNA(scores)) {
    stop("PFscores malformed or contain NA after extraction.")
  }

  # --- Select cases based on direction ---------------------------------------
  ord_local <- if (direction == "worse") order(scores) else order(-scores)
  sel_local <- head(ord_local, n_cases)
  respID    <- kept_index[sel_local]

  # --- Build PRFs using PRFplot ----------------------------------------------
  prf <- local({
    tmp <- tempfile(fileext = ".pdf")
    grDevices::pdf(file = tmp)
    on.exit({
      try(grDevices::dev.off(), silent = TRUE)
      try(unlink(tmp), silent = TRUE)
    }, add = TRUE)
    suppressWarnings(suppressMessages(
      PerFit::PRFplot(
        matrix   = data,
        respID   = respID,
        VarBands = FALSE,
        message  = FALSE,
        ...                    # pass additional arguments
      )
    ))
  })

  # --- Extract and process PRF data ------------------------------------------
  # Default N.FPts is 101 in PRFplot if not specified in ...
  dots <- list(...)
  N.FPts <- if ("N.FPts" %in% names(dots)) dots$N.FPts else 101

  grid  <- seq(0, 1, length.out = N.FPts)
  y_hat <- fda::eval.fd(grid, prf$PRF.FDO)
  if (is.null(dim(y_hat))) y_hat <- matrix(y_hat, ncol = 1)

  # Align columns to selected IDs
  if (ncol(y_hat) != length(respID)) {
    if (ncol(y_hat) == nrow(data)) {
      y_hat <- y_hat[, respID, drop = FALSE]
    } else {
      y_hat <- y_hat[, seq_along(respID), drop = FALSE]
    }
  }

  # Clamp to [0,1]
  prf_vec <- as.vector(y_hat)
  prf_vec <- pmin(pmax(prf_vec, 0), 1)

  # --- Create long-format data for plotting ----------------------------------
  df <- tibble::tibble(
    x   = rep(grid, times = ncol(y_hat)),
    id  = rep(respID, each = length(grid)),
    prf = prf_vec
  )
  df <- df[is.finite(df$prf) & is.finite(df$x), , drop = FALSE]

  # --- Generate plot ----------------------------------------------------------
  method_label <- if (method == "parametric") "Parametric" else "Nonparametric"
  pfs_label    <- if (pfs == "lz") "lz" else "Ht"
  title_txt    <- sprintf("Person Response Functions for top %d %s cases",
                          length(unique(df$id)), direction)
  subtitle_txt <- sprintf("(%s method: %s statistic)", method_label, pfs_label)

  # Conditional x-axis label based on method
  x_label <- if (method == "parametric") {
    "Item difficulty (IRT-based, 0-1)"
  } else {
    "Empirical difficulty rank (0-1)"
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = prf, group = id)) +
    ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
    ggplot2::facet_wrap(~ id, ncol = facet_ncol) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
    ggplot2::labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = x_label,
      y = "Smoothed P(correct)"
    ) +
    ggplot2::theme_classic(base_size = 12)

  # --- Return results ---------------------------------------------------------
  list(
    plot = p,
    data = df,
    id   = respID,
    raw  = cbind(orig_id = respID, data[respID, , drop = FALSE])
  )
}
