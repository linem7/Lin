#' Compute CFA Reliability Table
#'
#' Fits one or more lavaan‐style CFA specifications, then for each latent:
#' - Cronbach’s α (tau-equivalent)
#' - Composite reliability (CR, tau-unequal)
#' - Average variance extracted (AVE)
#'
#' If your spec has more than one factor, it also builds a single‐factor “Total”
#' model (all indicators) and returns its α/CR/AVE.
#'
#' @param ...   Unnamed character strings of lavaan model syntax (each may
#'              contain one or more “=~” lines).
#' @param data  A data.frame of observed variables.
#' @param digits Number of decimal places; uses \code{\link{sprintf}} so all
#'               outputs are character.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{Variable}{The object name you passed (e.g. “spec_sf”).}
#'     \item{Domain}{Each latent name (and “Total” if multi-factor).}
#'     \item{Alpha}{Formatted Cronbach’s α.}
#'     \item{CR}{Formatted composite reliability.}
#'     \item{AVE}{Formatted AVE.}
#'   }
#'
#' @examples
#' ## Single-factor CFA
#' spec_sf <- 'F1 =~ x1 + x2 + x3 + x4'
#' rel_sf  <- reliability_table(spec_sf, data = my_data, digits = 3)
#'
#' ## Multi-factor CFA
#' spec_mf <- paste0(
#'   'F1 =~ x1 + x2 + x3', '\\n',
#'   'F2 =~ y1 + y2 + y3 + y4'
#' )
#' rel_mf <- reliability_table(spec_mf, data = my_data, digits = 2)
#'
#' ## Pass multiple specs at once
#' rel_all <- reliability_table(spec_sf, spec_mf, data = my_data)
#'
#' @export
reliability_table <- function(..., data, digits = 3) {
  # capture all specs and the names you used when calling the function
  specs     <- list(...)
  var_names <- as.character(substitute(list(...)))[-1]

  # for each spec + its name, build a little tibble
  results <- purrr::map2_dfr(specs, var_names, function(spec, var_name) {
    # --- 1) break your model text into lines, keep only the "=~" lines
    lines        <- stringr::str_split(spec, pattern = '\\n')[[1]]
    domain_lines <- lines[grepl('=~', lines)]

    # --- 2) split each "latent =~ indicators" line into LHS & RHS
    split_mat    <- stringr::str_split_fixed(
      string  = domain_lines,
      pattern = '=~',
      n       = 2
    )
    latent_names <- stringr::str_trim(split_mat[,1])
    rhs_parts    <- stringr::str_trim(split_mat[,2])

    # --- 3) fit the CFA once
    fit <- lavaan::cfa(spec, data = data)

    # --- 4) compute α, CR, AVE for each latent
    rel_te   <- semTools::compRelSEM(fit, tau.eq   = TRUE)
    rel_tu   <- semTools::compRelSEM(fit, tau.eq   = FALSE)
    ave_vals <- semTools::AVE(fit)

    alpha_v <- rel_te[latent_names]
    cr_v    <- rel_tu[latent_names]
    ave_v   <- ave_vals[latent_names]

    # --- 5) format as strings
    fmt       <- function(x) sprintf(paste0('%.', digits, 'f'), x)
    alpha_fmt <- fmt(alpha_v)
    cr_fmt    <- fmt(cr_v)
    ave_fmt   <- fmt(ave_v)

    # --- 6) assemble domain‐level tibble
    df <- tibble::tibble(
      Variable = var_name,
      Domain   = latent_names,
      Alpha    = alpha_fmt,
      CR       = cr_fmt,
      AVE      = ave_fmt
    )

    # --- 7) if multi‐factor, also do a “Total” one-factor CFA
    if (length(latent_names) > 1) {
      # gather all indicators
      items_all <- unlist(strsplit(rhs_parts, split = ' \\+ '))
      items_all <- stringr::str_trim(items_all)
      items_all <- items_all[items_all != '']

      # build & fit Total model
      total_spec <- paste0('Total =~ ', paste(items_all, collapse = ' + '))
      fit_tot    <- lavaan::cfa(total_spec, data = data)

      # compute its reliabilities
      rel_te_tot  <- semTools::compRelSEM(fit_tot, tau.eq = TRUE)
      rel_tu_tot  <- semTools::compRelSEM(fit_tot, tau.eq = FALSE)
      ave_tot     <- semTools::AVE(fit_tot)

      alpha_tot   <- rel_te_tot['Total']
      cr_tot      <- rel_tu_tot['Total']
      ave_tot_val <- ave_tot['Total']

      alpha_tot_f <- fmt(alpha_tot)
      cr_tot_f    <- fmt(cr_tot)
      ave_tot_f   <- fmt(ave_tot_val)

      df_tot <- tibble::tibble(
        Variable = var_name,
        Domain   = 'Total',
        Alpha    = alpha_tot_f,
        CR       = cr_tot_f,
        AVE      = ave_tot_f
      )
      df <- dplyr::bind_rows(df, df_tot)
    }

    df
  })

  # --- 8) blank out repeated Variable names for readability
  results %>%
    dplyr::mutate(Variable = ifelse(duplicated(Variable), '', Variable))
}
