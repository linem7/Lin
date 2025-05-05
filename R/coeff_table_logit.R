#' Create a formatted coefficient table for logistic and multinomial models
#'
#' `coeff_table_logit()` takes a fitted model—either a `glm(..., family=binomial)` or
#' a `nnet::multinom()`—and returns an HTML table of coefficients (b), standard errors (SE),
#' odds ratios (OR), p-values, optional confidence intervals (CI), and significance stars.
#' Categorical predictors are indented; very small numbers display as “< 0.01” or “< 0.001”.
#'
#' @param model A fitted model object.  Must be either:
#'   - A `glm` with `family = binomial(link = "logit")`, or
#'   - A `multinom` object from **nnet**.
#' @param digits Integer ≥ 0. Number of decimal places for b, SE, and OR.
#'   Values nonzero but smaller than `1/10^digits` are shown as `< 0.xxx`.
#' @param p_digits Integer ≥ 0. Decimal places for p-values; values smaller than
#'   `1/10^p_digits` are shown as `< 0.xxx`.
#' @param ci Logical; if `TRUE`, include 95% confidence intervals for the OR
#'   computed via the Wald formula on the log-odds scale (SPSS style).
#' @param stars Logical; if `TRUE`, append `*` for p<.05, `**` for p<.01,
#'   and `***` for p<.001 to the b values.
#' @param notes Optional character.  A footnote printed below the table.
#' @param title Optional character.  The table caption.
#' @param ... Additional arguments passed to `knitr::kable()`.
#'
#' @details
#' By default, R uses the *first* factor level as the reference for both the outcome
#' (in `glm`) and any categorical predictors (`contr.treatment`).  SPSS by default
#' uses the *last* level.  As a result, the intercept (constant) may differ
#' numerically between R and SPSS.  To match SPSS’s intercept, relevel your outcome
#' and predictors in R, e.g.
#' ```r
#' df$Y <- relevel(df$Y, ref = "1")        # if SPSS models “1” as success
#' df$gender <- relevel(df$gender, ref = "Female")
#' ```
#'
#' @return An HTML `kable` object (styled via **kableExtra**) showing b, SE, OR,
#'   p value, optional CI, and significance stars, with one block of columns per
#'   non-baseline outcome level.
#'
#' @examples
#' # Binary logistic regression example (mtcars) ------------------------------
#' library(dplyr)
#' data(mtcars)
#' mtcars$am <- factor(mtcars$am,
#'                     levels = c(0,1),
#'                     labels = c("Automatic","Manual"))
#' glm_fit <- glm(am ~ mpg + factor(cyl) + hp,
#'                data   = mtcars,
#'                family = binomial(link = "logit"))
#' coeff_table_logit(
#'   glm_fit,
#'   digits    = 2,
#'   p_digits  = 3,
#'   ci        = TRUE,
#'   stars     = TRUE,
#'   notes     = "Standard errors in parentheses.",
#'   title     = "Logistic Regression on mtcars"
#' )
#'
#' # Multinomial logistic regression example (iris) ----------------------------
#' library(nnet)
#' library(dplyr)
#' data(iris)
#' iris2 <- iris %>%
#'   mutate(WidthCat = cut(Sepal.Width,
#'                         breaks = 3,
#'                         labels = c("Small","Medium","Large")))
#' # match SPSS reference for Species = "setosa"
#' iris2$Species <- relevel(iris2$Species, ref = "setosa")
#' multinom_fit <- multinom(Species ~ Sepal.Length + Petal.Width + WidthCat,
#'                          data  = iris2,
#'                          trace = FALSE)
#' coeff_table_logit(
#'   multinom_fit,
#'   digits    = 3,
#'   p_digits  = 3,
#'   ci        = TRUE,
#'   stars     = FALSE,
#'   notes     = "Odds ratios and 95% CIs via Wald method.",
#'   title     = "Iris Species Multinomial Model"
#' )
#'
#' @importFrom broom tidy
#' @importFrom dplyr select mutate left_join case_when across
#' @importFrom tidyr pivot_wider crossing
#' @importFrom knitr kable
#' @importFrom kableExtra kable_classic add_header_above footnote
#' @export
coeff_table_logit <- function(model,
                              digits     = 3,
                              p_digits   = 3,
                              ci         = FALSE,
                              stars      = TRUE,
                              notes      = NULL,
                              title      = NULL,
                              ...) {
  # 0. allow multinom or glm(binomial)
  ok1 <- inherits(model, "multinom")
  ok2 <- inherits(model, "glm") && family(model)$family == "binomial"
  if (!(ok1 || ok2)) {
    stop("`model` must be a multinom object or a glm(binomial) object.")
  }

  # 1. recover data & model.frame
  df_data <- eval(model$call$data, envir = parent.frame())
  mf      <- model.frame(formula(model), data = df_data)

  # 2. outcome levels (drops baseline)
  out_lvls <- levels(stats::model.response(mf))[-1]

  # 3. formatting helpers
  fmt_num <- function(x, digs) {
    th     <- 1 / (10^digs)
    th_str <- formatC(th, format="f", digits=digs)
    vapply(x, function(v) {
      if (is.na(v)) return("")
      if (v != 0 && abs(v) < th) return(paste0("< ", th_str))
      formatC(v, format="f", digits=digs)
    }, "")
  }
  fmt_p <- function(x, digs) {
    th     <- 1 / (10^digs)
    th_str <- formatC(th, format="f", digits=digs)
    vapply(x, function(v) {
      if (is.na(v)) return("")
      if (v < th)   return(paste0("< ", th_str))
      formatC(v, format="f", digits=digs)
    }, "")
  }

  # 4. tidy coefficient table
  est_tbl <- broom::tidy(model, exponentiate = FALSE, conf.int = FALSE)
  if (!"y.level" %in% names(est_tbl)) {
    est_tbl <- est_tbl %>% mutate(y.level = out_lvls)
  }

  # 5. tidy OR + CI table
  or_tbl <- broom::tidy(model,
                        exponentiate = TRUE,
                        conf.int     = TRUE,
                        conf.level   = 0.95) %>%
    rename(or_raw  = estimate,
           CI_low  = conf.low,
           CI_high = conf.high)
  if (!"y.level" %in% names(or_tbl)) {
    or_tbl <- or_tbl %>% mutate(y.level = out_lvls)
  }

  # 6. merge & format stats
  stats_df <- est_tbl %>%
    select(y.level, term,
           b_raw  = estimate,
           se_raw = std.error,
           p_raw  = p.value) %>%
    left_join(or_tbl %>% select(y.level, term, or_raw, CI_low, CI_high),
              by = c("y.level", "term")) %>%
    mutate(
      b_fmt  = fmt_num(b_raw,  digits),
      se_fmt = fmt_num(se_raw, digits),
      or_fmt = fmt_num(or_raw, digits),
      p_fmt  = fmt_p(p_raw,  p_digits),
      CI_fmt = if (ci) {
        lo <- formatC(CI_low,  format="f", digits=digits)
        hi <- formatC(CI_high, format="f", digits=digits)
        vapply(seq_along(lo), function(i) {
          if (is.na(lo[i])||is.na(hi[i])) return("")
          paste0("[", lo[i], ", ", hi[i], "]")
        }, "")
      } else NULL,
      star   = case_when(
        !stars       ~ "",
        is.na(p_raw) ~ "",
        p_raw < .001 ~ "***",
        p_raw < .01  ~ "**",
        p_raw < .05  ~ "*",
        TRUE         ~ ""
      ),
      b_disp = ifelse(b_fmt %in% c("", "-"), b_fmt, paste0(b_fmt, star))
    ) %>%
    { # dynamically select only existing columns
      cols <- c("y.level", "term", "b_disp", "se_fmt", "or_fmt", "p_fmt")
      if (ci) cols <- c(cols, "CI_fmt")
      select(., all_of(cols))
    }

  # 7. build display structure: intercept + predictors
  terms_all    <- names(attr(terms(model), "dataClasses"))[-1]
  xlevels      <- model$xlevels

  # intercept row
  intercept_df <- tibble(
    display = "(Intercept)",
    term    = "(Intercept)",
    type    = "est"
  )

  # predictor rows
  predictor_df <- bind_rows(lapply(terms_all, function(var) {
    lab  <- attr(mf[[var]], "label", exact = TRUE)
    disp <- if (!is.null(lab)) lab else tools::toTitleCase(var)
    if (var %in% names(xlevels)) {
      levs <- xlevels[[var]]
      tibble(
        display = c(disp, paste0("\u00A0\u00A0", levs)),
        term    = c(NA, NA, paste0(var, levs[-1])),
        type    = c("label", "ref", rep("est", length(levs)-1))
      )
    } else {
      tibble(display = disp, term = var, type = "est")
    }
  }))

  struct <- bind_rows(intercept_df, predictor_df) %>%
    mutate(ord = row_number())

  # 8. merge struct + stats, hyphens/blanks
  stat_cols <- c("b_disp", "se_fmt", "or_fmt", "p_fmt")
  if (ci) stat_cols <- c(stat_cols, "CI_fmt")

  combo <- tidyr::crossing(struct, y.level = out_lvls) %>%
    left_join(stats_df, by = c("y.level", "term")) %>%
    mutate(across(all_of(stat_cols), ~ case_when(
      type == "label" ~ "",
      type == "ref"   ~ "\u2013",
      TRUE            ~ .
    ))) %>%
    arrange(ord, factor(y.level, levels = out_lvls))

  # 9. pivot to wide
  wide <- combo %>%
    select(-type, -term, -ord) %>%
    pivot_wider(
      id_cols     = display,
      names_from  = y.level,
      values_from = all_of(stat_cols),
      names_glue  = "{y.level}_{.value}"
    ) %>%
    rename(Predictor = display)

  # 10. reorder columns
  suffixes <- c("_b_disp", "_se_fmt", "_or_fmt", "_p_fmt")
  if (ci) suffixes <- c(suffixes, "_CI_fmt")
  col_seq  <- unlist(lapply(out_lvls, function(l) paste0(l, suffixes)))
  wide     <- wide %>% select(Predictor, any_of(col_seq))

  # 11. headers & render
  stats_names <- c("b", "SE", "OR", "p value")
  if (ci) stats_names <- c(stats_names, "CI")
  bottom      <- c("Predictor", rep(stats_names, length(out_lvls)))
  top         <- c(" " = 1,
                   setNames(rep(length(stats_names), length(out_lvls)),
                            out_lvls))

  wide %>%
    knitr::kable(
      col.names = bottom,
      caption   = title,
      align     = c("l", rep("c", ncol(.)-1)),
      escape    = FALSE,
      ...
    ) %>%
    kableExtra::kable_classic(
      full_width = FALSE,
      html_font  = "Cambria",
      position   = "left"
    ) %>%
    kableExtra::add_header_above(top) %>%
    kableExtra::footnote(
      general           = notes,
      footnote_as_chunk = TRUE,
      escape            = FALSE
    )
}

