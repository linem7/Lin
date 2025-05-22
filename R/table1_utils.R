#' @title Statistical Test Formatter for table1()
#' @description
#' Formats and returns the test statistic and degrees of freedom (with significance stars)
#' for comparisons across groups in **table1**.
#'
#' @details
#' \enumerate{
#'   \item \strong{Variance checks:}
#'     \itemize{
#'       \item For two groups, runs \code{var.test(y ~ g)} to decide equal-variance vs. Welch t-test.
#'       \item For more than two groups, runs \code{car::leveneTest(y ~ g)}; if p < 0.05 uses \code{oneway.test()} (Welch ANOVA), otherwise \code{aov()}.
#'     }
#'   \item \strong{Categorical variables:}
#'     \itemize{
#'       \item Detects unused factor levels and suggests \code{droplevels()}.
#'       \item Performs \code{chisq.test()} on the contingency table.
#'     }
#'   \item \strong{Missing data & overall column:}
#'     \itemize{
#'       \item Drops \code{NA} values and any group that becomes empty.
#'       \item Strips any stratum whose name begins with “All” or “Overall…”.
#'     }
#'   \item \strong{Significance stars:}
#'     \itemize{
#'       \item \code{***} for p < 0.001
#'       \item \code{**}  for p < 0.01
#'       \item \code{*}   for p < 0.05
#'     }
#' }
#'
#' @note
#' To include a total column you must write:
#' \code{overall = "Overall"}
#' so that these helpers properly detect and ignore that column.
#'
#' @examples
#' \dontrun{
#' library(table1)
#' table1(
#'   ~ mpg | cyl,
#'   data      = mtcars,
#'   overall   = "Overall",
#'   extra.col = list(
#'     "P-value"    = table1_pvalue,
#'     "t / F / χ2" = table1_stat,
#'     "Effect size"= table1_ef
#'   )
#' )
#' }
#' @export
table1_stat <- function(x, name, ...) {
  # 0. ---- Handle total/overall column if present ----
  #    The first element of `x` will be the overall column when overall=TRUE.
  #    We detect that by its name being "" or "Overall" and drop it.
  if (!is.null(names(x))) {
    # grep for names beginning with "All" or "Overall" (case‐insensitive)
    overall_idx <- grep("^(All|Overall)", names(x), ignore.case = TRUE)
    if (length(overall_idx) > 0) {
      x <- x[-overall_idx]    # remove those elements before testing
    }
  }

  # 1. ---- Drop only true NAs from each group ----
  vals <- lapply(x, function(v) v[!is.na(v)])

  # 2. ---- Remove any groups that became empty ----
  vals <- vals[lengths(vals) > 0]

  # 3. ---- If fewer than 2 groups remain, warn and exit ----
  if (length(vals) < 2) {
    return(
      c("",
        paste0(
          "Warning: only one non‐missing group for ‘",
          name,
          "’"
        )
      )
    )
  }

  # 4. ---- Flatten data and rebuild grouping factor ----
  y  <- unlist(vals)
  g  <- factor(
    rep(seq_along(vals), times = lengths(vals)),
    labels = names(vals)
  )
  ng <- nlevels(g)

  # 5. ---- Initialize placeholders ----
  stat <- NA_real_
  df   <- NULL
  pval <- NA_real_

  # 6. ---- Numeric branch: t‐test or ANOVA/Welch ----
  if (is.numeric(y)) {
    if (ng == 2) {
      # 6a. Two groups → variance check + t‐test
      vt <- try(var.test(y ~ g), silent = TRUE)
      if (!inherits(vt, "try-error")) {
        tt <- try(
          t.test(y ~ g, var.equal = (vt$p.value > 0.05)),
          silent = TRUE
        )
        if (!inherits(tt, "try-error")) {
          stat <- as.numeric(tt$statistic)
          df   <- as.integer(tt$parameter)
          pval <- tt$p.value
        }
      }
    } else {
      # 6b. >2 groups → Levene’s test then ANOVA or Welch
      lt <- try(car::leveneTest(y ~ g), silent = TRUE)
      if (!inherits(lt, "try-error")) {
        if (lt$`Pr(>F)`[1] < 0.05) {
          # Welch ANOVA
          an <- try(oneway.test(y ~ g), silent = TRUE)
          if (!inherits(an, "try-error")) {
            stat <- as.numeric(an$statistic)
            df   <- as.integer(an$parameter)
            pval <- an$p.value
          }
        } else {
          # Classic ANOVA
          an <- try(aov(y ~ g), silent = TRUE)
          if (!inherits(an, "try-error")) {
            s    <- summary(an)[[1]]
            stat <- s["g", "F value"]
            df   <- c(s["g", "Df"], s["Residual", "Df"])
            pval <- s["g", "Pr(>F)"]
          }
        }
      }
    }

  } else {
    # 7. ---- Categorical branch: drop phantom levels + χ² ----
    orig_lvls    <- levels(x[[1]])
    present_lvls <- unique(as.character(unlist(vals)))
    unused_lvls  <- setdiff(orig_lvls, present_lvls)

    if (length(unused_lvls) > 0) {
      return(
        c("",
          paste0(
            "Warning: unused level(s) in ‘", name, "’: ",
            paste(unused_lvls, collapse = ", "),
            ". Use droplevels(data$", name, ")"
          )
        )
      )
    }

    tbl <- table(y, g)
    ct  <- try(chisq.test(tbl), silent = TRUE)
    if (!inherits(ct, "try-error")) {
      stat <- as.numeric(ct$statistic)
      df   <- as.integer(ct$parameter)
      pval <- ct$p.value
    }
  }

  # 8. ---- If test failed or pval is NA, return blank ----
  if (is.na(stat) || is.null(df) || is.na(pval)) {
    return(c("", ""))
  }

  # 9. ---- Star‐coding based on p-value ----
  stars <- if      (pval < 0.001) "***"
  else if (pval < 0.01)  "**"
  else if (pval < 0.05)  "*"
  else                   ""

  # 10. ---- Format stat and df ----
  stat_txt <- if (abs(stat) < 0.001) "< 0.001" else sprintf("%.3f", stat)
  df_txt   <- if (length(df) == 2) {
    sprintf("(%d, %d)", df[1], df[2])
  } else {
    sprintf("(%d)", df)
  }

  # 11. ---- Return for table1 alignment ----
  c("", paste0(stat_txt, " ", df_txt, " ", stars))
}


#' @title P-Value Extractor for table1()
#' @description
#' Extracts and formats the p-value for group comparisons in **table1**.
#'
#' @details
#' \enumerate{
#'   \item \strong{Test selection:}
#'     \itemize{
#'       \item Two groups → \code{t.test()} (variance checked via \code{var.test()}).
#'       \item >2 groups → \code{car::leveneTest()} then \code{oneway.test()} or \code{aov()}.
#'       \item Categorical → \code{chisq.test()}.
#'     }
#'   \item \strong{Missing & overall:}
#'     \itemize{
#'       \item Drops \code{NA}s and groups with zero observations.
#'       \item Strips any column named “All…” or “Overall…”.
#'     }
#'   \item \strong{Warnings:}
#'     \itemize{
#'       \item If only one non-missing group remains, prints a warning.
#'       \item If phantom factor levels remain, suggests \code{droplevels()}.
#'     }
#' }
#'
#' @note
#' Use \code{overall = "Overall"} in \code{table1()} to have the overall column recognized and dropped.
#'
#' @examples
#' \dontrun{
#' table1(
#'   ~ am | cyl,
#'   data      = mtcars,
#'   overall   = "Overall",
#'   extra.col = list("P-value" = table1_pvalue)
#' )
#' }
#' @export
table1_pvalue <- function(x, name, ...) {

  # 0. ─── Drop any stratum whose name starts with “All” or “Overall” ─────
  if (!is.null(names(x))) {
    # grep for names beginning with "All" or "Overall" (case‐insensitive)
    overall_idx <- grep("^(All|Overall)", names(x), ignore.case = TRUE)
    if (length(overall_idx) > 0) {
      x <- x[-overall_idx]    # remove those elements before testing
    }
  }

  # 1. Remove only true NAs from each group
  vals <- lapply(x, function(v) v[!is.na(v)])

  # 2. Drop groups that became empty
  vals <- vals[lengths(vals) > 0]

  # 3. If fewer than 2 groups remain, warn and exit
  if (length(vals) < 2) {
    return(
      c("",
        paste0(
          "Warning: only one non‐missing group for ‘",
          name, "’"
        )
      )
    )
  }

  # 4. Flatten data and rebuild grouping factor
  y <- unlist(vals)
  g <- factor(
    rep(seq_along(vals), times = lengths(vals)),
    labels = names(vals)
  )

  # 5. Initialize p
  p <- NA_real_

  # 6. Numeric branch: 2‐sample t or ANOVA/Welch depending on group count
  if (is.numeric(y)) {
    ng <- nlevels(g)

    if (ng == 2) {
      # 6a. Two groups → check variance
      vt <- try(var.test(y ~ g), silent = TRUE)
      if (!inherits(vt, "try-error")) {
        tt <- try(
          t.test(y ~ g, var.equal = (vt$p.value > 0.05)),
          silent = TRUE
        )
        if (!inherits(tt, "try-error")) p <- tt$p.value
      }

    } else {
      # 6b. >2 groups → Levene’s test then ANOVA or Welch
      lt <- try(car::leveneTest(y ~ g), silent = TRUE)
      if (!inherits(lt, "try-error")) {
        if (lt$`Pr(>F)`[1] < 0.05) {
          # Unequal variances → Welch ANOVA
          an <- try(oneway.test(y ~ g), silent = TRUE)
          if (!inherits(an, "try-error")) p <- an$p.value
        } else {
          # Homogeneous → classic ANOVA
          an <- try(aov(y ~ g), silent = TRUE)
          if (!inherits(an, "try-error")) {
            p <- summary(an)[[1]][ "g", "Pr(>F)" ]
          }
        }
      }
    }

  } else {
    # 7. Categorical branch: detect and warn on unused factor levels
    #    This checks the first group's levels attribute for unused levels
    orig_lvls    <- levels(x[[1]])
    present_lvls <- unique(as.character(unlist(vals)))
    unused_lvls  <- setdiff(orig_lvls, present_lvls)

    if (length(unused_lvls) > 0) {
      # 7a. Warn about phantom levels and suggest droplevels()
      return(
        c("",
          paste0(
            "Warning: unused level(s) in ‘", name, "’: ",
            paste(unused_lvls, collapse = ", "),
            ". Use droplevels(data$", name, ")"
          )
        )
      )
    }

    # 7b. Proceed with χ² test on present levels only
    tbl <- table(y, g)
    ct  <- try(chisq.test(tbl), silent = TRUE)
    if (!inherits(ct, "try-error")) p <- ct$p.value
  }

  # 8. If the test failed / p remains NA, return blank
  if (is.na(p)) return(c("", ""))

  # 9. Format p to 3 decimals (or "< 0.001")
  pv <- if (p < 0.001) "< 0.001" else sprintf("%.3f", p)

  # 10. Return two‐element vector for table1 alignment
  c("", pv)
}


#' @title Effect-Size Calculator for table1()
#' @description
#' Computes and formats effect sizes for group comparisons in **table1**.
#'
#' @details
#' \enumerate{
#'   \item \strong{Numeric comparisons:}
#'     \itemize{
#'       \item Two groups → Cohen’s \emph{d} (pooled SD).
#'       \item >2 groups → partial \emph{η}\^2 from ANOVA sums of squares.
#'     }
#'   \item \strong{Categorical comparisons:}
#'     \itemize{
#'       \item Cramer’s \emph{V} from \code{chisq.test()}.
#'       \item Warns on unused factor levels and suggests \code{droplevels()}.
#'     }
#'   \item \strong{Missing & overall:}
#'     \itemize{
#'       \item Drops \code{NA}s and empty groups.
#'       \item Strips any column named “All…” or “Overall…”.
#'     }
#' }
#'
#' @note
#' If you request \code{overall = "Overall"}, the helper will detect and remove it before computing effect sizes.
#'
#' @examples
#' \dontrun{
#' # Cohen’s d for mpg across 4 vs 6 cylinders
#' df <- subset(mtcars, cyl %in% c(4,6))
#' table1(
#'   ~ mpg | cyl,
#'   data      = df,
#'   overall   = "Overall",
#'   extra.col = list("Effect size" = table1_ef)
#' )
#' }
#' @export
table1_ef <- function(x, name, ...) {
  # 0. ---- Drop any “All”/“Overall” column (e.g. “Overall (N=200)”) ----
  if (!is.null(names(x))) {
    overall_idx <- grep("^(All|Overall)", names(x), ignore.case = TRUE)
    if (length(overall_idx) > 0) {
      x <- x[-overall_idx]
    }
  }

  # 1. ---- Remove only true NAs from each group ----
  vals <- lapply(x, function(v) v[!is.na(v)])

  # 2. ---- Drop groups that became empty ----
  vals <- vals[lengths(vals) > 0]

  # 3. ---- If fewer than 2 groups remain, warn and exit ----
  if (length(vals) < 2) {
    return(
      c("",
        paste0("Warning: only one non-missing group for ‘", name, "’")
      )
    )
  }

  # 4. ---- Flatten data and rebuild grouping factor ----
  y  <- unlist(vals)
  g  <- factor(
    rep(seq_along(vals), times = lengths(vals)),
    labels = names(vals)
  )
  ng <- nlevels(g)   # number of groups
  n  <- length(y)    # total non-missing N

  # 5. ---- Compute effect size ----
  es <- NA_real_
  if (is.numeric(y)) {
    if (ng == 2) {
      # 5a. Two groups → Cohen’s d (pooled SD)
      m   <- tapply(y, g, mean)
      s   <- tapply(y, g, sd)
      n1  <- sum(g == levels(g)[1])
      n2  <- sum(g == levels(g)[2])
      sd_p <- sqrt(((n1 - 1)*s[1]^2 + (n2 - 1)*s[2]^2) / (n1 + n2 - 2))
      es   <- abs(diff(m)) / sd_p
    } else {
      # 5b. >2 groups → partial η² from one-way ANOVA
      an <- try(aov(y ~ g), silent = TRUE)
      if (!inherits(an, "try-error")) {
        s    <- summary(an)[[1]]
        ss_g <- s["g", "Sum Sq"]
        ss_e <- s["Residual", "Sum Sq"]
        es   <- ss_g / (ss_g + ss_e)
      }
    }
  } else {
    # 5c. Categorical → Cramer’s V
    orig_lvls    <- levels(x[[1]])
    present_lvls <- unique(as.character(unlist(vals)))
    unused_lvls  <- setdiff(orig_lvls, present_lvls)
    if (length(unused_lvls) > 0) {
      return(
        c("",
          paste0(
            "Warning: unused level(s) in ‘", name, "’: ",
            paste(unused_lvls, collapse = ", "),
            ". Use droplevels(data$", name, ")"
          )
        )
      )
    }
    tbl <- table(y, g)
    ct  <- try(chisq.test(tbl), silent = TRUE)
    if (!inherits(ct, "try-error")) {
      chi <- as.numeric(ct$statistic)
      k   <- min(nrow(tbl), ncol(tbl))
      es  <- sqrt(chi / (n * (k - 1)))
    }
  }

  # 6. ---- If effect size is NA (test failed), return blank ----
  if (is.na(es)) return(c("", ""))

  # 7. ---- Format to three decimals (or "< 0.001") ----
  es_txt <- if (abs(es) < 0.001) "< 0.001" else sprintf("%.3f", es)

  # 8. ---- Return for table1() alignment ----
  c("", es_txt)
}
