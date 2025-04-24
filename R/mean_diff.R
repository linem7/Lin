#' Mean-Difference Table (APA Style)
#'
#' \code{mean_diff()} evaluates mean differences for one or more dependent variables
#' \code{dv} across one or more grouping variables \code{iv}.  If a grouping variable
#' has exactly two levels the function runs an independent‐samples \strong{t test};
#' otherwise it runs a \strong{one-way ANOVA}.  Results are returned as a named list
#' of data frames—one data frame per dependent variable—ready to paste into a
#' manuscript or report.
#'
#' Each data-frame contains one row per category, displaying M (SD), the test statistic, the exact p value, an optional effect size,
#' and (for significant ANOVAs) Bonferroni or Holm post-hoc contrasts written
#' as “A > B”.  To avoid repetition, the IV name appears only in the first row
#' of its block.
#'
#' @param dv Character vector of dependent-variable names.
#' @param iv Character vector of independent-variable (grouping) names.
#' @param data Data frame that contains all specified `dv` and `iv` columns.
#' @param stars Logical.  Add significance symbols to the test statistic,
#' (`*p < .05, ** p < .01, *** p < .001`).  Default is `FALSE`.
#' @param effect Logical.  Add effect sizes (Cohen’s d for t tests,
#'   partial η² for ANOVAs).  Default is `FALSE`.
#' @param digits_desc Integer.  Number of decimal places for means and SDs.
#'   Default is 2.
#' @param digits_stat Integer.  Decimal places for test statistics, p values,
#'   and effect sizes.  Default is 3.
#' @param p.adjust.method Character.  p-value adjustment for post-hoc tests;
#'   one of "bonferroni (default)", "holm", or "none".
#'
#' @return A named list of data frames.  Each data frame includes:
#' \describe{
#'   \item{IV}{Grouping variable (blank after the first row).}
#'   \item{Category}{Factor level.}
#'   \item{N}{Sample size.}
#'   \item{Mean (SD)}{Descriptive statistic.}
#'   \item{t/F}{Test statistic (with stars if \code{stars = TRUE}).}
#'   \item{P}{Exact p value (“< 0.001” when appropriate).}
#'   \item{Effect Size}{Cohen’s d or partial η² (absent if \code{effect = FALSE}).}
#'   \item{PostHoc}{Significant pairwise contrasts for ANOVA.}
#' }
#'
#' @importFrom stats aggregate t.test aov complete.cases
#' @importFrom car leveneTest
#' @importFrom effsize cohen.d
#'
#' @examples
#' ## Example using the mock dataset `test_data`
#' data(test_data)
#'
#' out <- mean_diff(
#'   dv = c("var1", "var2"),
#'   iv = c("gender", "grade", "area"),
#'   data  = test_data,
#'   stars = TRUE,
#'   effect = TRUE,
#'   p.adjust.method = "holm"
#' )
#'
#' out$var1    # results for var1
#' out$var2    # results for var2
#'
#' @export
mean_diff <- function(
    dv, iv, data,
    stars = FALSE, effect = FALSE,
    digits_desc = 2, digits_stat = 3,
    p.adjust.method = c("bonferroni", "holm", "none")
) {
  p.adjust.method <- match.arg(p.adjust.method)
  stopifnot(is.data.frame(data))
  stopifnot(all(dv %in% names(data)), all(iv %in% names(data)))

  num_fmt <- function(x, d) sprintf(paste0("%0.", d, "f"), x)
  star_fun <- function(p) {
    if (p < .001) "***" else if (p < .01) "**" else if (p < .05) "*" else ""
  }
  p_fmt <- function(p) if (p < .001) "< 0.001" else num_fmt(p, digits_stat)

  out <- setNames(vector("list", length(dv)), dv)

  for (d in dv) {
    rows_list <- list()

    for (g in iv) {
      sub <- data[, c(d, g)]
      sub <- sub[complete.cases(sub), ]
      sub[[g]] <- factor(sub[[g]])
      k <- nlevels(sub[[g]])
      if (k < 2) {
        warning(sprintf("Skipping '%s' by '%s': fewer than 2 groups.", d, g))
        next
      }

      # Descriptives ---------------------------------------------------------
      desc <- aggregate(sub[[d]], list(Category = sub[[g]]), function(x) {
        c(N = length(x), M = mean(x), SD = sd(x))
      })
      desc <- do.call(data.frame, desc)
      names(desc) <- c("Category", "N", "Mean", "SD")
      desc$`Mean (SD)` <- sprintf(paste0("%0.", digits_desc, "f (%0.", digits_desc, "f)"),
                                  desc$Mean, desc$SD)
      desc$Mean <- desc$SD <- NULL

      # Tests ---------------------------------------------------------------
      if (k == 2) {
        lev_p <- car::leveneTest(sub[[d]] ~ sub[[g]])[1, "Pr(>F)"]
        var_equal <- lev_p > .05
        t_res <- t.test(sub[[d]] ~ sub[[g]], var.equal = var_equal)
        stat_val <- t_res$statistic
        p_val   <- t_res$p.value
        effect_val <- if (effect) effsize::cohen.d(sub[[d]] ~ sub[[g]], hedges.correction = FALSE)$estimate else NA
        post_str <- ""
      } else {
        aov_res <- aov(sub[[d]] ~ sub[[g]], data = sub)
        tab <- summary(aov_res)[[1]]
        stat_val <- tab[1, "F value"]
        p_val   <- tab[1, "Pr(>F)"]
        ss_eff  <- tab[1, "Sum Sq"]; ss_err <- summary(aov_res)[[1]][2, "Sum Sq"]
        effect_val <- if (effect) ss_eff / (ss_eff + ss_err) else NA
        post_str <- ""
        if (p_val < .05) {
          pw <- pairwise.t.test(sub[[d]], sub[[g]], p.adjust.method = p.adjust.method)
          sig <- which(pw$p.value < .05, arr.ind = TRUE)
          if (nrow(sig)) {
            comps <- apply(sig, 1, function(idx) {
              l1 <- rownames(pw$p.value)[idx[1]]; l2 <- colnames(pw$p.value)[idx[2]]
              m1 <- mean(sub[[d]][sub[[g]]==l1]); m2 <- mean(sub[[d]][sub[[g]]==l2])
              if (m1 > m2) paste0(l1, " > ", l2) else paste0(l2, " > ", l1)
            })
            post_str <- paste(unique(comps), collapse = "; ")
          }
        }
      }

      stat_str <- paste0(num_fmt(stat_val, digits_stat), if (stars) star_fun(p_val) else "")
      p_str    <- p_fmt(p_val)
      eff_str  <- if (effect) num_fmt(effect_val, digits_stat) else ""

      # Assemble rows -------------------------------------------------------
      desc$IV <- g
      desc$`t/F` <- ""; desc$P <- ""; desc$`Effect Size` <- ""; desc$PostHoc <- ""
      desc$`t/F`[1]      <- stat_str
      desc$P[1]          <- p_str
      desc$`Effect Size`[1] <- eff_str
      desc$PostHoc[1]    <- post_str
      desc$IV[-1] <- ""   # blank repeated IV entries

      # Drop Effect Size column if not requested
      if (!effect) desc$`Effect Size` <- NULL

      desc <- desc[, c("IV","Category","N","Mean (SD)","t/F","P",
                       if (effect) "Effect Size" else NULL,
                       "PostHoc")]
      rows_list[[g]] <- desc
    }

    if (length(rows_list)) {
      out[[d]] <- do.call(rbind, rows_list)
      row.names(out[[d]]) <- NULL
    }
  }
  return(out)
}
