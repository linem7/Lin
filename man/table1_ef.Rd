% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/table1_utils.R
\name{table1_ef}
\alias{table1_ef}
\title{Effect-Size Calculator for table1()}
\usage{
table1_ef(x, name, ...)
}
\description{
Computes and formats effect sizes for group comparisons in **table1**.
}
\details{
\enumerate{
  \item \strong{Numeric comparisons:}
    \itemize{
      \item Two groups → Cohen’s \emph{d} (pooled SD).
      \item >2 groups → partial \emph{η}\^2 from ANOVA sums of squares.
    }
  \item \strong{Categorical comparisons:}
    \itemize{
      \item Cramer’s \emph{V} from \code{chisq.test()}.
      \item Warns on unused factor levels and suggests \code{droplevels()}.
    }
  \item \strong{Missing & overall:}
    \itemize{
      \item Drops \code{NA}s and empty groups.
      \item Strips any column named “All…” or “Overall…”.
    }
}
}
\note{
If you request \code{overall = "Overall"}, the helper will detect and remove it before computing effect sizes.
}
\examples{
\dontrun{
# Cohen’s d for mpg across 4 vs 6 cylinders
df <- subset(mtcars, cyl \%in\% c(4,6))
table1(
  ~ mpg | cyl,
  data      = df,
  overall   = "Overall",
  extra.col = list("Effect size" = table1_ef)
)
}
}
