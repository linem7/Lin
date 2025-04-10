% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/apa_corr_matrix.R
\name{apa_corr_matrix}
\alias{apa_corr_matrix}
\title{Correlation table}
\usage{
apa_corr_matrix(
  data,
  vars,
  digits_stat = 2,
  digits_corr = 3,
  stats_position = c("columns", "rows")
)
}
\arguments{
\item{data}{A data.frame}

\item{vars}{variables included in correlation analysis, use \code{select} from dplyr to choose columns}

\item{digits_stat}{digits in the mean and sd analysis}

\item{digits_corr}{digits in the correlation coefficients}

\item{stats_position}{should be mean and standard deviation placed in the first two columns or the bottom of the table}
}
\value{
a data.frame
}
\description{
Function for creating correlation matrix including Mean and SD.
}
\examples{
cor_data <- select(data, gender, age, ses, var1, var2, var3)

}
