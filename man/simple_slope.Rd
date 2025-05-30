% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simple_slope.R
\name{simple_slope}
\alias{simple_slope}
\title{Simple Slope Plot Based on MplusAutomation Output}
\usage{
simple_slope(
  model,
  x,
  w,
  int,
  y,
  slp_l,
  slp_h,
  mean_x = 0,
  mean_w = 0,
  legend_position = "right",
  legend_name = "Moderator Level",
  fontsize = 14
)
}
\arguments{
\item{model}{An MplusAutomation model object containing parameters (typically output from `readModels()`).}

\item{x}{Name of the independent variable in Mplus output.}

\item{w}{Name of the moderator variable.}

\item{int}{Name of the interaction term (x*w).}

\item{y}{Name of the dependent variable.}

\item{slp_l}{Name of the parameter for the simple slope at low moderator value under model constraint (M - 1SD).}

\item{slp_h}{Name of the parameter for the simple slope at high moderator value  under model constraint (M + 1SD).}

\item{mean_x}{Mean value of the independent variable (default = 0).}

\item{mean_w}{Mean value of the moderator (default = 0).}

\item{legend_position}{Position of the legend. Options are:
\itemize{
  \item `"UL"` = Upper Left inside plot,
  \item `"UR"` = Upper Right inside plot,
  \item `"BL"` = Bottom Left inside plot,
  \item `"BR"` = Bottom Right inside plot,
  \item `"right"` = Default, outside plot on right side.
}}

\item{legend_name}{The name of legend.}

\item{fontsize}{Font size for the plot text (default = 14).}
}
\value{
A ggplot2 object displaying the simple slope plot.
}
\description{
This function extracts model parameters from an MplusAutomation object and creates
a simple slope plot for a moderated regression model:
}
\details{
\deqn{y = b0 + b1*x + b2*w + b3*(x*w) + \varepsilon}

It automatically retrieves regression coefficients, variances, and simple slopes
at different moderator levels (mean, mean - 1SD, mean + 1SD), and generates a plot
with annotated legend entries showing the slope estimates, standard errors, and p-values.


The function extracts:
\itemize{
  \item The intercept from the \code{Intercepts} section.
  \item The coefficients for the independent variable (x), moderator (w), and interaction (x*w) from the \code{y.ON} section.
  \item The variances of x and w from the \code{Variances} section.
  \item The simple slopes at M - 1SD and M + 1SD from the \code{New.Additional.Parameters} section.
}

The legend labels dynamically include the simple slope estimates, standard errors, and p-values, separated by a new line (`\\n`).
}
\examples{
\dontrun{
# Example usage:
simple_slope(model,
             x = "x",
             w = "w",
             int = "int",
             y = "y",
             slp_h = "slp_hi",
             slp_l = "slp_lo",
             legend_position = "UL")
}

}
