% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_mplus_inp.R
\name{create_mplus_inp}
\alias{create_mplus_inp}
\title{Create Mplus .inp File}
\usage{
create_mplus_inp(
  syntax_string,
  filename,
  target_dir = ".",
  overwrite = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{syntax_string}{Character. A single string containing newlines or a character vector of lines defining the Mplus syntax.}

\item{filename}{Character. Desired name of the .inp file. Can include or omit the `.inp` extension.}

\item{target_dir}{Character. Path to the directory where the file should be written. Will be created if it doesn't exist.}

\item{overwrite}{Logical. Should an existing file be overwritten? Defaults to `FALSE`.}

\item{verbose}{Logical. Should a confirmation message be printed? Defaults to `TRUE`.}
}
\value{
Invisibly returns the full path to the created .inp file.
}
\description{
Write a vector of Mplus syntax lines or a single string (with newlines) to an .inp file.
}
\examples{
# Single string example with embedded newlines
syntax <- "TITLE: CFA model;
DATA: FILE = 'data.dat';
VARIABLE: NAMES = x1-x5;
MODEL: f1 BY x1-x3;
OUTPUT: STANDARDIZED;"
# Write to ./Mplus/cfa.inp, overwriting if exists, and print message
create_mplus_inp(syntax_string = syntax, filename = "cfa", "./Mplus", overwrite = TRUE)
}
