% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_cfa_model.R
\name{create_cfa_model}
\alias{create_cfa_model}
\title{Create lavaan model}
\usage{
create_cfa_model(lav_name, pattern, var_index)
}
\arguments{
\item{lav_name}{names of latent variable}

\item{pattern}{pattern of the model structure}

\item{var_index}{which item should be counted into the model}
}
\value{
a string contains the model structure
}
\description{
A more easy way to create lavaan style cfa model syntax
}
\examples{

create_cfa_model("pasu_f", "pasu{i}_f", 5)

# For multiple factors, combine them with paste after place into a vector:
paste(c(create_cfa_model("pasu_f", "pasu{i}_f", 5),
create_cfa_model("cmu_f", "cmu{i}_f", 3),
create_cfa_model("att_f", "att{i}_f", 4)), collapse = "\n")

}
