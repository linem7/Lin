
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Lin

<!-- badges: start -->

[![R-CMD-check](https://github.com/linem7/Lin/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/linem7/Lin/actions/workflows/R-CMD-check.yaml)

<!-- badges: end -->

The goal of Lin is to provide some useful functions

## Installation

``` r
# Install from GitHub (if the repo is public)
install.packages("devtools")
devtools::install_github("linem7/Lin")
```

## Examples

### Packages and working directory

``` r
library(MplusAutomation)
library(knitr)
library(kableExtra)
library(Lin)
```

### Mplus result demostration

``` r
models <- readModels("./inst/extdata") # Using the data from Mplus official site.
```

``` r
# Print model fit for comparison
fit_table(models, ref = "last", indices = c("AIC")) %>% apa(title = "Table 1. Model comparison.")
```

![](./inst/extdata/table.png)

``` r
# Print coefficients
coeff_table(models$ex5.11.out, models$ex5.12.out, models$ex5.13.out)
```

![](./inst/extdata/coeff_1.png)

### Retieve data lables

This is a basic example which shows you how to solve a common problem:

``` r
# Load the package (assuming it's installed or you're in development mode)

library(tibble)

df <- tibble(
  id = 1:3,
  score = c(10, 9, 8)
)
attr(df$id, "label") <- "Participant ID"
attr(df$score, "label") <- "Test Score"

extract_labels(df)
#> # A tibble: 2 × 2
#>   Variable Label         
#>   <chr>    <chr>         
#> 1 id       Participant ID
#> 2 score    Test Score
```
