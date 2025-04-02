
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Lin

<!-- badges: start -->
<!-- badges: end -->

The goal of Lin is to provide some useful functions

## Installation

``` r
# Install from GitHub (if the repo is public)
# install.packages("devtools")
# devtools::install_github("linem7/Lin")
```

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# Load the package (assuming it's installed or you're in development mode)
library(Lin)

# Example 1: Use extract_labels on a data frame with labeled columns
data <- data.frame(
  subject = 1:3,
  score   = c(10, 8, 9)
)
attr(data$subject, "label") <- "Subject ID"
attr(data$score, "label") <- "Test Score"

# Extract labels
label_df <- extract_labels(data)
print(label_df)
#> NULL
#> # A tibble: 2 x 2
#>   variable label      
#>   <chr>    <chr>      
#> 1 subject  Subject ID 
#> 2 score    Test Score
```
