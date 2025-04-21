
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
library(bruceR)
library(MplusAutomation)
library(knitr)
library(Lin)

set.wd()
#> ✔ Set working directory to "D:/Package/Lin"
```

### Mplus result demostration

``` r
models <- readModels("./inst/extdata") # Using the data from Mplus offcial site.
```

``` r
# Print model fit for comparison
fit_table(models, ref = "last", indices = c("AIC")) %>% apa(title = "Table 1. Model comparison.")
```

<table class=" lightable-classic" style="font-family: Cambria; ">
<caption>
Table 1. Model comparison.
</caption>
<thead>
<tr>
<th style="text-align:left;">
Models
</th>
<th style="text-align:left;">
Chisq
</th>
<th style="text-align:right;">
df
</th>
<th style="text-align:left;">
chisq/df
</th>
<th style="text-align:left;">
CFI
</th>
<th style="text-align:left;">
TLI
</th>
<th style="text-align:left;">
RMSEA
</th>
<th style="text-align:left;">
SRMR
</th>
<th style="text-align:left;">
AIC
</th>
<th style="text-align:left;">
ΔChisq
</th>
<th style="text-align:left;">
p_ΔChisq
</th>
<th style="text-align:left;">
ΔCFI
</th>
<th style="text-align:left;">
ΔRMSEA
</th>
<th style="text-align:left;">
ΔAIC
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
ex5.11.out
</td>
<td style="text-align:left;">
53.704
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:left;">
1.074
</td>
<td style="text-align:left;">
0.997
</td>
<td style="text-align:left;">
0.997
</td>
<td style="text-align:left;">
0.012
</td>
<td style="text-align:left;">
0.027
</td>
<td style="text-align:left;">
19373.920
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
</tr>
<tr>
<td style="text-align:left;">
ex5.12.out
</td>
<td style="text-align:left;">
53.704
</td>
<td style="text-align:right;">
50
</td>
<td style="text-align:left;">
1.074
</td>
<td style="text-align:left;">
0.997
</td>
<td style="text-align:left;">
0.997
</td>
<td style="text-align:left;">
0.012
</td>
<td style="text-align:left;">
0.027
</td>
<td style="text-align:left;">
19373.920
</td>
<td style="text-align:left;">
0.000 (0)
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
0.000
</td>
<td style="text-align:left;">
0.000
</td>
<td style="text-align:left;">
0.000
</td>
</tr>
<tr>
<td style="text-align:left;">
ex5.13.out
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:right;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
19426.131
</td>
<td style="text-align:left;">
-50.210 (1)
</td>
<td style="text-align:left;">
1.000
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
52.211
</td>
</tr>
<tr>
<td style="text-align:left;">
ex5.20.out
</td>
<td style="text-align:left;">
4.090
</td>
<td style="text-align:right;">
10
</td>
<td style="text-align:left;">
0.409
</td>
<td style="text-align:left;">
1.000
</td>
<td style="text-align:left;">
1.000
</td>
<td style="text-align:left;">
0.000
</td>
<td style="text-align:left;">
0.012
</td>
<td style="text-align:left;">
7831.898
</td>
<td style="text-align:left;">
11546.232 (-24)
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
NA
</td>
<td style="text-align:left;">
-11594.233
</td>
</tr>
</tbody>
</table>

``` r
# Print coefficients
# coeff_table(models, model_names = paste0("Models ", 1:4))
```

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
