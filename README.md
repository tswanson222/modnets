
<!-- README.md is generated from README.Rmd. Please edit that file -->

# modnets

<!-- badges: start -->
<!-- badges: end -->

The goal of modnets is to …

## Installation

You can install the released version of modnets from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("modnets")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(modnets)

data <- na.omit(psychTools::msq[, c('hostile', 'lonely', 'nervous', 'sleepy', 'depressed')])

fit <- fitNetwork(data, moderators = 'depressed')

plot(fit, threshold = TRUE)
```

<img src="man/figures/README-example-1.png" width="100%" />
