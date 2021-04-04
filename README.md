
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hbal

<!-- badges: start -->

[![Travis build status](https://travis-ci.com/xuyiqing/hbal.svg?branch=main)](https://travis-ci.com/xuyiqing/hbal)
<!-- badges: end -->

hbal performs hierarchically regularized entropy balancing such that the
covariate moments of the control group match those of the treatment
group. hbal automatically expands the covariate space to include higher
order terms and uses cross-validation to select variable penalties for
the balancing conditions.

## Installation

You can install the released version of hbal from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("hbal")
```

Or the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("xuyiqing/hbal", ref="main")
```

## Getting started

``` r
library(hbal)
```

hbal provides two main functions:

  - `hbal()`, which performs hierarchically regularized entropy
    balancing.

  - `att()`, which calculates the average treatment effect on the
    treated (ATT) from an `hbalobject` returned by `hbal()`.

## Example

For examples on how to use the package, see
[here](https://yiqingxu.org/software/hbal/hbal.html)

## Reference

For a detailed description of the method see: [Xu and Yang (2021)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3807620)
