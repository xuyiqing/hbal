
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hbal

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![License:
MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**R** package for implementing hierarchically regularized entropy
balancing. It is suitable for estimating average treatment effect on the
treated with binary treatments under strict ignorability. **hbal** is an
extension of entropy balancing: it expands the feature space by
including higher-order terms (such as squared and cubic terms and
interactions) of covariates and then achieves approximate balance on the
expanded features using ridge penalties with a hierarchical structure.

**Authors:** [Yiqing Xu](http://yiqingxu.org/) (Stanford); [Eddie
Yang](https://www.eddieyang.net/) (UCSD)

**Date:** March 27, 2022

**Repo:** [GitHub](https://github.com/xuyiqing/hbal)

**Examples:** R code used in the
[tutorial](https://yiqingxu.org/packages/hbal/articles/tutorial.html)
can be downloaded from [here](hbal_examples.R).

**Reference:** Eddie Yang & Yiqing Xu (2021). [Hierarchically
Regularized Entropy
Balancing](https://papers.ssrn.com/abstract=3807620). *Political
Analysis*, forthcoming.

------------------------------------------------------------------------

## Installation

You can install the **hbal** package from CRAN:

``` r
install.packages('hbal') # Not on CRAN yet.
```

You can also install the up-to-date development version from Github:

``` r
install.packages('devtools', repos = 'http://cran.us.r-project.org') # if not already installed
devtools::install_github('xuyiqing/hbal', ref="main")
```

**hbal** depends on the following packages, which will be installed
automatically when **hbal** is being installed; you can also install
them manually:

``` r
require(estimatr)  
require(glmnet) 
require(ggplot2)
require(gridExtra)
require(gtable)
require(nloptr)
require(Rcpp)
require(RcppEigen)
require(stringr)
```

------------------------------------------------------------------------

### Notes on installation failures

1.  Mac users who have updated to MacOS BigSur or Monterey will likely
    encounter compilation problems. See
    [here](http://yiqingxu.org/public/BigSurError.pdf) for a potential
    solution.
2.  Windows users please consider upgrading R to 4.0.0 or higher and
    installing the [latest
    Rtools](https://cran.r-project.org/bin/windows/Rtools/) to avoid
    C++17 complier errors when installing fastplm.
3.  For Rcpp, RcppArmadillo and MacOS “-lgfortran” and “-lquadmath”
    error, click
    [here](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/)
    for details.
4.  Installation failure related to OpenMP on MacOS, click
    [here](http://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/)
    for a solution.
5.  To fix these issues, try installing gfortran from
    [here](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS%20clang4%20R%20Binaries%20from%20https://github.com/coatless/r-macos-clang).

## Report bugs

Please report bugs to **z5yang \[at\] ucsd.edu** with your sample code
and data file. Much appreciated!
