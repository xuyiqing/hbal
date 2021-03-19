#' hbal
#' 
#' hbal performs hierarchically regularized entropy balancing such that the covariate distributions of the control group match those of the treatment group. hbal automatically expands the covariate space to include higher order terms and uses cross-validation to select variable penalties for the balancing conditions.
#' 
#' @docType package
#' @author Yiqing Xu <yiqingxu@stanford.edu>, Eddie Yang <z5yang@ucsd.edu>
#' @importFrom Rcpp evalCpp
#' @useDynLib hbal, .registration = TRUE
#' @name hbal
NULL