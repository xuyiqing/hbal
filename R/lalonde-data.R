#' Data from Hazlett (2020)
#' @description
#' Data on the treated units is from Dehejia and Wahba (1999), containing 185 individuals; data on the control units is from Panel Study
#' of Income Dynamics (PSID-1), containing 2,490 individuals.
#'
#' @name lalonde
#' @docType data
#' 
#' @format A data frame with 2675 rows and 13 columns.
#' \describe{
#'   \item{nsw}{treatment indicator of whether an individual participated in the National Supported Work (NSW) program}
#'   \item{age}{}
#'   \item{educ}{years of education}
#'   \item{black}{demographic indicator variables for Black}
#'   \item{hisp}{idemographic indicator variables for Hispanic}
#'   \item{married}{demographic indicator variables for married}
#'   \item{re74}{real earnings in 1974}
#'   \item{re75}{real earnings in 1975}
#'   \item{re78}{real earnings in 1978, outcome}
#'   \item{u74}{unemployment indicator for 1974}
#'   \item{u75}{unemployment indicator for 1975}
#'   \item{u78}{unemployment indicator for 1978}
#'   \item{nodegr}{indicator for no high school degree}
#' }
#' @source
#' National Supported Work (NSW) demonstration data from LaLonde (1986), with
#' the treated units as re-analyzed by Dehejia and Wahba (1999) and PSID-1
#' controls. Taken from the replication archive of Xu and Yang (2022): Xu, Y.,
#' and Yang, E. Replication Data for: Hierarchically Regularized Entropy
#' Balancing. Harvard Dataverse, \doi{10.7910/DVN/QI2WP9}. That archive is
#' released under the Creative Commons CC0 1.0 Universal Public Domain
#' Dedication.
#' @references
#' \itemize{
#'   \item Dehejia, R. H., and Wahba, S. (1999). Causal effects in nonexperimental studies: Reevaluating the evaluation of training programs. Journal of the American statistical Association, 94(448), 1053-1062.
#'   \item Hazlett, C. (2020). KERNEL BALANCING. Statistica Sinica, 30(3), 1155-1189.
#' }
#' @keywords datasets
NULL