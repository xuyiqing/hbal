#' Data from Black \& Owens (2016)
#' @description
#' Data on the contender judges from Black \& Owens (2016): Courting the president: how circuit court judges alter their behavior for promotion to the Supreme Court
#' This dataset includes 10,171 period-judge observations for a total of 68 judges.
#' The treatment variable of interest is \code{treatFinal0}, which indicates whether there was a vacancy in the Supreme Court
#' The outcome of interest is ideological alignment of judges' votes with the sitting President (\code{presIdeoVote}).
#' The remaining variables are characteristics of the judges and courts, to be used as controls.
#' @name contenderJudges
#' @docType data
#' 
#' @usage data(contenderJudges)
#' 
#' @format A data frame with 10171 rows and 10 columns.
#' \describe{
#'   \item{presIdeoVote}{ideological alignment of judges' votes with the sitting President (outcome)}
#'   \item{treatFinal0}{treatment indicator for vacancy period}
#'   \item{judgeJCS}{judge’s Judicial Common Space (JCS)score}
#'   \item{presDist}{Ideological distribution of the sitting President}
#'   \item{panelDistJCS}{ideological composition of the panel with whom the judge sat}
#'   \item{circmed}{median JCS score of the circuit judges}
#'   \item{sctmed}{JCS score of the median justice on the Supreme Court}
#'   \item{coarevtc}{indicator for whether the case decision was reversed by the circuit court}
#'   \item{casepub}{indicator for the publication status of thecourt’s opinion}
#'   \item{judge}{name of the judge}
#' }
#' @references 
#' \itemize{
#'   \item Black, R. C., & Owens, R. J. (2016). Courting the president: how circuit court judges alter their behavior for promotion to the Supreme Court. American Journal of Political Science, 60(1), 30-43.
#' }
#' @keywords datasets
NULL