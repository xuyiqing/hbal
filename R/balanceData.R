#' @title Balance Statistics from an \code{hbal} Object as a Data Frame
#' @aliases balanceData
#' @description \code{balanceData} returns the covariate balance statistics of an
#'   \code{hbal} object in long form, with a column giving each term's covariate group.
#'   This is the data \code{plot} draws for \code{type = "balance"}, in a shape a
#'   \pkg{ggplot2} call can consume directly.
#' @usage balanceData(hbalobject)
#' @param hbalobject  an object of class \code{hbal} as returned by \code{hbal}.
#' @details The returned frame has two rows per column of \code{hbalobject$mat}: one
#'   before weighting and one after. The values come from \code{hbalobject$bal.tab} and
#'   are therefore rounded to two decimals, exactly as \code{summary} and \code{plot}
#'   report them. \code{covar.group} repeats each group label of
#'   \code{hbalobject$grouping} as many times as that group has columns of \code{mat},
#'   in the order the groups appear.
#' @return A data frame with \eqn{2p} rows, where \eqn{p} is the number of columns of
#'   \code{hbalobject$mat}, and the columns
#'   \describe{
#'     \item{term}{character; the column name of \code{mat} the row refers to.}
#'     \item{covar.group}{factor; the covariate group of that term. The levels are the
#'     names of \code{hbalobject$grouping}, in that order.}
#'     \item{adjustment}{factor with levels \code{"before"} and \code{"after"}, in that
#'     order.}
#'     \item{std.diff}{numeric; the standardized difference in means between the treated
#'     and the control group, before or after weighting.}
#'     \item{tr.mean, co.mean, w.co.mean}{numeric; the treated mean, the unweighted
#'     control mean and the weighted control mean of the term. Each is repeated on both
#'     of the term's rows.}
#'   }
#' @author Yiqing Xu, Eddie Yang
#' @examples
#' #EXAMPLE
#' set.seed(1984)
#' N <- 500
#' X1 <- rnorm(N)
#' X2 <- rbinom(N, size = 1, prob = .5)
#' treat <- rbinom(N, 1, prob = 0.5)
#' Y <- 0.5 * treat + X1 + X2 + rnorm(N)
#' dat <- data.frame(treat = treat, X1 = X1, X2 = X2, Y = Y)
#' out <- hbal(Treat = 'treat', X = c('X1', 'X2'), Y = 'Y', data = dat)
#' balanceData(out)
#' @export

balanceData <- function(hbalobject){
	if (!inherits(hbalobject, "hbal")) {
		stop("hbalobject must be an hbal object from a call to hbal()")
	}
	bal.tab <- hbalobject$bal.tab
	grouping <- hbalobject$grouping
	term <- colnames(hbalobject$mat)
	out <- data.frame(
		term = rep(term, 2),
		covar.group = factor(rep(rep(names(grouping), grouping), 2),
							 levels = names(grouping)),
		adjustment = factor(rep(c("before", "after"), each = length(term)),
							levels = c("before", "after")),
		std.diff = c(bal.tab[, "Std.Diff.(O)"], bal.tab[, "Std.Diff.(W)"]),
		tr.mean = rep(bal.tab[, "Tr.Mean"], 2),
		co.mean = rep(bal.tab[, "Co.Mean"], 2),
		w.co.mean = rep(bal.tab[, "W.Co.Mean"], 2),
		stringsAsFactors = FALSE)
	rownames(out) <- NULL
	out
}
