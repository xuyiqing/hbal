#' @title Match Column Names to be Excluded
#' @aliases covarExclude
#' @description Internal function called by \code{hbal} to serially expand covariates.
#' @param colname              column name.
#' @param exclude              list of covariate name pairs or triplets to be excluded.
#' @return Logical
#' @author Yiqing Xu, Eddie Yang

covarExclude <- function(colname, exclude){
	out <- FALSE
	for (i in 1:length(exclude)){
		drop.name <- exclude[[i]]
		drop.logical <- rep(NA, length(drop.name))
		for (k in 1:length(drop.name)){
			drop.logical[k] <- grepl(drop.name[k], colname)
		}
		if (sum(drop.logical)==length(drop.logical)){
			out <- TRUE
		}
	}
	return(out)
}

#covarExclude <- function(colname, exclude){
#	out <- FALSE
#	for (i in 1:length(exclude)){
#		drop.element <- combinat::permn(exclude[[i]])
#		drop.name <- sapply(drop.element, paste0, sep="", collapse=".")	
#		if (colname%in%drop.name){
#			out <- TRUE
#		}
#	}
#	return(out)
#}