#' @title Serial Expansion of Covariates
#' @aliases covarExpand
#' @description Internal function called by \code{hbal} to serially expand covariates.
#' @param X                    matrix of covariates.
#' @param exp.degree           the degree of the polynomial. 
#' @param treatment            treatment indicator
#' @param exclude              list of covariate name pairs or triplets to be excluded.
#' @return A matrix of serially expanded covariates
#' @author Yiqing Xu, Eddie Yang
#' @importFrom stats poly
#' @importFrom stringr str_trunc
#' @export

covarExpand <- function(X, exp.degree=3, treatment=NULL, exclude=NULL){

	f1 <- function(a, b){
    return(sum(a==b))
	}

	if (exp.degree!=2 && exp.degree!=3){
		stop("covarExpand currently onyl supports series expansion to the second or thrid degree")
	}

	X.P <- poly(X, degree = exp.degree, raw = TRUE)
	degree <- attributes(X.P)$degree
	acending.degree <- order(degree)
	col.names <- colnames(X.P)

	if (exp.degree==2){
		second.start <- sum(degree==1)+1
		term.order <- list(
			linear = 1:sum(degree==1),
			second = second.start:length(degree)
			)

		linear.terms <- acending.degree[term.order[["linear"]]]
		sq.terms <- acending.degree[term.order[["second"]]]

		sq.terms.sq <- sq.terms[grepl("2", col.names[sq.terms])]
		sq.terms.interact <- sq.terms[!grepl("2", col.names[sq.terms])]
		pos <- list(linear.terms, sq.terms.sq, sq.terms.interact)

		X.P <- X.P[,unlist(pos)]
		rank <- qr(X.P[treatment==0,])
		X.P <- X.P[,rank$pivot[1:rank$rank]]
		split.name <- strsplit(colnames(X.P), ".", fixed=TRUE)
		new.name <- colnames(X)
		new.name <- str_trunc(new.name, 5, side="right", ellipsis="")
		for (nam in 1:ncol(X.P)){
			sp <- as.numeric(split.name[[nam]])
			colnames(X.P)[nam] <- paste0(rep(new.name[sp!=0], sp[sp!=0]), collapse = ".")
		}

		if (!is.null(exclude)){
			to.drop <- sapply(colnames(X.P), covarExclude, exclude=exclude)
			X.P <- X.P[,!to.drop]
		} else{
			to.drop <- rep(FALSE, ncol(X.P))
		}

		order.seq <- rep(1:3, times=unlist(lapply(pos, length)))
		order.seq <- order.seq[rank$pivot[1:rank$rank]]
		order.seq <- order.seq[!to.drop]
		grouping <- sapply(1:3, f1, order.seq)
		out <- list(mat=X.P, grouping=grouping)
	}

	if (exp.degree==3){
		# starting position of each order of terms
		second.start <- sum(degree==1)+1
		third.start <- sum(degree==1)+sum(degree==2)+1
		term.order <- list(
			linear = 1:sum(degree==1),
			second = second.start:(third.start-1),
			third = (third.start):length(degree)
			)

		linear.terms <- acending.degree[term.order[["linear"]]]
		sq.terms <- acending.degree[term.order[["second"]]]
		cube.terms <- acending.degree[term.order[["third"]]]

		sq.terms.sq <- sq.terms[grepl("2", col.names[sq.terms])]
		sq.terms.interact <- sq.terms[!grepl("2", col.names[sq.terms])]
		cube.terms.cube <- cube.terms[grepl("3", col.names[cube.terms])]
		cube.terms.interact.sq <- cube.terms[grepl("2", col.names[cube.terms])]	
		cube.terms.interact <- cube.terms[!cube.terms%in%c(cube.terms.cube, cube.terms.interact.sq)]
		pos <- list(linear.terms, sq.terms.interact, sq.terms.sq,
			cube.terms.interact, cube.terms.interact.sq, cube.terms.cube)

		X.P <- X.P[,unlist(pos)]
		rank <- qr(X.P[treatment==0,])
		X.P <- X.P[,rank$pivot[1:rank$rank]]
		split.name <- strsplit(colnames(X.P), ".", fixed=TRUE)
		new.name <- colnames(X)
		new.name <- str_trunc(new.name, 5, side="right", ellipsis="")
		for (nam in 1:ncol(X.P)){
			sp <- as.numeric(split.name[[nam]])
			colnames(X.P)[nam] <- paste0(rep(new.name[sp!=0], sp[sp!=0]), collapse = ".")
		}

		if (!is.null(exclude)){
			to.drop <- sapply(colnames(X.P), covarExclude, exclude=exclude)
			X.P <- X.P[,!to.drop]
		} else{
			to.drop <- rep(FALSE, ncol(X.P))
		}

		order.seq <- rep(1:6, times=unlist(lapply(pos, length)))
		order.seq <- order.seq[rank$pivot[1:rank$rank]]
		order.seq <- order.seq[!to.drop]
		grouping <- sapply(1:6, f1, order.seq)
		out <- list(mat=X.P, grouping=grouping)
	}
	return(out)
}