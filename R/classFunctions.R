#' @title Plotting Covariate Balance from an \code{hbal} Object
#' @aliases plot
#' @description This function plots the covariate difference between the control and treatment groups in standardized means before and after weighting.
#' @usage \method{plot}{hbal}(x, ...)
#' @param x       an object of class \code{hbalobject} as returned by \code{hbal}.
#' @param ...     Further arguments to be passed to \code{plot.hbal()}.
#' @return A matrix of ggplots of covariate balance by group
#' @import ggplot2
#' @import gridExtra
#' @importFrom stringr str_length str_pad
#' @importFrom gtable gtable_filter
#' @importFrom stats sd
#' @author Yiqing Xu, Eddie Yang
#' @export
#' @method plot hbal
plot.hbal <- function(x, ...){
	plots <- list()
	out.sub <- list()
	width <- list()
	groups <- c("linear terms", "two-way interactions", "square terms", "three-way interactions", "linear-squre interactions", "cubic terms")
	groups <- groups[1:length(x$group.assignment)]
	if(length(which(x$group.assignment==0))!=0){
		groups <- groups[-which(x$group.assignment==0)]
		x$group.assignment <- x$group.assignment[-which(x$group.assignment==0)]
	}
	var.names <- colnames(x$mat)
	T <- x$Treatment # treatment
	treat <- x$mat[T==1,] # treated
	control <- x$mat[T==0,] * c(x$weights) * sum(T==0) # control
	denom <- apply(rbind(treat, control), 2, sd)
	std.diff.after <- (apply(treat, 2, mean) - apply(control, 2, mean))/denom
	denom <- apply(rbind(treat, x$mat[T==0,]), 2, sd)
	std.diff.before <- (apply(treat, 2, mean) - apply(x$mat[T==0,], 2, mean))/denom
	out <- data.frame(val=c(std.diff.before, std.diff.after),
					  y=rep(var.names, 2),
					  group=factor(rep(c("before adjustment", "after adjustment"), each=length(std.diff.before))),
					  covar.group=rep(rep(groups, x$group.assignment),2))
	start <- 1
	max_length <- max(str_length(var.names))
	for (i in 1:length(groups)){
		end <- start+x$group.assignment[i]-1
		out.sub[[i]] <- rbind(out[start:end,], out[(start+length(std.diff.before)):(end+length(std.diff.before)),])
		l <- max(abs(out.sub[[i]]$val))
		plots[[i]] <- ggplot(aes_string(x="val", y="y"), data=out.sub[[i]]) + geom_point(size=3, shape = 21, colour = "black", aes_string(fill="group")) + 
						scale_y_discrete(limits=rev(var.names[start:end]), labels=str_pad(rev(var.names[start:end]),max_length, side = "left")) + xlim(-l, l) + scale_fill_manual(values=c("black", "white")) +
						geom_vline(xintercept = -0.1, lty=2) + geom_vline(xintercept = 0.1, lty=2) + theme_bw() + labs(x="Standardized Difference", y="") + 
						theme(legend.title = element_blank()) + ggtitle(groups[i]) + theme(axis.text.y = element_text(family = "mono")) +theme(legend.position="bottom")
		
		start <- end+1
	}
	legend <- gtable_filter(ggplot_gtable(ggplot_build(plots[[1]])), "guide-box")
	for (i in 1:length(groups)){
		plots[[i]] <- plots[[i]] + theme(legend.position="none")
	}
	grid.arrange(do.call("arrangeGrob", c(plots, ncol=ceiling(length(plots)/2))),  nrow=2, legend, heights=c(1.1, 0.1))

}


#' @title Summarizing from an \code{hbal} Object
#' @aliases summary
#' @description This function prints a summary from an \code{hbal} Object.
#' @usage \method{summary}{hbal}(object, ...)
#' @param object  an object of class \code{hbalobject} as returned by \code{hbal}.
#' @param ...     Further arguments to be passed to \code{summary.hbal()}.
#' @return a summary table
#' @importFrom stats sd
#' @author Yiqing Xu, Eddie Yang
#' @export
#' @method summary hbal
summary.hbal <- function(object, ...){
	if (!is.null(object$call)) cat("\nCall:", deparse(object$call), sep = "\n")
	groups <- c("linear", "two-way interact", "square", "three-way interact", "linear-squre interact", "cubic")
	groups <- groups[1:length(object$group.assignment)]
	if(length(which(object$group.assignment==0))!=0){
		groups <- groups[-which(object$group.assignment==0)]
		object$group.assignment <- object$group.assignment[-which(object$group.assignment==0)]
	}

	att <- summary(hbal::att(object))
	cat("\nAverage Treatment Effect on the Treated:\n")
	print(round(att$coefficients["Treat",], 3))
	
	cat("\nCovariates included in adjustment:\n")
	start <- 1
	for (i in 1:length(groups)){
		end <- start+object$group.assignment[i]-1
		msg <- paste0(colnames(object$mat)[start:end], collapse=", ")
		cat(groups[i], ":\n", msg, "\n\n")
		start <- end+1
	}	
	balance.out <- list()
	row.names <- colnames(object$mat)
	balance.tab <- matrix(NA, length(row.names), 3)
	rownames(balance.tab) <- row.names

	if(!is.null(object$penalty)){
		cat("penalty for each group of covariates selected:\n")
		penalty <- matrix(NA, 1, length(groups))
		colnames(penalty) <- groups
		rownames(penalty) <- "penalty value"
		penalty[1,] <- object$penalty
		print(round(penalty, 2))
	}
	cat("\nSummary of Balance for Matched Data:\n")
	T <- object$Treatment # treatment
	treat <- object$mat[T==1,] # treated
	control <- object$mat[T==0,] * c(object$weights) * sum(T==0) # control
	balance.tab[,1] <- apply(treat, 2, mean)
	balance.tab[,2] <- apply(control, 2, mean)
	denom <- apply(rbind(treat, control), 2, sd)
	std.diff.after <- (apply(treat, 2, mean) - apply(control, 2, mean))/denom
	balance.tab[,3] <- std.diff.after
	colnames(balance.tab) <- c("Means, Treated", "Means, Control", "Std. Mean Diff.")
	balance.tab <- round(balance.tab, 2)
	start <- 1
	for (i in 1:length(object$group.assignment)){
		end <- start+object$group.assignment[i]-1
		balance.out[[groups[i]]] <- balance.tab[start:end,]
		start <- end+1
	}
	balance.out

}