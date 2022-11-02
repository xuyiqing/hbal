#' @title Plotting Covariate Balance from an \code{hbal} Object
#' @aliases plot
#' @description This function plots the covariate difference between the control and treatment groups in standardized means before and after weighting.
#' @usage \method{plot}{hbal}(x, type = 'balance', ...)
#' @param x       an object of class \code{hbalobject} as returned by \code{hbal}.
#' @param type    type of graph to plot.
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
plot.hbal <- function(x, 
	type = 'balance',
	...){
	if (type == 'weight'){
		cat('sum(weights) normalized to the number of control units: ', length(x$weights), '\n')
		dat <- data.frame(x = x$weights * length(x$weights))
		ggplot(data=dat, aes_string(x="x", y="..scaled..")) + geom_density() + theme_bw() + labs(y='Density', x='Weights')
	}else{
		plots <- list()
		out.sub <- list()
		width <- list()
		groups <- names(x$grouping)
		var.names <- colnames(x$mat)
		Treatment <- x$Treatment # treatment
		treat <- x$mat[Treatment==1,] # treated
		control <- x$mat[Treatment==0,] * c(x$weights) * sum(Treatment==0) # control
		denom <- apply(rbind(treat, control), 2, sd)
		std.diff.after <- (apply(treat, 2, mean) - apply(control, 2, mean))/denom
		denom <- apply(rbind(treat, x$mat[Treatment==0,]), 2, sd)
		std.diff.before <- (apply(treat, 2, mean) - apply(x$mat[Treatment==0,], 2, mean))/denom
		out <- data.frame(val=c(std.diff.before, std.diff.after),
						  y=rep(var.names, 2),
						  group=factor(rep(c("before adjustment", "after adjustment"), each=length(std.diff.before))),
						  covar.group=rep(rep(groups, x$grouping),2))
		start <- 1
		max_length <- max(str_length(var.names))
		for (i in 1:length(groups)){
			end <- start+x$grouping[i]-1
			out.sub[[i]] <- rbind(out[start:end,], out[(start+length(std.diff.before)):(end+length(std.diff.before)),])
			l <- max(abs(out.sub[[i]]$val), 0.15)
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
	out <- list()
	if (!is.null(object$call)) out[['call']] <- object$call
	
	est <- summary(hbal::att(object, ...))
	
	out[['att']] <- est$coefficients["Treat",]

	out[['group.penalty']] <- object$group.penalty

	out[['term.penalty']] <- object$term.penalty
	
	row.names <- colnames(object$mat)
	balance.tab <- matrix(NA, length(row.names), 3)
	rownames(balance.tab) <- row.names

	Treatment <- object$Treatment # treatment
	treat <- object$mat[Treatment==1,] # treated
	control <- object$mat[Treatment==0,] * c(object$weights) * sum(Treatment==0) # control
	balance.tab[,1] <- apply(treat, 2, mean)
	balance.tab[,2] <- apply(control, 2, mean)
	denom <- apply(rbind(treat, control), 2, sd)
	std.diff.after <- (apply(treat, 2, mean) - apply(control, 2, mean))/denom
	balance.tab[,3] <- std.diff.after
	colnames(balance.tab) <- c("Mean (Treated)", "Mean (Control)", "Std. Mean Diff.")
	balance.tab <- round(balance.tab, 2)
	
	out[['balance']] <- balance.tab

	return(out)
}