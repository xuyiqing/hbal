#' @title Plotting Covariate Balance from an \code{hbal} Object
#' @aliases plot
#' @description This function plots the covariate difference between the control and treatment groups in standardized means before and after weighting.
#' @usage \method{plot}{hbal}(x, type = 'balance', ...)
#' @param x       an object of class \code{hbalobject} as returned by \code{hbal}.
#' @param type    type of graph to plot.
#' @param log     log scale for the weight plot
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
	log = TRUE,
	base_size = 10,
	...){
	if (type == 'weight'){
		cat('sum(weights) normalized to the number of treated units\n')
		w <- x$weights.co
		if (log  == TRUE) {
			dat <- data.frame(x = log(w)); xlab <- 'Weights (log10)'					
		} else {
			dat <- data.frame(x = w); xlab <- 'Weights'			
		}
		ggplot(data=dat, aes_string(x="x")) + geom_histogram(aes(y=after_stat(density)), color="black", fill="white", bins = 50) + 
			labs(y='Density', x=xlab) + geom_density(alpha=.2, fill="#FF6666") + theme_classic() 			
	}else{
		plots <- list()
		out.sub <- list()
		width <- list()
		groups <- names(x$grouping)
		var.names <- colnames(x$mat) # remove outcome & treat
		std.diff.before <- x$bal.tab[,4]
		std.diff.after <- x$bal.tab[,5]
		out <- data.frame(val=c(std.diff.before, std.diff.after),
						  y=rep(var.names, 2),
						  group=factor(rep(c("before adjustment", "after adjustment"), each=length(std.diff.before)), levels = c("before adjustment", "after adjustment")),
						  covar.group=rep(rep(groups, x$grouping),2))
		start <- 1
		max_length <- max(str_length(var.names))
		for (i in 1:length(groups)){
			end <- start+x$grouping[i]-1
			out.sub[[i]] <- rbind(out[start:end,], out[(start+length(std.diff.before)):(end+length(std.diff.before)),])
			l <- max(abs(out.sub[[i]]$val), 0.15)
			plots[[i]] <- ggplot(aes_string(x="val", y="y"), data=out.sub[[i]]) + geom_point(size=3, shape = 21, colour = "black", aes_string(fill="group")) +
							scale_fill_manual(values = c("white","black")) + 
							scale_y_discrete(limits=rev(var.names[start:end]), labels=str_pad(rev(var.names[start:end]),max_length, side = "left")) + xlim(-l, l) + 
							geom_vline(xintercept = -0.1, lty=2) + geom_vline(xintercept = 0.1, lty=2) + theme_bw(base_size = base_size) + labs(x="Std. Diff.", y="") + 
							theme(legend.title = element_blank(), legend.position="bottom") + ggtitle(groups[i]) + theme(axis.text.y = element_text(family = "mono")) 			
			start <- end+1
		}
		legend <- gtable_filter(ggplot_gtable(ggplot_build(plots[[1]])), "guide-box")
		for (i in 1:length(groups)){
			plots[[i]] <- plots[[i]] + theme(legend.position="none")
		}
		grid.arrange(do.call("arrangeGrob", c(plots, ncol=ceiling(length(plots)/2))),  nrow=2, legend, heights=c(1.1, 0.1))
	}
}

########################################################
########################################################
########################################################
########################################################

#' @title Summarizing from an \code{hbal} Object
#' @aliases summary
#' @description This function prints a summary from an \code{hbal} Object.
#' @usage \method{summary}{hbal}(object, ...)
#' @param object  an object of class \code{hbalobject} as returned by \code{hbal}.
#' @param print.level  level of detials to be printed
#' @param ...     Further arguments to be passed to \code{summary.hbal()}.
#' @return a summary table
#' @importFrom stats sd
#' @author Yiqing Xu, Eddie Yang
#' @export
#' @method summary hbal
summary.hbal <- function(object, print.level = 0, ...){
	
	# out <- list()
	# if (!is.null(object$call)) out[['call']] <- object$call

	# print(object$call)
	# cat("\n")

	# #est <- summary(hbal::att(object, ...))
	#out[['att']] <- est$coefficients["Treat",]
	cat("Call:\n ")
	print(object$call)
	cat("\n")

	# out[['group.penalty']] <- object$group.penalty
	# out[['term.penalty']] <- object$term.penalty
	Treatment <- object$Treatment # treatment
	nTr <- sum(Treatment == 1)
	nCo <- sum(Treatment == 0)

	# no. of treated and controls
	tr <- rev(table(Treatment))
	names(tr) <- c("Treated", "Controls")
	print(tr)
	cat(" Co/Tr Ratio =", round(nCo/nTr, 2), "\n\n")
	
	## groups and penalties
	cat(" Groups\n")
	print(cbind.data.frame("#Terms" = object$grouping, "Penalty" = round(object$group.penalty,1)))
	cat("\n")
	## terms and penalties
	if (print.level > 0){
		# terms 
		cat(" Terms\n")
		terms.dat <- data.frame("Group" = rep(names(object$grouping), times = (object$grouping), 
		"Penalty" = round(object$term.penalty,1)))
		rownames(terms.dat) <- colnames(object$mat)
		print(terms.dat)
		cat("\n")
	}

	# balance talbe
	cat(" Balance Table\n")
	print(object$bal.tab)	

	# return(out)
}