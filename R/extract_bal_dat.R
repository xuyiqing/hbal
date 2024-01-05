#' @title Extracting balancing statistics
#' @aliases extract_bal_dat
#' @description \code{extract_bal_dat} For more flexibility with making  \code{hbal}. 
#' @usage extract_bal_dat(hbalobject)
#' @param hbalobject  an object of class \code{hbal} as returned by \code{hbal}.
#' @details this is a simple function to return a data.frame from the bal.tab . 
#' @return A data table of balance statistitics
#' @importFrom stringr str_detect
#' @importFrom purrr pluck
#' @importFrom dplyr filter rename case_when
#' @author Josh Allen
#' @examples
#' #EXAMPLE 1
#' set.seed(1984)
#' N <- 500
#' X1 <- rnorm(N)
#' X2 <- rbinom(N,size=1,prob=.5)
#' X <- cbind(X1, X2)
#' treat <- rbinom(N, 1, prob=0.5) # Treatment indicator
#' y <- 0.5 * treat + X[,1] + X[,2] + rnorm(N) # Outcome
#' dat <- data.frame(treat=treat, X, Y=y)
#' out <- hbal(Treat = 'treat', X = c('X1', 'X2'), Y = 'Y', data=dat)
#' dout <- bal_plot_dat(out)
#' ggplot(dout, aes(x = value, y = term, shape = adjust, color = adjust)) +
#' geom_point(alpha = 0.5) +
#' geom_vline(xintercept = 0, alpha = 0.5, linetype = "dotted") +
#' facet_wrap(vars(covar.group))
#' @export



extract_bal_dat = function(hbal_object){
  if(!inherits(hbal_object, "hbal")){
    stop(paste0("This needs to be an hbal object"))
  } else{
    
    bal_plot_dat = hbal_object |> 
      pluck("bal.tab") |> 
      as.data.frame.table() |> 
      rename(term = Var1,
             measure = Var2,
             value = Freq) |> 
      filter(str_detect(measure, "Std")) |> 
      mutate(adjust = case_when(measure == "Std.Diff.(O)" ~ "Before Adjustment",
                                measure == "Std.Diff.(W)" ~ "After Adjustment"),
             covar.group = rep(rep(names(third_order$grouping), third_order$grouping), 2))
    
    
    return(bal_plot_dat)  
  }
}

