#' @title Subsidiary hbal Function
#' @description Function to load package description.
#' @param lib   libname
#' @param pkg   package name
#' @importFrom utils packageDescription
#' @references Xu, Y., & Yang, E. (2021). Hierarchically Regularized Entropy Balancing.

.onAttach <- function(lib, pkg) {
  mylib <- dirname(system.file(package = pkg))
  title <- packageDescription(pkg, lib.loc = mylib)$Title
  ver <- packageDescription(pkg, lib.loc = mylib)$Version
  #author <- packageDescription(pkg, lib.loc = mylib)$Author
  packageStartupMessage(pkg, ": ", title, "\nVersion: ", ver, "\nReference: ", 
    "Xu, Y., & Yang, E. (2021). Hierarchically Regularized Entropy Balancing. Political Analysis, Forthcoming.", "\n")
}