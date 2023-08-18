
#' @importFrom utils packageDescription
#' @importFrom Rcpp evalCpp
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats uniroot
#' @importFrom Matrix bdiag

.onAttach = function(libname, pkgname)
{
    temp = packageDescription("Rmodule")
    msg = paste(temp$Package, ": ", temp$Title, "\n", "Version ", temp$Version, " created on ", temp$Date, ".\n", sep = "")
    msg = paste(msg, "Copyright (c) 2023 John Hughes\n", sep = "")
    msg = paste(msg, 'For citation information, type citation("Rmodule").\n', sep = "")
    msg = paste(msg, 'Type help(package = Rmodule) to get started.\n', sep = "")
    packageStartupMessage(msg)
}
