#' A package to manipulate DNA metabarcoding data.
#' 
#' This package was written as a following of the OBITools. 
#' 
#' \tabular{ll}{
#'  Package: \tab ROBITools\cr
#'  Type: \tab Package\cr
#'  Version: \tab 0.1\cr
#'  Date: \tab 2013-06-27\cr
#'  License: \tab CeCILL 2.0\cr
#'  LazyLoad: \tab yes\cr
#'}
#' 
#' @name ROBITools-package
#' @aliases ROBITools
#' @docType package
#' @title  A package to manipulate DNA metabarcoding data.
#' @author Frederic Boyer
#' @author Aurelie Bonin
#' @author Lucie Zinger
#' @author Eric Coissac
#' 
#' @references http://metabarcoding.org/obitools
#' 
NA

.onLoad <- function(libname, pkgname) { 
  
  packageStartupMessage( "ROBITools package" )
  #print(getwd())
  
}

