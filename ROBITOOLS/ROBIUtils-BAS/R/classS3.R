#' @include ROBIUtils.R
NULL

#' Adds a class into the class hierarchie attribute.
#' 
#' \code{addS3Class} adds a new class name to the vector
#' of class associated to the object. This the way to
#' assign an object to an S3 class. \code{addS3Class} add
#' the new class name in front of the class vector
#' 
#' @param object the object to modify
#' @param classname the name of the new class
#' 
#' @return the object given as parametter casted to the new
#'         class
#' 
#' @examples
#' x = c(1,3,2,5)
#' x = addS3Class(x,"my.vector")
#' class(x)
#' 
#' @seealso \code{\link{rmS3Class}}
#' 
#' @note for efficiency purpose no check is done on the input
#'       parametters
#'       
#' @keywords system function
#' 
#' @author Eric Coissac
#' @export
#' 
addS3Class = function(object,classname) {
  class(object) = c(classname,class(object))
  return(object)
}

#' Removes a class from the class hierarchie attribute.
#' 
#' \code{rmS3Class} removes a class name from the vector
#' of class associated to the object. This the way to
#' remove the association between an object and a S3 class. 
#' 
#' @param object the object to modify
#' @param classname the name of the class to remove
#' 
#' @return the object given as parametter.
#' 
#' @examples
#' x = c(1,3,2,5)
#' x = addS3Class(x,"my.vector")
#' class(x)
#' x = rmS3Class(x,"my.vector")
#' class(x)
#' 
#' @seealso \code{\link{addS3Class}}
#' 
#' @note for efficiency purpose no check is done on the input
#'       parametters
#'       
#' @keywords system function
#' 
#' @author Eric Coissac
#' @export
#' 
rmS3Class = function(object,classname) {
  c = class(object)
  if (! is.null(c))
    index = match(classname,c)
  class(object)=c[-index]
  return(object)
}

#' create basic functions to manipulate a new S3 class
#' 
#' createS3Class function create in the \code{package:ROBITools}
#' environment an \code{is.xxx} function and an \code{as.xxx} function
#' allowing to test if an abject belong the class \code{xxx} and to add
#' the class \code{xxx} to the class list of an object. \code{xxx} is a 
#' generic class name that is specified through the \code{classname}
#' argument of the function.
#' 
#' @param classname a \code{character string} indicating the name
#'         of the new class.
#'         
#' @examples
#' 
#' # Create a new S3 class named mynewclass
#' createS3Class('mynewclass')
#' 
#' #create a new vector object
#' x=c(1,4,6)
#' 
#' # test if it belongs the new class, that is false
#' is.mynewclass(x)
#' 
#' # Associate x to the new class
#' as.mynewclass(x)
#' 
#' # test again if x belongs the new class, that is now true
#' is.mynewclass(x)
#' 
#' @seealso \code{\link{rmS3Class}}
#' 
#' @note Take care that the new functions are created in the 
#' \code{package:ROBITools} environment.
#' 
#' @keywords system function
#' 
#' @author Eric Coissac
#' @export
#' 
createS3Class = function(classname,envir=NULL) {
  if (is.null(envir))
    envir=parent.frame()
  
  is.class = function(object) any(class(object)==classname)
  as.class = function(object) return(addS3Class(object,classname))
   
  assign(paste('is',classname,sep="."),is.class,envir=envir)
  assign(paste('as',classname,sep="."),as.class,envir=envir)
  
}



