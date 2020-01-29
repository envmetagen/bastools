#' @include taxonomy.R
NULL


#
#
# Manage le loading of the default taxonomy
#
#

.__default__taxonomy__ = NULL

#' Returns the default taxonomy
#' 
#' Returns a \code{\linkS4class{obitools.taxonomy}} instance corresponding
#' to a NCBI taxonomy included by default in the \pkg{\link{ROBITaxonomy}} package.
#' 
#' @return a \code{\linkS4class{obitools.taxonomy}} instance.
#' 
#' @examples
#' 
#' # Load the default taxonomy
#' taxo = default.taxonomy()
#' 
#' # and use it for requesting a scientific name
#' scientificname(taxo,7742)
#' 
#' @seealso \code{\linkS4class{obitools.taxonomy}}
#' 
#' @author Eric Coissac
#' @keywords taxonomy
#' @export
#' 
default.taxonomy = function() {
  if (is.null(get(".__default__taxonomy__",envir = environment())))
    assign(".__default__taxonomy__", 
           read.taxonomy(paste(system.file("extdata", 
                                           package="ROBITaxonomy"),
                                           'ncbitaxo',
                                        sep='/')), 
           envir=globalenv())
  
  return(get(".__default__taxonomy__",envir = globalenv()))
}


#' @export
#' 
is.obitools.taxonomy = function(taxonomy) {
  class(t)[1] == "obitools.taxonomy"
}

