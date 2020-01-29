#' @include 02_class_metabarcoding.data.R
NULL

#
#
# Managment of layers
#
# Layers a matrix or factors with the same dimension
# than the read matrix
#

# get motus data.frames

#' @export
setGeneric("layer.names", function(obj) {
  return(standardGeneric("layer.names"))
})

#' Returns the names of all the layers 
#' 
#' \code{layer.names} extracts the list of all the layer
#' names attached to a \code{\link{metabarcoding.data}} instance.
#' 
#' @param   obj a \code{\link{metabarcoding.data}} instance
#' @return  a vector of type \code{character} containing the
#'          list of all the layer names.
#' 
#' @docType methods
#' @rdname layer.names-methods
#' @aliases layer.names-methods,metabarcoding.data
#' 
setMethod("layer.names", "metabarcoding.data", function(obj) {
  return(names(obj@layers))
})


#' Returns the a layer associated to a \code{\link{metabarcoding.data}}
#' 
#' [[ operator Extracts a layer
#' attached to a \code{\link{metabarcoding.data}} instance.
#' 
#' @usage \method{[[}{unmutable}(x,i)
#'
#' @param   x a \code{\link{metabarcoding.data}} instance
#' @return  matrix or a factor.
#' 
#' @docType methods
#' @rdname double-open-brace-methods
#' @aliases double-open-brace-methods,metabarcoding.data
#' @method [[
#' @export
#' 
setMethod("[[", "metabarcoding.data", 
          function(x, i, j, ...) {
            
            if (! is.character(i))
              stop('Just named index must be used')  
            
            if (i=="reads")
              return(x@reads)
            
            if (i=="samples")
              return(x@samples)
            
            if (i=="motus")
              return(x@motus)
            
            if (i=="reads")
              return(x@reads)
            
            return(x@layers[[i,exact=TRUE]])
          })

#' @method $
#' @export
setMethod("$", "metabarcoding.data", 
          function(x, name) {
            return(x[[name]])
          })


# set one data layer data.frames

#' @method [[<-
#' @export
setMethod("[[<-","metabarcoding.data", 
          function(x, i, j, ...,value) {
            
            if (any(dim(value)!=c(x@scount,x@mcount)))
              stop("data dimmension are not coherent with this metabarcoding.data")
            
            if (hasArg('j'))
              stop('Just one dimension must be specified')
            
            if (! is.character(i))
              stop('Just named index must be used')
            
            if (i=='reads')
              stop('you cannot change the reads layer by this way')
            
            if (i=='motus' | i=='samples')
              stop('layers cannot be names motus or samples')
            
            value = as.factor.or.matrix(value)
            rownames(value)=rownames(x@reads)
            colnames(value)=colnames(x@reads)
            x@layers[[i]]=value
            
            return(x)
          })

#' @method $<-
#' @export
setMethod("$<-","metabarcoding.data", 
          function(x, name, value) {
            
            x[[name]]=value
            return(x)
          })
