#' @include 02_class_metabarcoding.data.R
NULL

# TODO: Add comment
# 
# Author: coissac
###############################################################################

#' @export
setGeneric("extracts.obiclean", function(obj) {
			return(standardGeneric("extracts.obiclean"))
		})

#' Extracts the obiclean results
#' 
#' The method \code{extracts.obiclean} of the class \code{\linkS4class{metabarcoding.data}}
#' extracts \code{obiclean} results from the MOTUs descriptions include in the 
#' \code{\linkS4class{metabarcoding.data}} instance. 
#' When an \code{obitab} file is imported using the \code{\link{import.metabarcoding.data}}
#' if \code{obiclean} results are present in the file they are stored in the 
#' \code{motu} data.frame. By calling this methods, MOTU descriptors describing 
#' the \code{obiclean} status are moved to a set of layers.
#' 
#' @param obj the \code{\linkS4class{metabarcoding.data}} to analyze
#' 
#' @return the modified \code{\linkS4class{metabarcoding.data}} instance
#' 
#' @examples
#' 
#' # load termite data set from the ROBITools sample data
#' data(termes)
#' 
#' # shows the initial list of layer names
#' layer.names(t)
#' 
#' # extracts the obiclean status
#' termes = extracts.obiclean(termes)
#' 
#' # shows the name of the newly created layers
#' layer.names(t)
#' 
#' 
#' 
#' @seealso \code{\linkS4class{metabarcoding.data}}, \code{\link{threshold.mask}}, \code{\link{normalize}}
#'
#' @docType methods
#' @rdname extracts-obiclean-methods
#' @aliases extracts.obiclean-methods,metabarcoding.data
#' @author Eric Coissac
#' 


setMethod("extracts.obiclean", "metabarcoding.data", function(obj) {
  
			pat = "^obiclean_status:.*$"
			cols = colnames(obj@motus)
			cleancols = grep(pat,cols)
			clean.names=cols[cleancols]
      p = grep(pat,cols)
      d = t(as.factor.or.matrix(obj@motus[,p]))
      n = sapply(strsplit(cols[p],':'),function(y) y[[2]])
      rownames(d)=n
      d = d[rownames(obj@reads),]
      obj[["obiclean_status"]]=d
      
      newmotus = obj@motus[-cleancols]
      
			pat = "^obiclean_count:.*$"
			cols = colnames(newmotus)
			cleancols = grep(pat,cols)
			clean.names=cols[cleancols]
			p = grep(pat,cols)
			d = t(as.factor.or.matrix(newmotus[,p]))
			n = sapply(strsplit(cols[p],':'),function(y) y[[2]])
			rownames(d)=n
			d = d[rownames(obj@reads),]
			obj[["obiclean_count"]]=d
			
			newmotus = newmotus[-cleancols]
			
			pat = "^obiclean_cluster:.*$"
			cols = colnames(newmotus)
			cleancols = grep(pat,cols)
			clean.names=cols[cleancols]
			p = grep(pat,cols)
			d = t(as.factor.or.matrix(newmotus[,p]))
			n = sapply(strsplit(cols[p],':'),function(y) y[[2]])
			rownames(d)=n
			d = d[rownames(obj@reads),]
			obj[["obiclean_cluster"]]=d
			
			newmotus = newmotus[-cleancols]
			
			newdata = copy.metabarcoding.data(obj,motus=newmotus)
      
			return(newdata)
		})


#' @export
setGeneric("extracts.obiclean_cluster", function(obj) {
  return(standardGeneric("extracts.obiclean_cluster"))
})

setMethod("extracts.obiclean_cluster", "metabarcoding.data", function(obj) {
 
    obiclean = extracts.obiclean(obj)
    obihead  = obiclean[,! is.na(obiclean$motus$obiclean_head)]
    obihead$obiclean_count[is.na(obihead$obiclean_count)]=0
    reads = obihead$obiclean_count
    
    l = obihead@layers[layer.names(obihead) != "obiclean_count"]
    
    newdata = copy.metabarcoding.data(obihead,reads=reads,layers=l)
    
    return(newdata)
}
)