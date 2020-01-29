#' @include 02_class_metabarcoding.data.R
NULL

#' Detects contaminants in metabarcoding data
#' 
#' Detects sequences/motus in a \code{\link{metabarcoding.data}} object
#' for which frequencies over the entire dataset are maximum in negative controls and 
#' hence, most likely to be contaminants. 
#' 
#' 
#' @param x a \code{\link{metabarcoding.data}} object
#' @param controls a vector of samples names where conta are suspected to be detected 
#'                 (typically negative control names).
#' @param clust a vector for grouping sequences. Default set to \code{NULL}.
#'
#' @return a vector containing the names of sequences identified as contaminants
#'
#' @examples
#' 
#' data(termes)
#' termes.ok = termes[,colSums(termes$reads)>0]
#' neg = rownames(termes.ok)[grep("r",rownames(termes.ok))]
#' 
#' #finds contaminants based on neg samples
#' contaslayer(termes.ok, neg)
#' 
#' # extanding contamininant detection with grouping factor, 
#' # typically obiclean/sumatra cluster or taxonomy membership
#' contaslayer(termes.ok, neg, termes.ok$motus$scientific_name)
#'   
#' @seealso \code{\link{threshold}} for further trimming
#' @author Lucie Zinger
#' @export

contaslayer = function(x,controls,clust=NULL){
  
  x.fcol = normalize(x, MARGIN=2)$reads
  x.max = rownames(x.fcol[apply(x.fcol, 2, which.max),])
  conta = colnames(x)[!is.na(match(x.max,controls))]
  
  if (length(clust)!=0) {
    agg = data.frame(conta.id=colnames(x.fcol), clust)
    conta.ext = agg$conta.id[which(!is.na(match( agg$clust, agg$clust[match(conta,agg$conta.id)])))]
    return(as.vector(conta.ext))
  } 
  else {
    return(conta)
  }
}  
