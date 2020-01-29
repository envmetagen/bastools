#' @include taxonomy.R
NULL

#' @export
setGeneric("longest.path", function(taxonomy,taxid) {
  return(standardGeneric("longest.path"))
})

#' Returns the longuest path from a taxon.
#' 
#' The method \code{longest.path} returns the length of the
#' path linking a taxid to the farest leaf belonging this taxid.
#' 
#' @param taxonomy the \code{\linkS4class{obitools.taxonomy}} to use.
#' 
#' @param taxid an \code{integer} vector containing the list of taxids.
#' 
#' @return an \code{integer} vector containing the list length.
#' 
#' @examples
#' # loads the default taxonomy database
#' taxo=default.taxonomy()
#' 
#' # returns the longest path in the taxonomy (from the root node)
#' longest.path(taxo,1)
#' 
#' 
#' @seealso \code{\linkS4class{obitools.taxonomy}}
#' 
#' @author Eric Coissac
#' @keywords taxonomy
#' @docType methods
#' @rdname longest.path-method
#' @aliases longest.path,obitools.taxonomy
#' 
setMethod("longest.path", "obitools.taxonomy",
          function(taxonomy,taxid) {
            getp = function(t)   {	
              if (is.na(t))
                return(NA)
              else
                return(.Call('R_longest_path',
                             taxonomy,
                             t,
                             PACKAGE="ROBITaxonomy"))
            }
            
            taxid = as.integer(taxid)
            sapply(taxid,getp)
          })


#' @export
setGeneric("distance.taxonomy", function(taxonomy,taxid1,taxid2=NULL,name=F) {
  return(standardGeneric("distance.taxonomy"))
})


#' Computes a distance matrix between taxids
#' 
#' The method \code{taxonomy.distance} computes a distance matrix between a
#' set of taxids. The distance between two taxa is based on the topology of
#' the taxonomomy tree.
#' 
#' \deqn{ d(Taxon_A,Taxon_B) = \frac{longest.path(lca(Taxon_A,Taxon_B))}{max(longest.path(Taxon_A),longest.path(Taxon_B))}}
#'                                  { longest.path(lca(Taxon_A,Taxon_B)) / max(longest.path(Taxon_A),longest.path(Taxon_B)) }
#' 
#' 
#' @param taxonomy the \code{\linkS4class{obitools.taxonomy}} to use.
#' 
#' @param taxid1 an \code{integer} vector containing a list of taxids.
#' 
#' @param taxid2 an \code{integer} vector containing a list of taxids.
#'               If \code{taxid2} is set to \code{NULL} (it's default value)
#'               then the \code{taxid2} list is considered as equal to 
#'               \code{taxid1} list. 
#' @param name  A logical value \code{TRUE} or \code{FALSE} indicating 
#'          if the method return distance matrix annotated by taxids or 
#'          by scientific names.
#'          
#' @return the distance matrix between taxids specified in the \code{taxid1}
#'         set and the \code{taxid2} set.
#'         
#' @examples
#' # loads the default taxonomy database
#' taxo=default.taxonomy()
#' 
#' # build a vector of 6 taxids corresponding to species
#' sp.taxid=c(7000,7004,7007,7009,7010,7011)
#' 
#' # computes the distance matrix between taxids
#' distance.taxonomy(taxo,sp.taxid)
#' 
#' # Same thing but the matrix is annotated by scientific names
#' distance.taxonomy(taxo,sp.taxid,name=TRUE)
#' 
#' @seealso \code{\link{longest.path}}
#' 
#' @author Eric Coissac
#' @keywords taxonomy
#' @docType methods
#' @rdname distance.taxonomy-method
#' @aliases taxonomy.distance,obitools.taxonomy
#' 
setMethod("distance.taxonomy", "obitools.taxonomy",
          function(taxonomy,taxid1,taxid2=NULL,name=F) {
            taxdist = function(r)
            {	
              t1=r[1]
              t2=r[2]
              if (is.na(t1) | is.na(t2))
                return(NA)
              
              p1 = path(taxonomy,t1)
              p2 = path(taxonomy,t2)
              
              minp = min(length(p1),length(p2))
              common = sum(p1[1:minp] == p2[1:minp])
              lca = p1[common]
              lp = longest.path(taxonomy,lca)
              return(lp/(lp+common))
            }
            
            multitaxdist=function(t1,t2) {
              apply(data.frame(t1,t2),1,taxdist)
            }
            
            taxid1 = taxid1[! is.na(validate(taxonomy,taxid1))]
            t1 = path(taxonomy,taxid1)
            
            same = is.null(taxid2)
            
            if (same)
            {
              ntaxon = length(taxid1) 
              t2 = t1[unlist(sapply(2:ntaxon,
                                    function(x) x:ntaxon))]
              t1 = t1[rep(1:(ntaxon-1),(ntaxon-1):1)]
            }
            else
            {
              taxid2 = taxid2[! is.na(validate(taxonomy,taxid2))]
              t2 = path(taxonomy,taxid2)
              nt1 = length(taxid1)
              nt2 = length(taxid2)
              t1 = t1[rep(1:nt1,nt2)]
              t2 = t2[rep(1:nt2,rep(nt1,nt2))]
            }
            
            lmin = mapply(function(a,b) min(length(a),length(b)),
                          t1,
                          t2)
            
            llca = mapply(function(x,y,l) sum(x[1:l]==y[1:l]),
                          t1,
                          t2,
                          lmin)
            
            lb   = longest.path(taxonomy,mapply(function(x,y) x[y],t1,llca))
            d    = as.double(lb / (lb + llca))
            
            if (same) {
              attr(d, "Size") <- ntaxon
              if (name)
                attr(d, "Labels") <- scientificname(taxonomy,taxid1)
              else
                attr(d, "Labels") <- as.character(taxid1)
              attr(d, "Diag") <- FALSE
              attr(d, "Upper") <- FALSE
              attr(d, "method") <- NULL
              attr(d, "call") <- match.call()
              class(d) <- "dist"
            }
            else {
              if (name)
                d = matrix(d,nt1,nt2,
                           dimnames=list(scientificname(taxonomy,taxid1),
                                         scientificname(taxonomy,taxid2)))
              else
                d = matrix(d,nt1,nt2,
                           dimnames=list(as.character(taxid1),
                                         as.character(taxid2)))
              
            }
            
            return(d)
          })


