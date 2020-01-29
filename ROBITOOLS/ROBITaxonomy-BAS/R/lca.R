#' @include taxonomy.R
NULL

#' @export
setGeneric("lowest.common.ancestor", function(taxonomy,taxid,threshold=1.0,error=0,name=FALSE) {
  return(standardGeneric("lowest.common.ancestor"))
})

#' Computes the lowest common ancestor in the taxonomy tree between a set of taxa
#' 
#' The \code{lowest.common.ancestor} function in package \pkg{ROBITaxonomy} computes
#' the lowest common ancestor of a set of taxids. The lowest common ancestor (LCA)
#' is the most precise taxonomic group shared by all the considered taxa. Tha
#' \code{lowest.common.ancestor} function implemented in the \pkg{ROBITaxonomy}
#' package, considers a fuzzy definition of the LCA as the most precise
#' taxonomic group shared by a quorum of the considered taxa.
#' 
#' @param taxonomy an instance of \code{\linkS4class{obitools.taxonomy}}
#' @param taxid an integer value or a vector of integer representing NCBI 
#'   taxonomic identifiers.
#' @param threshold a numeric value between 0.0 and 1.0 indicating the minimum 
#' quorum of taxid that must belong the LCA.
#' @param error an integer value indicating the maximum count of taxids that 
#' have not to belong the returned taxid. A \code{threshold} below 1.0 have
#'  priority on the \code{error} parameter.
#' @param name A logical value \code{TRUE} or \code{FALSE} indicating if the
#' method return a \emph{taxid} or a scientific name.
#' 
#' @return Depending on the value of the \code{name} argument, set by default 
#'         to \code{FALSE} the method returns :
#'    \describe{  
#'           \item{If \code{name==FALSE}}{ the taxid of the taxon corresponding 
#'             to the LCA as an integer value}
#'           \item{If \code{name==TRUE}}{ the scientific name of the taxon 
#'                 corresponding to the LCA as a string}
#'          }
#'          
#' @examples
#' require(ROBITaxonomy)
#'
#' \dontshow{# switch the working directory to the data package directory}
#' \dontshow{setwd(system.file("extdata", package="ROBITaxonomy"))}
#'   
#' # read the taxonomy database
#'   
#' taxo=read.taxonomy('ncbitaxo')
#' 
#' # build a vector of 6 taxids corresponding to species
#' 
#' sp.taxid=c(7000,7004,7007,7009,7010,7011)
#' 
#' # look for the lowest common ancestor taxids
#' 
#' lowest.common.ancestor(taxo,sp.taxid)
#' 
#' # same thing but returns results as a vector of scientific names
#' lowest.common.ancestor(taxo,sp.taxid,name=TRUE)
#' 
#' # If we accept than 2 or 1 taxa do not belong the LCA
#' lowest.common.ancestor(taxo,sp.taxid,name=TRUE,error=2)
#' lowest.common.ancestor(taxo,sp.taxid,name=TRUE,error=1)
#' 
#' # Partial LCA can also be speciefied as the minimal frequency of
#' # taxa belonging the LCA 
#' lowest.common.ancestor(taxo,sp.taxid,name=TRUE,threshold=0.8)
#' 
#' @seealso class \code{\linkS4class{obitools.taxonomy}},
#' and methods \code{\link{path}}, \code{\link{parent}},
#' 
#' @author Eric Coissac
#' @keywords taxonomy
#' @docType methods
#' @rdname lowest.common.ancestor-method
#' @aliases lowest.common.ancestor,obitools.taxonomy
#' 
setMethod("lowest.common.ancestor", "obitools.taxonomy",
          function(taxonomy,taxid,threshold=1.0,error=0,name=FALSE) {
            
            if (threshold != 1.0)
              error=as.integer(floor(length(taxid) * (1-threshold)))
            
            
            #
            # Remove nod valid taxid
            #
            
            taxid = validate(taxonomy,taxid)
            if (any(is.na(taxid)))
              return(NA)
                      
            ntaxid=length(taxid)
            nok = ntaxid - error
            if (ntaxid==1)
              return(taxid)
            
            allpath  = path(taxonomy,taxid)
            minlength= min(vapply(allpath,length,0))
                        
            lca=NA
            for (i in 1:minlength) {
              
              n = vapply(allpath,function(x) x[i],0)
              nt = table(n)
              mt = max(nt)
              if (mt >= nok) {
                p = nt[nt==mt]
                if (length(p)==1)
                  lca=as.integer(names(p)[1])
                else
                  break
              }
              else
                  break
            }
            
            if (name)
              return(scientificname(taxonomy,lca))
            else
              return(lca)
          })


