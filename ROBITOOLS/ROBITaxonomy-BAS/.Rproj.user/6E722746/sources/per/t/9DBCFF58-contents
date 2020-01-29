#' @include taxonomy.R
NULL

#' @export
setGeneric("scientificname", function(taxonomy,taxid) {
  return(standardGeneric("scientificname"))
})

#' Returns the scientific name corresponding to a \emph{NCBI taxid}
#' 
#' \code{scientificname} function in package \pkg{\link{ROBITaxonomy}} returns the 
#' scientific name corresponding to a \emph{NCBI taxid}.
#' 
#' @param   taxonomy a \code{\link{obitools.taxonomy}} instance
#' @param   taxid an integer value or a vector of integer representing NCBI 
#'          taxonomic identifiers.
#' @return  The scientific name of the corresponding taxon as a string or a 
#'          vector of string if the \code{taxid} argument is itself a vector
#' 
#' @examples
#' # load the default taxonomy database include in the ROBITaxonomy library
#' taxo=default.taxonomy()
#' 
#' # build a vector of 6 taxids corresponding to species
#' sp.taxid=c(7000,7004,7007,7009,7010,7011)
#' 
#' # look for the scientific names correponding to these taxids
#' scientificname(taxo,sp.taxid)
#' 
#' @seealso class \code{\linkS4class{obitools.taxonomy}}
#'          
#'          
#' @docType methods
#' @rdname scientificname-methods
#' @aliases scientificname-methods,obitools.taxonomy
#' @author Eric Coissac
#' 
setMethod("scientificname", "obitools.taxonomy",function(taxonomy,taxid) {
  getscname = function(t)  {	
    if (is.na(t))
      return(NA)
    else
      return(	.Call('R_get_scientific_name',
                    taxonomy,
                    t,
                    PACKAGE="ROBITaxonomy"))
  }
  
  taxid = as.integer(taxid)
  sapply(taxid,getscname)
})


######################################################################
######################################################################


#' @export
setGeneric("parent", function(taxonomy,taxid,name=FALSE) {
  return(standardGeneric("parent"))
})

#' Returns the parent taxon corresponding to a \emph{NCBI taxid}
#' 
#' \code{parent} function in package \pkg{\link{ROBITaxonomy}} returns the 
#' parent taxon corresponding to a \emph{NCBI taxid}.
#' 
#' @param   taxonomy a \code{\link{obitools.taxonomy}} instance
#' @param   taxid an integer value or a vector of integer representing NCBI 
#'          taxonomic identifiers.
#' @param   name A logical value \code{TRUE} or \code{FALSE} indicating if the 
#'          method returns a taxid or a scientific name.
#'                            
#' @return  \describe{ \item{If \code{name==FALSE}}{the taxid of the 
#'                            parent taxon as an integer or a vector of 
#'                            integers if the \code{taxid} argument is itself 
#'                            a vector} 
#'                      \item{If \code{name==TRUE}}{the scientific name of the
#'                            parent taxon as a string or a vector of 
#'                            string if  the \code{taxid} argument is itself a 
#'                            vector} }
#' 
#' @examples
#' # load the default taxonomy database include in the ROBITaxonomy library
#' taxo=default.taxonomy()
#' 
#' # build a vector of 6 taxids corresponding to species
#' sp.taxid=c(7000,7004,7007,7009,7010,7011)
#' 
#' # look for the parent taxa correponding to these taxids
#' parent(taxo,sp.taxid)
#' 
#' # same things but scientific names are returned
#' parent(taxo,sp.taxid,TRUE)
#' 
#' @seealso class \code{\linkS4class{obitools.taxonomy}}
#'          
#'          
#' @docType methods
#' @rdname parent-methods
#' @aliases parent-methods,obitools.taxonomy
#' @author Eric Coissac
#' 
setMethod("parent", "obitools.taxonomy",function(taxonomy,taxid,name=FALSE) {
  getp = function(t)	 {	
    if (is.na(t))
      return(NA)
    else
      return(.Call('R_get_parent',
                   taxonomy,
                   as.integer(t),
                   name,
                   PACKAGE="ROBITaxonomy"))
  }
  
  taxid = as.integer(taxid)
  name  = as.logical(name[1])
  sapply(taxid,getp)
})



######################################################################
######################################################################



#' @export
setGeneric("taxid.list", function(taxonomy) {
  return(standardGeneric("taxid.list"))
})


#' Returns the list of all taxids belonging the taxonomy.
#'  
#' \code{taxid.list} returns the list of all taxids included in the
#' instance of the class \code{\linkS4class{obitools.taxonomy}}
#' 
#' @param taxonomy the \code{\linkS4class{obitools.taxonomy}} to use.
#' 
#' @return an \code{integer} vector containing the list of taxids.
#' 
#' @examples
#' # loads the default taxonomy database
#' taxo=default.taxonomy()
#' 
#' # returns the count of taxa described in the taxonomy
#' length(taxo)
#' 
#' # extracts the list of all valid taxids
#' good = taxid.list(taxo)
#' 
#' # returns the size of the returned list
#' length(good)
#' 
#' @seealso \code{\linkS4class{obitools.taxonomy}}
#' 
#' @author Eric Coissac
#' @keywords taxonomy
#' @docType methods
#' @rdname taxid.list-method
#' @aliases taxid.list
#' 
setMethod("taxid.list", "obitools.taxonomy",
          function(taxonomy) {
            return(.Call('R_taxid_list',
                         taxonomy,
                         PACKAGE="ROBITaxonomy"))
          })

######################################################################
######################################################################



#' Returns the count of taxa in the taxonomy.
#'  
#' \code{length} returns the count of taxa included in the
#' instance of the class \code{\linkS4class{obitools.taxonomy}}
#' 
#' @param x the \code{\linkS4class{obitools.taxonomy}} to use.
#' 
#' @return an \code{integer} corresponding to the count of taxa.
#' 
#' @examples
#' # loads the default taxonomy database
#' taxo=default.taxonomy()
#' 
#' # returns the count of taxa described in the taxonomy
#' length(taxo)
#' 
#' @seealso \code{\link{length}}, \code{\linkS4class{obitools.taxonomy}}
#' 
#' @author Eric Coissac
#' @keywords taxonomy
#' @export  length.obitools.taxonomy
#' 
length.obitools.taxonomy = function(x) 
{
  return(.Call('R_length_taxonomy',
               x,
               PACKAGE="ROBITaxonomy"))
}


######################################################################
######################################################################

setGeneric('max')  	

#' Returns the maximum taxid in the taxonomy.
#'  
#' \code{length} returns the maximum taxid included in the
#' instance of the class \code{\linkS4class{obitools.taxonomy}}
#' 
#' @param taxonomy the \code{\linkS4class{obitools.taxonomy}} to use.
#' @param na.rm included for compatibility purpose, this parameter as
#'              no effect on this implementation of \code{max}
#' 
#' @return an \code{integer} corresponding to the count of taxa.
#' 
#' @examples
#' # load the default taxonomy database
#' taxo=default.taxonomy()
#' 
#' # gets the larger taxid of the database
#' max(taxo)
#' 
#' @seealso \code{\link{max}}, \code{\linkS4class{obitools.taxonomy}}
#' 
#' @author Eric Coissac
#' @keywords taxonomy
#' @export  max.obitools.taxonomy
#' 
max.obitools.taxonomy=function(taxonomy,na.rm = FALSE) {
  return(.Call('R_max_taxid',
               taxonomy,
               PACKAGE="ROBITaxonomy"))
}

#' @export
setGeneric("ecofind", function(taxonomy,patterns,rank=NULL,alternative=FALSE) {
  return(standardGeneric("ecofind"))
})

#' Returns taxids associated to the names
#' 
#' Return the set of taxids having their name matching the given pattern.
#' 
#' @param taxonomy the \code{\linkS4class{obitools.taxonomy}} to use.
#' @param patterns one or several regular pattern used to select the the taxa.
#' @param   rank a \code{character} indicating a taxonomic rank. If not \code{NULL}
#'               only taxids correponding to this rank are returned.
#' @param   alternative  A logical value \code{TRUE} or \code{FALSE} indicating 
#'          if the function must only look for a scientific name.
#' 
#' @return if just one pattern is given, an integer vector is returned with the
#'         corresponding taxids. If a list of patterns is given, the function
#'         returns a list of integer vectors, each vector containing the taxids
#'         corresponding to a pattern. The returned list is in the same order 
#'         than the given patern list.
#' 
#' @examples
#' # load the default taxonomy database
#' taxo=default.taxonomy()
#' 
#' # retreives the Vertebrata taxid
#' taxid = ecofind(taxo,"Vertebrata")
#' 
#' taxid
#' scientificname(taxo,taxid)
#' 
#' 
#' taxid = ecofind(taxo,"^Vertebrata$")
#' 
#' taxid
#' scientificname(taxo,taxid)
#' 
#'
#' @author Eric Coissac
#' @keywords taxonomy
#' @docType methods
#' @rdname ecofind-method
#' @aliases ecofind,obitools.taxonomy
#' 
setMethod("ecofind", "obitools.taxonomy",function(taxonomy,patterns,rank=NULL,alternative=FALSE) {
  getp = function(t)   {	
    if (is.na(t))
      return(NA)
    else
      return(unique(.Call('R_ecofind',
                          taxonomy,
                          t,
                          rank,
                          alternative,
                          PACKAGE="ROBITaxonomy")))
  }
  
  patterns = as.character(patterns)
  taxid=lapply(patterns,getp)
  if (length(taxid)==1)
    taxid=taxid[[1]]
  
  return(taxid)
})



#' @export
setGeneric("validate", function(taxonomy,taxid) {
  return(standardGeneric("validate"))
})

#' Checks that a \emph{taxid} is really present in taxonomy
#' 
#' \code{validate} function in package \pkg{\link{ROBITaxonomy}} checks 
#' that a \emph{taxid} is declared in the considered taxonomy.
#' 
#' @param   taxonomy a \code{\link{obitools.taxonomy}} instance
#' @param   taxid an integer value or a vector of integer representing NCBI 
#'          taxonomic identifiers.
#'                            
#' @return  The taxid if it exists, NA otherwise. If the input taxid is a 
#'          vector of integer returns an integer vector composed of validated
#'          taxids and NA values.
#' 
#' @examples
#' # load the default taxonomy database include in the ROBITaxonomy library
#' taxo=default.taxonomy()
#' 
#' # build a vector of 101 taxids 
#' sp.taxid=c(7000:7100)
#' 
#' # checks the list of taxids
#' validate(taxo,sp.taxid)
#' 
#' @seealso class \code{\linkS4class{obitools.taxonomy}}
#'          
#'          
#' @docType methods
#' @rdname validate-methods
#' @aliases validate-methods,obitools.taxonomy
#' @author Eric Coissac
#' 
setMethod("validate", "obitools.taxonomy",function(taxonomy,taxid) {
  getp = function(t)   {  
    if (is.na(t))
      return(NA)
    else
      return(.Call('R_validate_taxid',
                   taxonomy,
                   t,
                   PACKAGE="ROBITaxonomy"))
  }
  
  taxid = as.integer(taxid)
  sapply(taxid,getp)
})

