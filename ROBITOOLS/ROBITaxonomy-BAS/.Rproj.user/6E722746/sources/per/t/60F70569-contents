#' @include taxonomy.R
NULL

#' @export
setGeneric("rank.list", function(taxonomy) {
  return(standardGeneric("rank.list"))
})

#' Returns the list of taxonomic ranks
#' 
#' The \code{rank.list} function returns the list of all taxonomic
#' ranks described in the taxonomy
#' 
#' @param   taxonomy a \code{\linkS4class{obitools.taxonomy}} instance
#' 
#' @return a vector of type \code{character} containing the taxonomic rank names
#'  
#' @examples
#' # read the taxonomy database
#' taxo=default.taxonomy()
#' 
#' # returns the taxonomic rank for all taxid between 1000 and 1020
#' rank.list(taxo)
#' 
#' @docType methods
#' @rdname rank.list-methods
#' @aliases rank.list-methods,obitools.taxonomy
#' @author Eric Coissac
#' @keywords taxonomy
#' 
setMethod("rank.list", "obitools.taxonomy",
          function(taxonomy) {
            return(.Call('R_rank_list',
                         taxonomy,
                         PACKAGE="ROBITaxonomy"))
          })


#' @export
setGeneric("taxonomicrank", function(taxonomy,taxid) {
  return(standardGeneric("taxonomicrank"))
})

#' Returns the taxonomic rank associated to a taxid
#' 
#' @param   taxonomy a \code{\linkS4class{obitools.taxonomy}} instance
#' @param   taxid a vector of taxid to analyse
#'
#' @return a vector of type \code{character} containing the taxonomic ranks
#'  
#' @examples
#' # read the taxonomy database
#' taxo=default.taxonomy()
#' 
#' # returns the taxonomic rank for all taxid between 1000 and 1020
#' taxonomicrank(taxo,1000:1020)
#' 
#' @docType methods
#' @rdname taxonomicrank-methods
#' @aliases taxonomicrank-methods,obitools.taxonomy
#' @author Eric Coissac
#' @keywords taxonomy
#' 
setMethod("taxonomicrank", "obitools.taxonomy",function(taxonomy,taxid) {  
  taxid = as.integer(taxid)
  return(.Call('R_get_rank',
               taxonomy,
               taxid,
               PACKAGE="ROBITaxonomy"))
})

#' @export
setGeneric("taxonatrank", function(taxonomy,taxid,rank,name=FALSE) {
  return(standardGeneric("taxonatrank"))
})

#' Extracts the taxid at a specified taxonomic rank.
#' 
#' The \code{taxonatrank} method of \code{\linkS4class{obitools.taxonomy}} class
#' returns the \emph{taxid} or the scientific name corresponding 
#' to a \emph{taxid}.at a specified taxonomic rank
#' 
#' @param   taxonomy a \code{\linkS4class{obitools.taxonomy}} instance
#' @param   taxid a vector of taxid to analyse
#' @param   rank a \code{character} indicating the desired rank
#' @param   name  A logical value \code{TRUE} or \code{FALSE} indicating 
#'          if the method return a taxid or a scientific name.
#'          
#' @return  \describe{
#'   \item{If \code{name==FALSE}}{the taxid of the corresponding 
#'                                taxon as an integer or a vector of integers 
#'                                if the \code{taxid} argument is itself 
#'                                a vector}
#'   \item{If \code{name==TRUE}}{the scientific name of the corresponding 
#'                               taxon as a string or a vector of string 
#'                               if the \code{taxid} argument is itself 
#'                               a vector}
#'  }
#'  
#' @examples
#' # read the taxonomy database
#' taxo=default.taxonomy()
#' 
#' # build a vector of 6 taxids corresponding to species
#' sp.taxid=c(7000,7004,7007,7009,7010,7011)
#' 
#' # look for the subfamily taxids
#' taxonatrank(taxo,sp.taxid,"subfamily")
#' 
#' # same thing but returns results as a vector of scientific names
#' taxonatrank(taxo,sp.taxid,"subfamily",TRUE)
#' 
#' @seealso class \code{\linkS4class{obitools.taxonomy}},
#' and methods \code{\link{species}},\code{\link{genus}},
#' \code{\link{family}},\code{\link{kingdom}},
#' \code{\link{superkingdom}}
#' 
#' @docType methods
#' @rdname taxonatrank-methods
#' @aliases taxonatrank-methods,obitools.taxonomy
#' @author Eric Coissac
#' @keywords taxonomy
#' 
setMethod("taxonatrank", "obitools.taxonomy",function(taxonomy,taxid,rank,name=FALSE) {
  getsp = function(t)   {	
    if (is.na(t[1]) | is.na(t[2]))
      return(NA)
    else
      return(.Call('R_findtaxonatrank',taxonomy,
                   as.integer(t[1]),
                   t[2],
                   name,
                   PACKAGE="ROBITaxonomy"))
  }
  
  rank  = as.character(rank)
  name  = as.logical(name[1])
  
  apply(data.frame(taxid,rank),1,getsp)
})


#' @export
setGeneric("species", function(taxonomy,taxid,name=FALSE) {
  return(standardGeneric("species"))
})

#' Extracts the species corresponding to a taxid
#' 
#' The \code{species} method of \code{\linkS4class{obitools.taxonomy}} class
#' returns the \emph{taxid} or the scientific name of the species corresponding 
#' to a \emph{taxid}.
#' 
#' @param   taxonomy a \code{\linkS4class{obitools.taxonomy}} instance
#' @param   taxid a vector of taxid to analyse
#' @param   name  A logical value \code{TRUE} or \code{FALSE} indicating 
#'          if the method return a taxid or a scientific name.
#'          
#' @return  \describe{
#'   \item{If \code{name==FALSE}}{the taxid of the corresponding 
#'                                taxon as an integer or a vector of integers 
#'                                if the \code{taxid} argument is itself 
#'                                a vector}
#'   \item{If \code{name==TRUE}}{the scientific name of the corresponding 
#'                               taxon as a string or a vector of string 
#'                               if the \code{taxid} argument is itself 
#'                               a vector}
#'  }
#'  
#' @examples
#' # read the taxonomy database
#' taxo=default.taxonomy()
#' 
#' # build a vector of 6 taxids corresponding to species
#' sp.taxid=c(7000,7004,7007,7009,7010,7011)
#' 
#' # look for the species taxids
#' species(taxo,sp.taxid)
#' 
#' # same thing but returns results as a vector of scientific names
#' species(taxo,sp.taxid,TRUE)
#' 
#' @seealso class \code{\linkS4class{obitools.taxonomy}},
#' and methods \code{\link{taxonatrank}},\code{\link{genus}},
#' \code{\link{family}},\code{\link{kingdom}},
#' \code{\link{superkingdom}}
#' 
#' @docType methods
#' @rdname species-methods
#' @aliases species-methods,obitools.taxonomy
#' @author Eric Coissac
#' @keywords taxonomy
#' 
setMethod("species", "obitools.taxonomy",function(taxonomy,taxid,name=FALSE) {
  getsp = function(t)	 {	
    if (is.na(t))
      return(NA)
    else
      return(.Call('R_get_species',
                   taxonomy,
                   t,
                   name,
                   PACKAGE="ROBITaxonomy"))
  }
  
  taxid = as.integer(taxid)
  name  = as.logical(name[1])
  sapply(taxid,getsp)
})

#' @export
setGeneric("genus", function(taxonomy,taxid,name=FALSE) {
  return(standardGeneric("genus"))
})

#' Extracts the genus corresponding to a taxid
#' 
#' The \code{genus} method of \code{\linkS4class{obitools.taxonomy}} class
#' returns the \emph{taxid} or the scientific name of the genus corresponding 
#' to a \emph{taxid}.
#' 
#' @param   taxonomy a \code{\linkS4class{obitools.taxonomy}} instance
#' @param   taxid a vector of taxid to analyse
#' @param   name  A logical value \code{TRUE} or \code{FALSE} indicating 
#'          if the method return a taxid or a scientific name.
#'          
#' @return  \describe{
#'   \item{If \code{name==FALSE}}{the taxid of the corresponding 
#'                                taxon as an integer or a vector of integers 
#'                                if the \code{taxid} argument is itself 
#'                                a vector}
#'   \item{If \code{name==TRUE}}{the scientific name of the corresponding 
#'                               taxon as a string or a vector of string 
#'                               if the \code{taxid} argument is itself 
#'                               a vector}
#'  }
#'  
#' @examples
#' # read the taxonomy database
#' taxo=default.taxonomy()
#' 
#' # build a vector of 6 taxids corresponding to species
#' sp.taxid=c(7000,7004,7007,7009,7010,7011)
#' 
#' # look for the genus taxids
#' genus(taxo,sp.taxid)
#' 
#' # same thing but returns results as a vector of scientific names
#' genus(taxo,sp.taxid,TRUE)
#' 
#' @seealso class \code{\linkS4class{obitools.taxonomy}},
#' and methods \code{\link{species}},\code{\link{taxonatrank}},
#' \code{\link{family}},\code{\link{kingdom}},
#' \code{\link{superkingdom}}
#' 
#' @docType methods
#' @rdname genus-methods
#' @aliases genus-methods,obitools.taxonomy
#' @author Eric Coissac
#' @keywords taxonomy
#' 
setMethod("genus", "obitools.taxonomy",function(taxonomy,taxid,name=FALSE) {
  getsp = function(t)	 {	
    if (is.na(t))
      return(NA)
    else
      return(.Call('R_get_genus',
                   taxonomy,
                   t,
                   name,
                   PACKAGE="ROBITaxonomy"))
  }
  
  taxid = as.integer(taxid)
  name  = as.logical(name[1])
  sapply(taxid,getsp)
})

#' @export
setGeneric("family", function(taxonomy,taxid,name=FALSE) {
  return(standardGeneric("family"))
})

#' Extracts the family corresponding to a taxid
#' 
#' The \code{family} method of \code{\linkS4class{obitools.taxonomy}} class
#' returns the \emph{taxid} or the scientific name of the family corresponding 
#' to a \emph{taxid}.
#' 
#' @param   taxonomy a \code{\linkS4class{obitools.taxonomy}} instance
#' @param   taxid a vector of taxid to analyse
#' @param   name  A logical value \code{TRUE} or \code{FALSE} indicating 
#'          if the method return a taxid or a scientific name.
#'          
#' @return  \describe{
#'   \item{If \code{name==FALSE}}{the taxid of the corresponding 
#'                                taxon as an integer or a vector of integers 
#'                                if the \code{taxid} argument is itself 
#'                                a vector}
#'   \item{If \code{name==TRUE}}{the scientific name of the corresponding 
#'                               taxon as a string or a vector of string 
#'                               if the \code{taxid} argument is itself 
#'                               a vector}
#'  }
#'  
#' @examples
#' # read the taxonomy database
#' taxo=default.taxonomy()
#' 
#' # build a vector of 6 taxids corresponding to species
#' sp.taxid=c(7000,7004,7007,7009,7010,7011)
#' 
#' # look for the family taxids
#' family(taxo,sp.taxid)
#' 
#' # same thing but returns results as a vector of scientific names
#' family(taxo,sp.taxid,TRUE)
#' 
#' @seealso class \code{\linkS4class{obitools.taxonomy}},
#' and methods \code{\link{species}},\code{\link{genus}},
#' \code{\link{taxonatrank}},\code{\link{kingdom}},
#' \code{\link{superkingdom}}
#' 
#' @docType methods
#' @rdname family-methods
#' @aliases family-methods,obitools.taxonomy
#' @author Eric Coissac
#' @keywords taxonomy
#' 
setMethod("family", "obitools.taxonomy",function(taxonomy,taxid,name=FALSE) {
  getsp = function(t)	 {	
    if (is.na(t))
      return(NA)
    else
      return(.Call('R_get_family',
                   taxonomy,
                   t,
                   name,
                   PACKAGE="ROBITaxonomy"))
  }
  
  taxid = as.integer(taxid)
  name  = as.logical(name[1])
  sapply(taxid,getsp)
})

#' @export
setGeneric("kingdom", function(taxonomy,taxid,name=FALSE) {
  return(standardGeneric("kingdom"))
})

#' Extracts the kingdom corresponding to a taxid
#' 
#' The \code{kingdom} method of \code{\linkS4class{obitools.taxonomy}} class
#' returns the \emph{taxid} or the scientific name of the kingdom corresponding 
#' to a \emph{taxid}.
#' 
#' @param   taxonomy a \code{\linkS4class{obitools.taxonomy}} instance
#' @param   taxid a vector of taxid to analyse
#' @param   name  A logical value \code{TRUE} or \code{FALSE} indicating 
#'          if the method return a taxid or a scientific name.
#'          
#' @return  \describe{
#'   \item{If \code{name==FALSE}}{the taxid of the corresponding 
#'                                taxon as an integer or a vector of integers 
#'                                if the \code{taxid} argument is itself 
#'                                a vector}
#'   \item{If \code{name==TRUE}}{the scientific name of the corresponding 
#'                               taxon as a string or a vector of string 
#'                               if the \code{taxid} argument is itself 
#'                               a vector}
#'  }
#'  
#' @examples
#' # read the taxonomy database
#' taxo=default.taxonomy()
#' 
#' # build a vector of 6 taxids corresponding to species
#' sp.taxid=c(7000,7004,7007,7009,7010,7011)
#' 
#' # look for the kingdom taxids
#' kingdom(taxo,sp.taxid)
#' 
#' # same thing but returns results as a vector of scientific names
#' kingdom(taxo,sp.taxid,TRUE)
#' 
#' @seealso class \code{\linkS4class{obitools.taxonomy}},
#' and methods \code{\link{species}},\code{\link{genus}},
#' \code{\link{family}},\code{\link{taxonatrank}},
#' \code{\link{superkingdom}}
#' 
#' @docType methods
#' @rdname kingdom-methods
#' @aliases kingdom-methods,obitools.taxonomy
#' @author Eric Coissac
#' @keywords taxonomy
#' 
setMethod("kingdom", "obitools.taxonomy",function(taxonomy,taxid,name=FALSE) {
  getsp = function(t)	 {	
    if (is.na(t))
      return(NA)
    else
      return(.Call('R_get_kingdom',
                   taxonomy,
                   t,
                   name,
                   PACKAGE="ROBITaxonomy"))
  }
  
  taxid = as.integer(taxid)
  name  = as.logical(name[1])
  sapply(taxid,getsp)
})

#' @export
setGeneric("superkingdom", function(taxonomy,taxid,name=FALSE) {
  return(standardGeneric("superkingdom"))
})

#' Extracts the superkingdom corresponding to a taxid
#' 
#' The \code{superkingdom} method of \code{\linkS4class{obitools.taxonomy}} class
#' returns the \emph{taxid} or the scientific name of the superkingdom corresponding 
#' to a \emph{taxid}.
#' 
#' @param   taxonomy a \code{\linkS4class{obitools.taxonomy}} instance
#' @param   taxid a vector of taxid to analyse
#' @param   name  A logical value \code{TRUE} or \code{FALSE} indicating 
#'          if the method return a taxid or a scientific name.
#'          
#' @return  \describe{
#'   \item{If \code{name==FALSE}}{the taxid of the corresponding 
#'                                taxon as an integer or a vector of integers 
#'                                if the \code{taxid} argument is itself 
#'                                a vector}
#'   \item{If \code{name==TRUE}}{the scientific name of the corresponding 
#'                               taxon as a string or a vector of string 
#'                               if the \code{taxid} argument is itself 
#'                               a vector}
#'  }
#'  
#' @examples
#' # read the taxonomy database
#' taxo=default.taxonomy()
#' 
#' # build a vector of 6 taxids corresponding to species
#' sp.taxid=c(7000,7004,7007,7009,7010,7011)
#' 
#' # look for the superkingdom taxids
#' superkingdom(taxo,sp.taxid)
#' 
#' # same thing but returns results as a vector of scientific names
#' superkingdom(taxo,sp.taxid,TRUE)
#' 
#' @seealso class \code{\linkS4class{obitools.taxonomy}},
#' and methods \code{\link{species}},\code{\link{genus}},
#' \code{\link{family}},\code{\link{kingdom}},
#' \code{\link{taxonatrank}}
#' 
#' @docType methods
#' @rdname superkingdom-methods
#' @aliases superkingdom-methods,obitools.taxonomy
#' @author Eric Coissac
#' @keywords taxonomy
#' 
setMethod("superkingdom", "obitools.taxonomy",function(taxonomy,taxid,name=FALSE) {
  getsp = function(t)	 {	
    if (is.na(t))
      return(NA)
    else
      return(.Call('R_get_superkingdom',
                   taxonomy,
                   t,
                   name,
                   PACKAGE="ROBITaxonomy"))
  }
  
  taxid = as.integer(taxid)
  name  = as.logical(name[1])
  sapply(taxid,getsp)
})


