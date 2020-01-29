#' @include ROBITaxonomy.R
#' @useDynLib ROBITaxonomy
NULL

#' Gives access to a taxonomy preformated by OBITools
#' 
#' A S4 class describing a taxonomy. It allows access to 
#' taxonomy formated for OBITools.
#'
#' @references \describe{
#'  \item{NCBI Taxonomy : }{\url{http://www.ncbi.nlm.nih.gov/taxonomy}}
#'  \item{OBITools : }{\url{http://metabarcoding/obitools/doc}}
#'  }
#'  
#' @seealso \code{\link{read.taxonomy}}
#'  
#' @name obitools.taxonomy
#' @rdname obitools-taxonomy-class
#' @keywords taxonomy
#' @author Eric Coissac
#' @exportClass obitools.taxonomy
#' 
setClass("obitools.taxonomy",
         
         
         #
         # Attribute declaration
         #
         
         # data.frame containing the counts of reads per samples
         #      1 samples per line
         #      1 sequence per column
         
         representation(
           
           # An external pointer structure to 
           # the C taxonomy structure
           
           pointer    = "externalptr",
           
           # the name of the database on the hard disk
           
           dbname     = 'character',
           
           # the working directory when the taxonomy
           # object is created. 
           # This inforation combined with bname allows
           # to reload taxonomy from disk
           
           workingdir = 'character',
           
           # Indicate if the taxonomy is saved in a file
           # Taxonomy created in R or modified in R are
           # not saved ==> This have to be take into 
           # consideration but how ???
           saved = 'logical'
         ),
         
         #
         # Check object structure 
         #
         
         validity = function(object) { 			
           return(TRUE)
         }
)



#' obitools.taxonomy constructor
#' 
#' --> this constructor have not to be called directly
#'     use the read.obitools.taxonomy function to
#'     create a new instance of taxonomy
#'
#' @docType methods
#' @rdname initialize-methods-obitools.taxonomy
#' @aliases initialize-methods,obitools.taxonomy
setMethod("initialize",
          "obitools.taxonomy",
          function(.Object, pointer,dbname,workingdir,saved) {
            .Object@pointer    <- pointer
            .Object@dbname     <- dbname
            .Object@workingdir <- workingdir
            .Object@saved      <- saved
            
            validObject(.Object) ## valide l'objet
            return(.Object)
          })


#' @exportClass obitools.taxonomyOrNULL
setClassUnion("obitools.taxonomyOrNULL",c("obitools.taxonomy","NULL"))


#' @export
setGeneric("path", function(taxonomy,taxid,name=FALSE) {
  return(standardGeneric("path"))
})



setMethod("path", "obitools.taxonomy",function(taxonomy,taxid,name=FALSE) {
  getp = function(t)	 {	
    if (is.na(t))
      return(NA)
    else
    {
      path=c()
      
      t=.Call('R_validate_taxid',
              taxonomy,
              as.integer(t),
              PACKAGE="ROBITaxonomy")
      
      if (is.na(t))
        return(NA)
      
      repeat {
        if (name)
          path = c(scientificname(taxonomy,t),path)
        else
          path = c(t,path)
        
        t = .Call('R_get_parent',
                  taxonomy,
                  t,
                  FALSE,
                  PACKAGE="ROBITaxonomy")
        if (is.na(t))
          break
      }
      
      return(path)
    }
  }
  
  taxid=as.integer(taxid)
  name=as.logical(name)
  
  p = lapply(taxid,getp)
  d = dim(p)
  
  if (!is.null(d))
    if (d[2]==1)
      p = as.vector(p)
  
  return(p)
})


#' @export
setGeneric("is.subcladeof", function(taxonomy,taxid,parent) {
  return(standardGeneric("is.subcladeof"))
})

setMethod("is.subcladeof", "obitools.taxonomy",function(taxonomy,taxid,parent) {
  taxid = as.integer(taxid)
  parent= as.integer(parent)
  return(.Call('R_is_under_taxon',
               taxonomy,
               taxid,
               parent,
               PACKAGE="ROBITaxonomy"))
})



build.taxonomy = function(pointer,dbname,workingdir,saved) {
  rd <- new('obitools.taxonomy',
            pointer=pointer,
            dbname=dbname,
            workingdir=workingdir,
            saved=saved
  )
  return(rd)
}


#' Reads a taxonomy
#' 
#' \code{read.taxonomy} reads a taxonomy formated by OBITools.
#' NCBI taxonomy can be download from the NCBI FTP site in taxdump format.
#' The taxdump must be formated using the obitaxonomy command from OBITools
#' before being used in R. A OBITools formated taxonomy is composed of 3 files 
#' with the same prefix name and suffixes .tdx, .rdx, .ndx, two extra files 
#' suffixed .adx and .ldx can also be present. 
#' 
#' @param  dbname A character string containing the file name of the database
#' 
#' @return an instance of the class \code{\linkS4class{obitools.taxonomy}}
#' 
#' @examples
#' 
#' \dontshow{# switch the working directory to the data package directory}
#' \dontshow{setwd(system.file("extdata", package="ROBITaxonomy"))}
#'
#' # read the taxonomy ncbi
#' ncbi = read.taxonomy("ncbitaxo")
#' 
#' # and use it for requesting a scientific name
#' scientificname(ncbi,7742)
#' 
#' @seealso \code{\linkS4class{obitools.taxonomy}}
#' 
#' @author Eric Coissac
#' @keywords taxonomy
#' @export

read.taxonomy = function(dbname) {
  t <- .Call('R_read_taxonomy',dbname,TRUE,PACKAGE="ROBITaxonomy")	
  
  return(build.taxonomy(t,dbname,getwd(),TRUE))
}


