#' @include ROBIBarcodes.R
#' @import  ROBITaxonomy
#' @import XML
#' @useDynLib ROBIBarcodes
#'
NULL


extractTaxa <-function(taxon){
  
  id=xmlAttrs(taxon)["ID"]
  name=xmlValue(xmlChildren(xmlChildren(taxon)$name)$text)
  rank=xmlValue(xmlChildren(xmlChildren(taxon)$rank)$text)
  partof=xmlValue(xmlChildren(xmlChildren(taxon)$partof)$text)
  
  p=list(id=id, name=name, rank=rank, partof=partof)
  
  return(p)
}

#' Builds taxa data frame from metabarcodedb
#' 
#' The \code{taxonomy.data.frame} function extracts all the taxon information
#' from the \code{metabarcodedb} database.
#' 
#' @param barcodedb a xml document containing a metabarcodedb.
#' 
#' @return a \code{data.frame} describing taxa.
#' 
#' @examples
#' # load the XML library
#' library(XML)
#' 
#' # load the example metabarcodedb database
#' db = xmlParseDoc(system.file("extdata/barcodedb.xml", package="ROBIBarcodes"))
#' 
#' # extracts the taxonomy table
#' taxonomy.data.frame(db)
#' 
#' @author Eric Coissac
#' @keywords metabarcodes
#' 
#' @export
taxonomy.data.frame = function(barcodedb) {
  p=getNodeSet(db, 
               path="/obi:obimetabarcodedb/obi:taxonomy/obi:taxon" , 
               namespaces=c(obi="http://metabarcoding.org/OBIMetabarcodes"))
  
  taxonomyTable=as.data.frame(do.call(rbind,lapply(p,extractTaxa)))
  
  
  rownames(taxonomyTable)=unlist(taxonomyTable$id)
  taxonomyTable=taxonomyTable[,-1]
  
  taxonomyTable$name=unlist(taxonomyTable$name)
  taxonomyTable$rank=unlist(taxonomyTable$rank)
  taxonomyTable$partof=unlist(taxonomyTable$partof)

  return(taxonomyTable)
  
}

#' Builds a \code{taxonomy.obitools} from a metabarcodedb
#' 
#' The \code{metabarcodedb.taxonomy} function extracts all the taxon information
#' from the \code{metabarcodedb} database and create a \code{taxonomy.obitools}
#' instance with them.
#' 
#' @param barcodedb a xml document containing a metabarcodedb.
#' 
#' @return a \code{taxonomy.obitools} instance.
#' 
#' @examples
#' # load the XML library
#' library(XML)
#' 
#' # load the example metabarcodedb database
#' db = xmlParseDoc(system.file("extdata/barcodedb.xml", package="ROBIBarcodes"))
#' 
#' # extracts the taxonomy table
#' barcodetaxo = metabarcodedb.taxonomy(db)
#' 
#' # Look for the Verbrata taxid 
#' ecofind(barcodetaxo,"vertebrata")
#' 
#' @author Eric Coissac
#' @keywords metabarcodes
#' 
#' @export
metabarcodedb.taxonomy = function(barcodedb) {
  
  table = taxonomy.data.frame(barcodedb)
  
  t <- .Call('R_buildbarcodetaxo',table,TRUE,PACKAGE="ROBIBarcodes")	
  
  return(ROBITools:::build.taxonomy.obitools(t,"barcodedb",getwd(),FALSE))
}
