#Commentaires lus par roxygen
#'@include ROBIBarcodes.R
#'@import XML
#'
NULL

#NULL termine le commentaire pour roxygen

extractPrimers <-function(primer){
  
  id=xmlAttrs(primer)["ID"]
  name=xmlValue(xmlChildren(xmlChildren(primer)$name)$text)
  sequence=xmlValue(xmlChildren(xmlChildren(primer)$sequence)$text)
  coding=as.logical(xmlValue(xmlChildren(xmlChildren(primer)$coding)$text))
  
  p=list(id=id, name=name, sequence=sequence, coding=coding)
  
  return(p)
}

#Export pour rendre publique la fonction

#' Builds primer data frame from metabarcodedb
#' 
#' The \code{primers.data.frame} function extracts all the primer information
#' from the \code{metabarcodedb} database.
#' 
#' @param barcodedb a xml document containing a metabarcodedb.
#' 
#' @return a \code{data.frame} describing primers.
#' 
#' @examples
#' # load the XML library
#' library(XML)
#' 
#' # load the example metabarcodedb database
#' db = xmlParseDoc(system.file("extdata/barcodedb.xml", package="ROBIBarcodes"))
#' 
#' # extracts the primer table
#' primers.data.frame(db)
#' 
#' @author Aurelie Bonin
#' @keywords metabarcodes
#' 
#' @export
primers.data.frame <-function(barcodedb){
  p=getNodeSet(db, 
               path="/obi:obimetabarcodedb/obi:primers/obi:primer" , 
               namespaces=c(obi="http://metabarcoding.org/OBIMetabarcodes"))

  primerTable=as.data.frame(do.call(rbind,lapply(p,extractPrimers)))

  rownames(primerTable)=primerTable$id
  primerTable=primerTable[,-1]
  
  return(primerTable)
}
