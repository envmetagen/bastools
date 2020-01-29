#'@include ROBIBarcodes.R
#'@import XML
#'
NULL

#' Creates a new empty metabarcode database.
#' 
#' @export
metabarcodedb = function() {
  emptyfile = paste(system.file("extdata", 
                                package="ROBIBarcodes"),
                    'empty.xml',
                    sep='/')

 empty = xmlParseDoc(emptyfile)

 return(empty)
}