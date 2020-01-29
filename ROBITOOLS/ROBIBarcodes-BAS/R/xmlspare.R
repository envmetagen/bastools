#'@include ROBIBarcodes.R
#'@import XML
#'
NULL

.__spare_tree__ = NULL

sparexmltree = function() {
  
  if (is.null(get(".__spare_tree__",envir = environment()))) {
  
    sparefile = paste(system.file("extdata", 
                                  package="ROBIBarcodes"),
                      'spare.xml',
                      sep='/')
    
    spare = xmlParseDoc(sparefile)
    
    assign(".__spare_tree__",spare, envir=globalenv())
  }
  
  return(get(".__spare_tree__",envir = globalenv()))
}
