#'@include xmlMods.R
#'@import XML
#'
NULL

#' Tests if the metabarcode database has at least one bibliography reference
#' 
#' @export
has.bibliography = function(barcodedb) {
  length(getNodeSet(barcodedb,
                    path='/obi:obimetabarcodedb/obi:bibliography',
                    c(obi="http://metabarcoding.org/OBIMetabarcodes")))>0
}

#' @export
add.reference.barcodedb  = function(barcodedb,bibfile,bibutils='bib2xml') {
  
  if (! has.bibliography(barcodedb)) {
    # We create the bibliography node
    
    metabarcode = getNodeSet(barcodedb,
                          path='/obi:obimetabarcodedb',
                          c(obi="http://metabarcoding.org/OBIMetabarcodes"))[[1]]
    
    spare = sparexmltree()
    
    bibliography = getNodeSet(spare,
                          path='/obi:obimetabarcodedb/obi:bibliography',
                          c(obi="http://metabarcoding.org/OBIMetabarcodes"))[[1]]
    
    bibliography = xmlClone(bibliography)
    
    addChildren(metadata,bibliography,
                at = NA)
    

  }

  bibliography=getNodeSet(barcodedb,
                          path='/obi:obimetabarcodedb/obi:bibliography',
                          c(obi="http://metabarcoding.org/OBIMetabarcodes"))
  
  ref = bib2mods(bibfile,bibutils)
  hidden=addChildren(bibliography[[1]],kids=ref) 
}