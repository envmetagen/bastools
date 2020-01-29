#'@import ROBITaxonomy
#'@include ROBIBarcodes.R
NULL


#'@export
resolution = function(taxonomy,ecopcr) {
  l = aggregate(ecopcr$taxid,
                by=list(barcode=ecopcr$sequence),
                function(x) lowest.common.ancestor(taxonomy,x))
  r = taxonomicrank(taxonomy,l$x)
  names(r)=as.character(l$barcode)
  
  return(r[as.character(ecopcr$sequence)])
}