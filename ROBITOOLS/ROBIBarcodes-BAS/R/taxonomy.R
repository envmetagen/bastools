#'@include xmlMods.R
#'@import XML
#'@import ROBITaxonomy
#'@include ROBIBarcodes.R
NULL


taxon.data.frame  = function(taxonomy,taxids,strict=TRUE,known.taxid=c()) {
  taxids = as.integer(sub("TX.","",as.character(taxids)))
  good.taxid =  validate(taxonomy,taxids)
  
  if (strict & any(is.na(good.taxid))) 
    stop(sprintf("The following taxids are absent from the taxonomy : %s",
                 toString(taxids[is.na(good.taxid)])))
  
  good.taxid = good.taxid[! is.na(good.taxid)]
  all.path   = path(taxonomy,good.taxid)
  all.taxid  = Reduce(union,all.path)
  all.taxid  = sort(union(all.taxid,known.taxid))[-1]
  all.parent = sprintf("TX.%d",parent(taxonomy,all.taxid))
  all.rank   = taxonomicrank(taxonomy,all.taxid)
  all.scientificname = scientificname(taxonomy,all.taxid)
  
  all.id = sprintf("TX.%d",all.taxid)
  
  rep = data.frame(taxid=all.id,
                   name=all.scientificname,
                   rank=all.rank,
                   partof=all.parent,
                   stringsAsFactors=FALSE)
  
  return(rep)
}

build.taxon.node = function(taxid,name,rank,partof) {
  nodes = lapply(sprintf('\n<taxon ID="%s"><name>%s</name><rank>%s</rank><partof>%s</partof></taxon>',
                     taxid,
                     name,
                     rank,
                     partof),
                 xmlParseString)
  
  
  return(nodes)
  
}

#'@export
add.taxon.barcodedb = function(barcodedb,
                                   taxonomy,
                                   taxids) {
  
  taxonomy.node = getNodeSet(barcodedb,
                        path='/obi:obimetabarcodedb/obi:taxonomy',
                        c(obi="http://metabarcoding.org/OBIMetabarcodes"))[[1]]
  
  known.taxid = as.character(
                  getNodeSet(
                    taxonomy.node,
                    path="./obi:taxon/@ID",
                    c(obi="http://metabarcoding.org/OBIMetabarcodes")))
    
  known.taxid = as.integer(sub("TX.","",known.taxid))
  
  taxon = taxon.data.frame(taxonomy,taxids,strict=TRUE,known.taxid)
  
  taxon.nodes = c(xmlChildren(taxonomy.node)$root,
                  build.taxon.node(taxon$taxid,
                                 taxon$name,
                                 taxon$rank,
                                 taxon$partof))
  spare = sparexmltree()
  new.taxonomy.node =  getNodeSet(spare,
                                  path='/obi:obimetabarcodedb/obi:taxonomy',
                                  c(obi="http://metabarcoding.org/OBIMetabarcodes"))[[1]]
  
  replaceNodes(taxonomy.node,new.taxonomy.node)
  
  hidden = addChildren(new.taxonomy.node,kids=taxon.nodes,append=FALSE)
}

