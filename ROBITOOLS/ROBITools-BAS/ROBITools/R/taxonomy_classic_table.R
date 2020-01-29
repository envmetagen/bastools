#' @import ROBITaxonomy
#' @include 02_class_metabarcoding.data.R
NULL

#' Get classical taxonomy format
#' 
#' Creates a table with the classical taxonomic description (from phylum to species)
#' 
#' @param x a \code{\link{metabarcoding.data}} object
#' @param taxonomy a instance of \code{\linkS4class{taxonomy.obitools}}
#' @param coltaxid a the name of the column containing taxids to be used for creating classical taxonomic description
#' 
#' @return returns a data.frame with the classical taxonomic description ("kingdom", "phylum", "class", "order", "family", "genus", "species"), as well as
#'         sequence taxonomic assignment rank and scientific name for each sequences stored in the \code{\link{metabarcoding.data}} object
#' 
#' @examples
#' 
#' data(termes)
#' 
#' taxo=default.taxonomy()
#' 
#' termes.taxo.table = get.classic.taxonomy(termes, taxo, "taxid")
#' head(termes.taxo.table)
#' 
#' attr(termes, "motus") = data.frame(termes$motus, termes.taxo.table)
#' 
#'   
#' @seealso \code{\linkS4class{taxonomy.obitools}}, and methods \code{\link{species}},\code{\link{genus}}, \code{\link{family}},\code{\link{kingdom}},
#'          \code{\link{superkingdom}},\code{\link{taxonatrank}}, \code{\link{taxonmicank}}
#'
#' @author Lucie Zinger
#' @keywords taxonomy
#' @export
#' 

get.classic.taxonomy = function(x, taxonomy, coltaxid) {
  
  classic.taxo = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  
  taxids = x$motus[,coltaxid]
  
  out = as.data.frame(do.call("cbind", lapply(classic.taxo, function(y) {
    scientificname(taxonomy, taxonatrank(taxonomy,taxids,y))
  })))
  
  colnames(out) = paste(classic.taxo, "_name_ok", sep="")
  rownames(out) = colnames(x)
  
  out$scientific_name_ok = scientificname(taxonomy, taxids)
  out$taxonomic_rank_ok = taxonomicrank(taxonomy, taxids)

  return(out)
}