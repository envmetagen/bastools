#' @import ROBITaxonomy
#' @include 02_class_metabarcoding.data.R
NULL

#' Dataset taxonomic resolution summary.
#' 
#' Summarizes the taxonomic relution of reads and MOTUs over the entire dataset
#' 
#' 
#' @param x          a \code{\link{metabarcoding.data}} object
#' @param colranks   a string indicating column name where ranks are stored in \code{x}
#' @param colscores  a string indicating column name where taxonomic identification scores are stored in \code{x}
#' @param thresh     a threshold for defining at which taxonomic identification scores a sequence can be considered as "not assigned". 
#'                   Default is \code{0.7}
#' 
#' @return returns a data.frame and piecharts of the number/proportion of MOTUs/reads assigned to each taxonomic levels
#' 
#' @examples
#' 
#' data(termes)
#' taxo=default.taxonomy()
#' 
#' termes.taxo.table = get.classic.taxonomy(termes, taxo, "taxid")
#' attr(termes, "motus") = data.frame(termes$motus, termes.taxo.table)
#' attr(termes, "motus")["count"] = colSums(termes$reads)
#' 
#' summary.taxores(termes, "taxonomic_rank_ok","best_identity")
#' 
#' @seealso \code{\linkS4class{taxonomy.obitools}}, and method \code{\link{taxonmicank}}
#'
#' @author Lucie Zinger
#' @keywords taxonomy
#' 
#' @export
#' 
summary.taxores = function(x,colranks,colscores, thresh=0.7){
  
  #vector encompassing all ranked possible taxonomic levels
  taxorank = c("superkingdom", "kingdom", "subkingdom", "superphylum", "phylum", "subphylum", "superclass", "class", "subclass", "infraclass",
               "superorder", "order", "suborder", "infraorder", "parvorder", "superfamily", "family", "subfamily", "supertribe", "tribe", 
               "subtribe", "supergenus", "genus", "subgenus", "species group", "species subgroup", "superspecies", "species", "subspecies",
               "varietas", "forma", "no rank", "not assigned")
  
  #settings if thresh
  ranks = as.vector(x$motus[,colranks])
  ranks[x$motus[,colscores]<thresh] =  "not assigned"
  
  #nb of otus
  tmp = table(ranks)
  taxores.otu = tmp[match(taxorank, names(tmp))]
  names(taxores.otu) = taxorank
  taxores.otu[is.na(taxores.otu)] = 0
  
  #nb of reads
  tmp = aggregate(x$motus$count, by=list(ranks), sum)
  taxores.reads = tmp[match(taxorank,tmp[,1]),2]
  names(taxores.reads) = taxorank
  taxores.reads[is.na(taxores.reads)] = 0
  
  #plot
  layout(matrix(c(1,2,1,3),2,2),heights=c(0.3,1))
  col.tmp = c(rainbow(length(taxorank)-2,start=0, end=0.5, alpha=0.6), "lightgrey", "darkgrey")
  par(mar=c(1,0,0,0), oma=c(0,0,2,0))
  frame()
  legend("top", taxorank, ncol=6, cex=0.8, fill=col.tmp)
  pie(taxores.otu, col=col.tmp, border="lightgrey", labels="", clockwise=T)
  mtext("OTUs", side=1, cex=1)
  pie(taxores.reads, col=col.tmp, border="lightgrey", labels="", clockwise=T)
  mtext("Reads", side=1, cex=1)
  
  #result
  out = data.frame(otu=taxores.otu, reads=taxores.reads)
  out
}
