#'@include 02_class_metabarcoding.data.R
#'@import ROBITaxonomy 

NULL

#' Construct a taxonomic tree from a list of taxa
#' 
#' Construct a graph from a table containing the taxonomic path of sequences
#' 
#' 
#' @param x a table containing the taxonomic path of the references. Typically an output from get.classic.taxonomy
#' 
#' @return g a directed graph displaying the taxonomy hierarchy of the input data. Stored in a \code{\link{igraph}} object 
#'         where the taxonomic ranks of the vertices are added to the vertices attributes
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
#' graph.tax.termes = dbtree(termes.taxo.table[,1:7])
#' library(igraph)
#' 
#' #plot the tree
#' coord = layout.reingold.tilford(graph.tax.termes, root=1, circular=F)
#' v.cex = as.factor(V(graph.tax.termes)$rank)
#' levels(v.cex) = match(levels(v.cex), colnames(termes.taxo.table))
#' plot(graph.tax.termes, vertex.size=1, vertex.label.cex=2*(as.numeric(as.vector(v.cex))^-1), edge.arrow.size=0, layout=coord)
#' 
#' 
#' #Vizualization with sequence counts
#' tax.count = log10(colSums(termes$reads)[match(as.vector(V(graph.tax.termes)$name), termes$motus$scientific_name)])
#' tax.count[is.na(tax.count)|tax.count<0] = 0.01
#' V(graph.tax.termes)$count = unname(tax.count)
#' 
#' plot(graph.tax.termes, vertex.size=V(graph.tax.termes)$count, vertex.label.cex=2*(as.numeric(as.vector(v.cex))^-1), edge.arrow.size=0, layout=coord)
#' 
#'   
#' @seealso \code{\link{get.classic.taxonomy}}
#' @author Lucie Zinger
#' @export

dbtree = function(x) {
  
  #dealing with noranks
  x2 = x
  for (i in 1:ncol(x2)) {
    x2[,i] = as.character(x[,i])
    if(length(which(is.na(x[,i])==T))!=0) {
      if(i==1) {
        x2[which(is.na(x[,i])==T),i] = "NR"
      } else {
        x2[which(is.na(x[,i])==T),i] = as.character(x2[,i-1])[which(is.na(x2[,i])==T)]
      }
    }
  }
      
  #prepare an edgelist
  edgelist = list()
  
  for (i in 1:(ncol(x2)-1)){
    out = x2[,c(i,i+1)]
    out2 = out[-which(duplicated(out)==T),]
    colnames(out2) = c("parent", "kid")
    edgelist[[i]] = out2[which(out2[,1]!=out2[,2]),]
  } 
  
  edgelist = do.call("rbind", edgelist)
  
  
  #construct the graph
  
  g = igraph::graph.edgelist(as.matrix(edgelist), directed=T)
  
  #get taxorank for each taxa
  ranks = do.call("rbind", lapply(1:ncol(x), function(y) {
    out = cbind(unique(as.character(x[,y])), colnames(x)[y])
    out
  }))
  
  #Assign nodes to taxorank
  igraph::V(g)$rank = ranks[match(igraph::V(g)$name, ranks[,1]),2]
    
  return(g)
}