#' Make taxon tables from obitab-megan merged table
#' @param readCount The name of the column containing readcounts. If no readcounts available, leave NULL
#' @param obimeg Output of \code{merge.tab.taxon}
#'
#' @return List of taxon tables as dataframes, one for each taxonomic level
#'@examples
#'
#'@export
taxtablise<-function(obimeg,readCount=NULL){
  taxonLevel<-c("SK","P","C","O","F","G","S")
  taxa.tables<-list()
  if(is.null(readCount)){
    #couldnt figure out loop! think its because of quoted characters in formula
  taxa.tables[["SK"]]<-data.table::dcast(data = obimeg,formula = SK~barcode,fun.aggregate = length,value.var = "read")
  taxa.tables[["P"]]<-data.table::dcast(data = obimeg,formula = P~barcode,fun.aggregate = length,value.var = "read")
  taxa.tables[["C"]]<-data.table::dcast(data = obimeg,formula = C~barcode,fun.aggregate = length,value.var = "read")
  taxa.tables[["O"]]<-data.table::dcast(data = obimeg,formula = O~barcode,fun.aggregate = length,value.var = "read")
  taxa.tables[["F"]]<-data.table::dcast(data = obimeg,formula = F~barcode,fun.aggregate = length,value.var = "read")
  taxa.tables[["G"]]<-data.table::dcast(data = obimeg,formula = G~barcode,fun.aggregate = length,value.var = "read")
  taxa.tables[["S"]]<-data.table::dcast(data = obimeg,formula = S~barcode,fun.aggregate = length,value.var = "read")
  }

  if(!is.null(readCount)){
    message("not tested!")
    taxa.tables[["SK"]]<-data.table::dcast(data = obimeg,formula = SK~barcode,fun.aggregate = sum,value.var = "read")
    taxa.tables[["P"]]<-data.table::dcast(data = obimeg,formula = P~barcode,fun.aggregate = sum,value.var = "read")
    taxa.tables[["C"]]<-data.table::dcast(data = obimeg,formula = C~barcode,fun.aggregate = sum,value.var = "read")
    taxa.tables[["O"]]<-data.table::dcast(data = obimeg,formula = O~barcode,fun.aggregate = sum,value.var = "read")
    taxa.tables[["F"]]<-data.table::dcast(data = obimeg,formula = F~barcode,fun.aggregate = sum,value.var = "read")
    taxa.tables[["G"]]<-data.table::dcast(data = obimeg,formula = G~barcode,fun.aggregate = sum,value.var = "read")
    taxa.tables[["S"]]<-data.table::dcast(data = obimeg,formula = S~barcode,fun.aggregate = sum,value.var = "read")
    }

    return(taxa.tables)
}

#' Plot taxon tables
#'@examples
#'
#'@export
basplot.taxtablise<-function(taxon.tables){
  melted<-list()
  melted<-lapply(X = taxon.tables,FUN = function(x)reshape2::melt(x))

  melted2<-list()
  melted2<-lapply(X=melted,FUN = function(x)x[!is.na(x[,1]),])

  kingdom.plot<-ggplot2::ggplot(data=melted2[[1]] , aes(y=value, x=variable, fill=SK))+
    geom_bar(stat = "identity")+
    theme(legend.title = element_text(size=5), legend.text=element_text(size=5),axis.text.x=element_text(size=5,angle=45, hjust=1))
  phylum.plot<-ggplot2::ggplot(data=melted2[[2]] , aes(y=value, x=variable, fill=P))+
    geom_bar(stat = "identity")+
    theme(legend.title = element_text(size=5), legend.text=element_text(size=5),axis.text.x=element_text(size=5,angle=45, hjust=1))
  class.plot<-ggplot2::ggplot(data=melted2[[3]] , aes(y=value, x=variable, fill=C))+
    geom_bar(stat = "identity")+
    theme(legend.title = element_text(size=5), legend.text=element_text(size=5),axis.text.x=element_text(size=5,angle=45, hjust=1))
  order.plot<-ggplot2::ggplot(data=melted2[[4]] , aes(y=value, x=variable, fill=O))+
    geom_bar(stat = "identity")+
    theme(legend.title = element_text(size=5), legend.text=element_text(size=5),axis.text.x=element_text(size=5,angle=45, hjust=1))
  family.plot<-ggplot2::ggplot(data=melted2[[5]] , aes(y=value, x=variable, fill=F))+
    geom_bar(stat = "identity")+
    theme(legend.title = element_text(size=5), legend.text=element_text(size=5),axis.text.x=element_text(size=5,angle=45, hjust=1))
  genus.plot<-ggplot2::ggplot(data=melted2[[6]] , aes(y=value, x=variable, fill=G))+
    geom_bar(stat = "identity")+
    theme(legend.title = element_text(size=5), legend.text=element_text(size=5),axis.text.x=element_text(size=5,angle=45, hjust=1))
  species.plot<-ggplot2::ggplot(data=melted2[[7]] , aes(y=value, x=variable, fill=S))+
    geom_bar(stat = "identity")+
    theme(legend.title = element_text(size=5), legend.text=element_text(size=5),axis.text.x=element_text(size=5,angle=45, hjust=1))

  plot.list<-list(kingdom.plot,phylum.plot,class.plot,order.plot,family.plot,genus.plot,species.plot)
  return(plot.list)
}
