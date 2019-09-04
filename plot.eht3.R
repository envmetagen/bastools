#' @export
basplot.clean.ecopcrhits<-function(clean.ecopcrhits,specify_taxon,specify_level=NULL){
  Hits_temp<-clean.ecopcrhits
  Hits_temp$phylum_name<-as.character(Hits_temp$phylum_name)
  Hits_temp$kingdom_name<-as.character(Hits_temp$kingdom_name)
  Hits_temp$genus_name<-as.character(Hits_temp$genus_name)
  Hits_temp$family_name<-as.character(Hits_temp$family_name)
  Hits_temp$species_name<-as.character(Hits_temp$species_name)
  if(is.null(specify_level)){
  #look up taxonlevel
  col.number<-grep(pattern = specify_taxon,x = Hits_temp,ignore.case = T)
  if(length(col.number)>1) stop("Taxon name exists in more than one level, try specifying taxon level")
  }
  if(!is.null(specify_level)){

    specify_level<-paste0(specify_level,"_name")
    col.number<-match(x = specify_level,table = colnames(Hits_temp))
  }

  b<-Hits_temp[Hits_temp[,col.number]==specify_taxon,]

  kingdomLength<-length(na.omit(unique(b$kingdom_name)))
  if(kingdomLength==0) stop("No taxa found")
  phylumLength<-length(na.omit(unique(b$phylum_name)))
  classLength<-length(na.omit(unique(b$class_name)))
  orderLength<-length(na.omit(unique(b$order_name)))
  familyLength<-length(na.omit(unique(b$family_name)))
  genusLength<-length(na.omit(unique(b$genus_name)))
  speciesLength<-length(na.omit(unique(b$species_name)))

  e<-b[!is.na(b$forward_match),]

  Kingdom<-kingdomLength*colSums(table(e$kingdom_name,e$total_mismatches))/sum(colSums(table(e$kingdom_name,e$total_mismatches)))
  Phylum<-phylumLength*colSums(table(e$phylum_name,e$total_mismatches))/sum(colSums(table(e$phylum_name,e$total_mismatches)))
  Class<-classLength*colSums(table(e$class_name,e$total_mismatches))/sum(colSums(table(e$class_name,e$total_mismatches)))
  Order<-orderLength*colSums(table(e$order_name,e$total_mismatches))/sum(colSums(table(e$order_name,e$total_mismatches)))
  Family<-familyLength*colSums(table(e$family_name,e$total_mismatches))/sum(colSums(table(e$family_name,e$total_mismatches)))
  Genus<-genusLength*colSums(table(e$genus_name,e$total_mismatches))/sum(colSums(table(e$genus_name,e$total_mismatches)))
  Species<-speciesLength*colSums(table(e$species_name,e$total_mismatches))/sum(colSums(table(e$species_name,e$total_mismatches)))
  GG<-as.data.frame(rbind(Kingdom,Phylum,Class,Order,Family,Genus,Species))
  GG$Tlevel<-rownames(GG)
  GG<-reshape2::melt(GG)
  colnames(GG)<-gsub(pattern = "variable",replacement = "Mismatch_No.",x = colnames(GG))

  ggplot2::ggplot(data=GG,ggplot2::aes(y=value, x=Tlevel,group=Mismatch_No.)) +
    ggplot2::geom_col(ggplot2::aes(fill = Mismatch_No.))+
    ggplot2::scale_x_discrete(limits=c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))+
    ggplot2::labs(y="No. of taxa amplified", x="Taxonomic level")+
    ggplot2::ggtitle(specify_taxon,subtitle = "No. of taxa amplified at each taxonomic level. Primer mismataches
indicate the total number of mismatches across both forward and reverse primers")
}

