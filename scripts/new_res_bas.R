ecopcroutput<-"COI.fullgb.ecopcrResults.txt.clean"

#########FINISH THIS



add.res.bas3<-function(ecopcroutput,ncbiTaxDir){
  message("Assumes 'sequence' is comprised only of insert")
  
  ecopcr<-data.table::fread(ecopcroutput,sep = "\t",data.table = F)
  
  colnames(ecopcr)<-gsub("taxid","taxids",colnames(ecopcr))
  
  ecopcr<-add.lineage.df(ecopcr,ncbiTaxDir)
  
  ecopcr$path<-paste(ecopcr$K,ecopcr$P,ecopcr$C,ecopcr$O,ecopcr$F,ecopcr$G,ecopcr$S,sep = ";")
  
  ecopcr$dummy<-1
  
  lca.step1<-aggregate(ecopcr$dummy,by=list(ecopcr$sequence),FUN=sum)
  
  lca.include<-lca.step1[lca.step1$x!=1,]
  
  ecopcr_to_do<-ecopcr[ecopcr$sequence %in% lca.include$Group.1,]
  
  #ecopcr_to_do<-ecopcr_to_do[order(ecopcr_to_do$sequence),]
  
  #test with 1
  # first if identical paths
  
  
  lca.bas<-function(paths){
    
    
  }