source("/home/tutorial/TOOLS/bastools/master_functions.R")
library(ggplot2)

fill_NA<-100
threshold=c(98,95,92)
simMat<-data.table::fread("/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/FILTURB/SARA_MSc/Alignments/North_Paper_DB_similarity_matrix.csv",header = T,data.table = F)
ncbiTaxDir = "TOOLS/DBS/ncbi_taxonomy/taxdump/"

a<-barcode_gap_analysis(simMat = simMat,threshold = threshold,fill_NA = 100,ncbiTaxDir="/home/tutorial/TOOLS/DBS/ncbi_taxonomy/taxdump/")
a[[1]]
a2<-cowplot::plot_grid(plotlist = a[[2]],labels = "AUTO") #list("species level","genus level","family level")
cowplot::save_plot("12SV51.cowplot.png",a2,base_width = 15,base_height = 10)
barcode_gap_analysis<-function(simMat,threshold=NULL,fill_NA=NULL,ncbiTaxDir){

#simMat is a dataframe matrix of identities (i.e. from an alignment), with the first column being taxids
#simMat<-data.table::fread("/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/FILTURB/SARA_MSc/Alignments/North_Paper_DB_similarity_matrix.csv",header = T,data.table = F)
#fill_NA is used to fill Nas with any value
  #threshold=c(98,95,92)
  
colnames(simMat)[1]<-"taxids"
simMat[is.na(simMat)]<-fill_NA
taxids<-simMat[,"taxids",drop=F]
lineage<-add.lineage.df(taxids,ncbiTaxDir = ncbiTaxDir)

taxlevel<-c("species","genus","family")
newtab.list<-list()
plot.list<-list()

for(j in 1:3){
  
  if(taxlevel[j]=="species") taxcol="S"
  if(taxlevel[j]=="genus") taxcol="G"
  if(taxlevel[j]=="family") taxcol="F"

  taxa<-lineage[match(simMat$taxids,lineage$taxids),taxcol]
  
  simMat2<-simMat
  simMat2[,1]<-taxa
  colnames(simMat2)<-c("taxids",taxa)
  long<-reshape2::melt(simMat2)
  
  #ranges 
  newtab<-data.frame(matrix(ncol=5,nrow = length(unique(long$taxids))))
  colnames(newtab)<-c("taxon","max_within","min_within","max_among","min_among")
  
  for(i in 1:length(unique(long$taxids))){
    taxoni<-unique(long$taxids)[i]
    longtemp<-long[long$taxids==taxoni,]
    newtab$taxon[i]=taxoni
    newtab$max_within[i]=max(longtemp[longtemp$variable==taxoni,"value"])
    newtab$min_within[i]=min(longtemp[longtemp$variable==taxoni,"value"])
    newtab$max_among[i]=max(longtemp[longtemp$variable!=taxoni,"value"])
    newtab$min_among[i]=min(longtemp[longtemp$variable!=taxoni,"value"])
  }

  newtab.list[[j]]<-newtab

  #plot
  newtablong<-reshape2::melt(newtab)
  newtablong$groups<-""
  newtablong$groups[newtablong$variable=="max_within" | newtablong$variable=="min_within"]<-paste("within",taxlevel[j],"range")
  newtablong$groups[newtablong$variable=="max_among" | newtablong$variable=="min_among"]<-paste("among",taxlevel[j],"range")
  
  plota<-ggplot(newtablong,aes(x = taxon,y = value,colour=groups,group=interaction(taxon,groups)))+geom_point()+
    theme(axis.text.x=element_text(size=8,angle=45, hjust=1),legend.title = element_blank(),
          panel.grid.major = element_blank(),panel.grid.minor = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+
    geom_line()+ylab("Percentage identity")+xlab("")
  
  if(!is.null(threshold)) plota<-plota+geom_hline(yintercept = threshold[j],linetype="dashed")
  
  plot.list[[j]]<-plota
}

out<-list(newtab.list,plot.list)

}