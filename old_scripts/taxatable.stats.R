#modelling taxtabs
binarise.taxatab<-function(taxatab){
  #transpose to have species as columns
  taxatab2<-as.data.frame(t(taxatab[,-1]))
  colnames(taxatab2)<-taxatab[,1]
  #make binary
  taxatab2[taxatab2>0]<-1
  return(taxatab2)
}

#make bray distance matrix
taxatab2bray<-function(binarised.taxatab){
  #convert binary matrix to distance matrix
  distance_matrix<-vegan::vegdist(binarised.taxatab,method = "bray",binary = T)
}