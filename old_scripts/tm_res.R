#calculating Tm
#calculate Tm for each primer binding site in vector
Tm.calc<-function(primerVec){
  Tms<-vector()
for (i in 1:length(primerVec)){
  Tms[i]<-2*(stringr::str_count(primerVec[i],"A")+stringr::str_count(primerVec[i],"T")) +
    4*(stringr::str_count(primerVec[i],"G")+stringr::str_count(primerVec[i],"C"))
}
  return(Tms)
}

#add taxonomic resolution to ecopcroutput
add.res.Bas<-function(ecopcroutput,obitaxdb){
  if (class(obitaxdb)[1]!="obitools.taxonomy") obitaxdb2=ROBITaxonomy::read.taxonomy(obitaxdb)
  if (class(obitaxdb)[1]=="obitools.taxonomy") obitaxdb2=obitaxdb
  res=ROBIBarcodes::resolution(taxonomy = obitaxdb2,ecopcr = ecopcroutput)
  a<-as.data.frame(res)
  a$res<-gsub(x = a$res,pattern = "infraclass",replacement = "class")
  a$res<-gsub(x = a$res,pattern = "subfamily",replacement = "family")
  a$res<-gsub(x = a$res,pattern = "tribe",replacement = "family")
  a$res<-gsub(x = a$res,pattern = "suborder",replacement = "order")
  a$res<-gsub(x = a$res,pattern = "infraorder",replacement = "order")
  a$res<-gsub(x = a$res,pattern = "superfamily",replacement = "order")
  a$res<-gsub(x = a$res,pattern = "subclass",replacement = "class")
  a$res<-gsub(x = a$res,pattern = "subgenus",replacement = "genus")
  a$res<-gsub(x = a$res,pattern = "species subgroup",replacement = "species")
  a$res<-gsub(x = a$res,pattern = "cohort",replacement = "order")

  #K,P,C,O,F,G,S

  b<-cbind(ecopcroutput,a)
}

#function to calculate % species with family or lower res
res.fam.or.better<-function(ecopcroutput){
  (length(ecopcroutput$res[ecopcroutput$res=="family"])+
     length(ecopcroutput$res[ecopcroutput$res=="genus"])+
     length(ecopcroutput$res[ecopcroutput$res=="species"])) /
    (length(ecopcroutput$res))
}
