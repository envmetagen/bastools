add.3pmms<-function(ecopcroutput,Pf,Pr){
#add 3' mismatches to ecopcroutput
f_mismatch_table<-mismatch.table(ecopcroutput,Pf,"f")
f_mismatches_3prime<-as.data.frame(rowSums(f_mismatch_table[,as.integer(nchar(Pf)/2):nchar(Pf)]))
colnames(f_mismatches_3prime)<-"f_mismatches_3prime"
r_mismatch_table<-mismatch.table(ecopcroutput, Pr,"r")
r_mismatches_3prime<-as.data.frame(rowSums(r_mismatch_table[,as.integer(nchar(Pr)/2):nchar(Pr)]))
colnames(r_mismatches_3prime)<-"r_mismatches_3prime"
ecopcroutput<-cbind(ecopcroutput,f_mismatches_3prime,r_mismatches_3prime)

#add 3' mismatches to ecopcroutput - 6bp
f_mismatches_3prime6<-as.data.frame(rowSums(f_mismatch_table[,as.integer(nchar(Pf)-5):nchar(Pf)]))
colnames(f_mismatches_3prime6)<-"f_mismatches_3prime6"
r_mismatches_3prime6<-as.data.frame(rowSums(r_mismatch_table[,as.integer(nchar(Pr)-5):nchar(Pr)]))
colnames(r_mismatches_3prime6)<-"r_mismatches_3prime6"
ecopcroutput<-cbind(ecopcroutput,f_mismatches_3prime6,r_mismatches_3prime6)
}

add.res.ecopcroutput<-function(ecopcroutput){
#add taxonomic resolution to ecopcroutput
ecopcroutput.res<-add.res.Bas(ecopcroutput,obitaxdb=obitaxoR)

#remove extra columns
ecopcroutput.res$genus=NULL
ecopcroutput.res$genus_name=NULL
ecopcroutput.res$species_name=NULL
ecopcroutput.res$species=NULL
ecopcroutput.res$forward_tm=NULL
ecopcroutput.res$reverse_tm=NULL

ecopcroutput.res
}

add.tm.ecopcroutput<-function(ecopcroutput){
  ecopcroutput$fTms<-Tm.calc(ecopcroutput$forward_match)
  ecopcroutput$rTms<-Tm.calc(ecopcroutput$reverse_match)
  ecopcroutput$fTms3primehalf<-Tm.calc(substr(x = ecopcroutput$forward_match,
                                              start = as.integer(nchar(as.character(ecopcroutput$forward_match))/2+1),
                                              stop = nchar(as.character(ecopcroutput$forward_match))))
  ecopcroutput$rTms3primehalf<-Tm.calc(substr(x = ecopcroutput$reverse_match,
                                              start = as.integer(nchar(as.character(ecopcroutput$reverse_match))/2+1),
                                              stop = nchar(as.character(ecopcroutput$reverse_match))))
  ecopcroutput$fTms3prime6<-Tm.calc(substr(x = ecopcroutput$forward_match,
                                           start = as.integer(nchar(as.character(ecopcroutput$forward_match))-5),
                                           stop = nchar(as.character(ecopcroutput$forward_match))))
  ecopcroutput$rTms3prime6<-Tm.calc(substr(x = ecopcroutput$reverse_match,
                                           start = as.integer(nchar(as.character(ecopcroutput$reverse_match))-5),
                                           stop = nchar(as.character(ecopcroutput$reverse_match))))
  ecopcroutput
  
}

add.gc.ecopcroutput<-function(ecopcroutput){
  ecopcroutput$fgc<-gc.calc(ecopcroutput$forward_match)
  ecopcroutput$rgc<-gc.calc(ecopcroutput$reverse_match)
  ecopcroutput
}

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

#calculating gc%
#calculate gc% for each primer binding site in vector
gc.calc<-function(primerVec){
  gc<-vector()
  for (i in 1:length(primerVec)){
    gc[i]<-(stringr::str_count(primerVec[i],"G")+stringr::str_count(primerVec[i],"C"))/
      (nchar(primerVec[i]))*100
  }
  return(gc)
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



# #gc clamp present?
# gc.clamp<-function(primerVec){
#   gc<-vector()
#   for (i in 1:length(primerVec)){
#     gc[i]<-(stringr::str_count(primerVec[i],"G")+stringr::str_count(primerVec[i],"C"))/
#       (nchar(primerVec[i]))*100
#   }
#   return(gc)
# }