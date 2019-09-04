#make level of interest hits
ecohitsDerepByLevel<-function(ecopcrhits,levelofinterest){
ecopcrhits$full.seq<-paste0(as.character(ecopcrhits$forward_match),as.character(ecopcrhits$sequence),
                           as.character(ecopcrhits$reverse_match))

#split by levelofinterest
levelofinterest2<-paste0(levelofinterest,"_name")
splitdf<-split(ecopcrhits, f = ecopcrhits[,levelofinterest2])

#derep
splitdf2<-list()
for(i in 1:length(splitdf)){
  splitdf2[[i]]<-splitdf[[i]][!duplicated(splitdf[[i]]$full.seq),]
}

#back into df
ecopcrhits2<-as.data.frame(bind_rows(splitdf2))

#remove full.seq
ecopcrhits2$full.seq=NULL

return(ecopcrhits2)
}