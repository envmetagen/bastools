#infasta<-"DB_MiFish12S_EUfwfish_04022019.temp.fasta"
add.lineage.fasta.BAS<-function(infasta,ncbiTaxDir,out){
message("Reminders: headers must contain \"species=species_name;\"")
  
fasta.table<-phylotools::read.fasta(infasta)
fasta.table$name<-gsub("_"," ",stringr::str_match(fasta.table$seq.name, "species=(.*?);")[,2])
fasta.table$taxids<-names2taxids(vector = fasta.table$name,ncbiTaxDir)
new_lineage<-add.lineage.df(fasta.table,ncbiTaxDir)

new_lineage$full<-paste0("kingdom=",new_lineage$K,"; phylum=",new_lineage$P,"; class=",new_lineage$C,
                         "; order=",new_lineage$O,"; family=",new_lineage$F,"; genus=",new_lineage$G,
                         "; species=",new_lineage$S,";")
new_lineage$full<-gsub("=;","=unknown;",new_lineage$full)

#make new headers: "seqid taxid path definition"
seqid<-gsub(" .*","",fasta.table$seq.name)
definition<-gsub("^(.*) species=.*;","",fasta.table$seq.name)
newname<-paste0(seqid," taxid=",new_lineage$taxids,"; ",new_lineage$full,definition)
#make df for outputting as fasta
fasta.table.winfo<-as.data.frame(cbind(newname,as.character(fasta.table$seq.text)))
colnames(fasta.table.winfo)<-c("seq.name","seq.text")
phylotools::dat2fasta(fasta.table.winfo,outfile = out)
}

names2taxids<-function(vector,ncbiTaxDir){
  message("Reminder - this function may require user input, DO NOT RUN LINES OF SCRIPT AFTER THIS FUNCTION")
  names_fileA<-paste0("names",as.numeric(Sys.time()),".txt")
  vector<-as.character(vector)
  #dont search names with sp. in the first round
  if(length(grep(" sp\\.",vector))>0){
  vector2<-vector[-grep(" sp\\.",vector)]} else {vector2<-vector}
  write.table(unique(vector2),file = names_fileA,row.names = F,col.names = F,quote = F) 
  taxids_fileA<-gsub("names","taxids",names_fileA)
  system2(command = "taxonkit", args = c("name2taxid",names_fileA,"-r","--data-dir",ncbiTaxDir),
                 stdout = taxids_fileA,stderr = "",wait = T)
  taxidsA<-read.table(taxids_fileA,sep = "\t")

  #for unfound taxids, repeat search using only first word
  namesB<-c(as.character(taxidsA$V1)[is.na(taxidsA$V2)],vector[grep(" sp\\.",vector)])
  if(length(namesB)>0){
  namesB<-gsub(" .*","",namesB)
  names_fileB<-paste0("names",as.numeric(Sys.time()),".txt")
  taxids_fileB<-gsub("names","taxids",names_fileA)
  write.table(unique(namesB),file = names_fileB,row.names = F,col.names = F,quote = F)
  system2(command = "taxonkit", args = c("name2taxid",names_fileB,"-r","--data-dir",ncbiTaxDir),
          stdout = taxids_fileB,stderr = "",wait = T)
  taxidsB<-read.table(taxids_fileB,sep = "\t")}

  #Note that if a name matches more than one taxid, taxonkit creates a new row and includes both taxids, 
  #allow user input for choices
  message("Should add choices for unknowns")
  
  
  #first for good results in taxidsA
  taxidsA2<-taxidsA[!is.na(taxidsA$V2),]
  #find names with multiple matches
  if(length(taxidsA2$V1[duplicated(taxidsA2$V1)])>0){
  
  taxidsA4<-taxidsA2[duplicated(taxidsA2$V1),]
  taxidsA5<-taxidsA2[taxidsA2$V1 %in% taxidsA4$V1,]
  colnames(taxidsA5)<-c("name","taxids","rank")
  #add lineage to make choice easier
  taxidsA5<-add.lineage.df(taxidsA5,ncbiTaxDir)
  taxidsA5$old_taxids=NULL
  #present choice
  message("The following ", length(unique(taxidsA5$name))," taxa had more than one taxid match. Please
  type the full taxid of your choice")
  choicesA_df<-as.data.frame(unique(taxidsA5$name))
  for(i in 1:length(unique(taxidsA5$name))){
  print(taxidsA5[taxidsA5$name %in%  unique(taxidsA5$name)[i],])
  choicesA_df[i,2] <- readline("Type full chosen taxid: ")  
  }
  choicesA_df$V3<-"NA"
  colnames(choicesA_df)<-c("V1","V2","V3")
  #combine good results
  taxidsA6<-taxidsA2[!taxidsA2$V1 %in% unique(taxidsA5$name),]
  taxidsA7<-rbind(taxidsA6,choicesA_df)
  } else {taxidsA7<-taxidsA2}
  
  
  if(length(namesB)>0){
  #for good results in taxidsB
  taxidsB2<-taxidsB[!is.na(taxidsB$V2),]
  
  #find names with multiple matches
  if(length(taxidsB2$V1[duplicated(taxidsB2$V1)])>0){
  
  taxidsB4<-taxidsB2[duplicated(taxidsB2$V1),]
  taxidsB5<-taxidsB2[taxidsB2$V1 %in% taxidsB4$V1,]
  colnames(taxidsB5)<-c("name","taxids","rank")
  #add lineage to make choice easier
  taxidsB5<-add.lineage.df(taxidsB5,ncbiTaxDir)
  taxidsB5$old_taxids=NULL
  #present choice
  message("The following ", length(unique(taxidsB5$name))," taxa had more than one taxid match. Please
  type the full taxid of your choice")
  choicesB_df<-as.data.frame(unique(taxidsB5$name))
  for(i in 1:length(unique(taxidsB5$name))){
    print(taxidsB5[taxidsB5$name %in%  unique(taxidsB5$name)[i],])
    choicesB_df[i,2] <- readline("Type full chosen taxid: ")  
  }
  choicesB_df$V3<-"NA"
  colnames(choicesB_df)<-c("V1","V2","V3")
  #combine good results
  taxidsB6<-taxidsB2[!taxidsB2$V1 %in% unique(taxidsB5$name),]
  taxidsB7<-rbind(taxidsB6,choicesB_df)
  } else {taxidsB7<-taxidsB2}
  
  }
  #combine taxidA and taxidB results
  outdf<-as.data.frame(vector)
  outdf$vec2<-gsub(" sp\\.","",outdf$vector)
  outdf$vec3<-gsub(" .*","",outdf$vector)
  outdf<-merge(outdf,taxidsA7,by.x = "vec2",by.y = "V1",all.x = T,all.y = F)
  if(length(namesB)>0){
  outdf<-merge(outdf,taxidsB7,by.x = "vec3",by.y = "V1",all.x = T,all.y = F)
  outdf$V2.z<-ifelse(!is.na(outdf$V2.x) & !is.na(outdf$V2.y) | is.na(outdf$V2.y),outdf$V2.z<-outdf$V2.x,NA)}

  #combine taxids into one column
  if(length(namesB)>0){
  outdf$taxids<-do.call(pmax, c(outdf[,c("V2.z","V2.y")], list(na.rm=TRUE)))} else {outdf$taxids<-outdf$V2}
  
  #remove files
  unlink(names_fileA)
  unlink(taxids_fileA)
  if(length(namesB)>0){
  unlink(names_fileB)
  unlink(taxids_fileB)}
  
  #remove extraneous columns
  colnames(outdf)<-gsub("vector","name",colnames(outdf))
  outdf<-outdf[match(vector, outdf$name),]
  outdf<-outdf[,"taxids"]
}

add.lineage.df<-function(df,ncbiTaxDir){
  if(is.null(df$taxids)) {stop("No column called taxids")}
  df$taxids<-as.integer(df$taxids) 
  #write taxids to file
  taxids_fileA<-paste0("taxids",as.numeric(Sys.time()),".txt")
  write.table(unique(df$taxids),file = taxids_fileA,row.names = F,col.names = F,quote = F)
  
  #get taxonomy from taxids and format in 7 levels
  taxids_fileB<-paste0("taxids",as.numeric(Sys.time()),".txt")
  system2(command = "taxonkit",args =  c("lineage","-r",taxids_fileA,"-c","--data-dir",ncbiTaxDir)
          ,stdout = taxids_fileB,stderr = "",wait = T)
  taxids_fileC<-paste0("taxids",as.numeric(Sys.time()),".txt")
  system2(command = "taxonkit",args =  c("reformat",taxids_fileB,"-i",3,"--data-dir",ncbiTaxDir)
          ,stdout = taxids_fileC,stderr = "",wait = T)
  lineage<-as.data.frame(data.table::fread(file = taxids_fileC,sep = "\t"))
  colnames(lineage)<-gsub("V1","taxids",colnames(lineage))
  colnames(lineage)<-gsub("V2","new_taxids",colnames(lineage))
  colnames(lineage)<-gsub("V5","path",colnames(lineage))
  
  #merge with df
  message("replacing taxids with updated taxids. Saving old taxids in old_taxids.")
  df<-merge(df,lineage[,c(1,2,5)],by = "taxids")
  df$old_taxids<-df$taxids
  df$taxids<-df$new_taxids
  df$new_taxids=NULL
  df<-cbind(df,do.call(rbind, stringr::str_split(df$path,";")))
  colnames(df)[(length(df)-6):length(df)]<-c("K","P","C","O","F","G","S")
  df$K<-as.character(df$K)
  df$P<-as.character(df$P)
  df$C<-as.character(df$C)
  df$O<-as.character(df$O)
  df$F<-as.character(df$F)
  df$G<-as.character(df$G)
  df$S<-as.character(df$S)
  #not sure why the following command wasnt working
  #df[,(length(df)-6):length(df)] <- sapply(df[,(length(df)-6):length(df)],as.character)
  df[,(length(df)-6):length(df)][df[,(length(df)-6):length(df)]==""]<- "unknown"
  df$path=NULL
  unlink(taxids_fileA)
  unlink(taxids_fileB)
  unlink(taxids_fileC)
  return(df)
}


taxids2names<-function(df,ncbiTaxDir){
  taxids_fileA<-paste0("taxids",as.numeric(Sys.time()),".txt")
  write.table(unique(df$taxid),file = taxids_fileA,row.names = F,col.names = F,quote = F)
  taxids_fileB<-paste0("taxids",as.numeric(Sys.time()),".txt")
  g<-process$new(command = "taxonkit", args = c("lineage","-r",taxids_fileA,"--data-dir",ncbiTaxDir),
                 echo_cmd = T,stdout = taxids_fileB)
  g$wait()
  lineage<-read.table(taxids_fileB,sep = "\t")

  lineage$V3<-gsub("infraclass","class",lineage$V3)
  lineage$V3<-gsub("subfamily","family",lineage$V3)
  lineage$V3<-gsub("tribe","family",lineage$V3)
  lineage$V3<-gsub("suborder","order",lineage$V3)
  lineage$V3<-gsub("infraorder","order",lineage$V3)
  lineage$V3<-gsub("superfamily","order",lineage$V3)
  lineage$V3<-gsub("subclass","class",lineage$V3)
  lineage$V3<-gsub("subgenus","genus",lineage$V3)
  lineage$V3<-gsub("species subgroup","species",lineage$V3)
  lineage$V3<-gsub("cohort","order",lineage$V3)

  df2<-merge(df,lineage[,c(1,3)],by.x = "taxid",by.y = "V1", all.y = T)
  colnames(df)<-gsub("V3","rank",colnames(df))
  df<-cbind(df,do.call(rbind, stringr::str_split(df$path,";")))
  colnames(df)[(length(df)-6):length(df)]<-c("K","P","C","O","F","G","S")
  df[,(length(df)-6):length(df)] <- sapply(df[,(length(df)-6):length(df)],as.character)
  df[,(length(df)-6):length(df)][df[,(length(df)-6):length(df)]==""]<- "unknown"

  unlink(taxids_fileA)
  unlink(taxids_fileB)

  return(df)
}
