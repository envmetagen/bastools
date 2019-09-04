add.lineage.fasta.BAS<-function(infasta,ncbiTaxDir,out){
message("Reminders: headers must contain \"species=species_name;\"")
  message("this NEEDS to be rewritten based on what I have learned writing add.lineage.df")
formatted.exShort<-phylotools::read.fasta(infasta)
formatted.exShort$name<-gsub("_"," ",stringr::str_match(formatted.exShort$seq.name, "species=(.*?);")[,2])
formatted.exShort$name<-gsub(" sp.$","",formatted.exShort$name)
write.table(formatted.exShort$name,file = "names.txt",row.names = F,col.names = F,quote = F)
f<-process$new(command = "taxonkit", args = c("name2taxid","names.txt","--data-dir",ncbiTaxDir),
               echo_cmd = T,stdout = "taxids1.txt")
f$wait()
taxids1<-read.table("taxids1.txt",sep = "\t")
colnames(taxids1)<-gsub("V2","taxid1",colnames(taxids1))

#for unfound taxids, repeat search using only first word
b<-taxids1[is.na(taxids1$taxid1),]
b$V1<-gsub(" .*","",b$V1)
write.table(b$V1,file = "names2.txt",row.names = F,col.names = F,quote = F)
f<-process$new(command = "taxonkit", args = c("name2taxid","names2.txt","--data-dir",ncbiTaxDir),
               echo_cmd = T,stdout = "taxids2.txt")
f$wait()
taxids2<-read.table("taxids2.txt",sep = "\t")
colnames(taxids2)<-gsub("V2","taxid2",colnames(taxids2))

#####Note that if a name matches more than one taxid, taxonkit creates a new row and includes both taxids, so using merge
d<-taxids1[!duplicated(taxids1),] #just for speed during merge
#also, if a "name" has more than one taxid, remove
d2<-as.data.frame(table(d$V1))
d3<-d2[d2$Freq>1,] #taxa with more than one taxid
d4<-d[!d$V1 %in% d3$Var1,]
message(c("The following taxa were removed due to ambiguous taxids: ",paste(d3$Var1,collapse = " ")))
e<-merge(formatted.exShort,d4,by.x = "name",by.y = "V1",all.x = T,all.y = F)

#repeat for taxids2
d<-taxids2[!duplicated(taxids2),] #just for speed during merge
#also, if a "name" has more than one taxid, remove
d2<-as.data.frame(table(d$V1))
d3<-d2[d2$Freq>1,] #taxa with more than one taxid
d4<-d[!d$V1 %in% d3$Var1,]
message(c("The following taxa were removed due to ambiguous taxids: ",paste(d3$Var1,collapse = " ")))

e$name2<-gsub(" .*","",e$name)
f<-merge(e,d4,by.x = "name2",by.y = "V1",all.x = T,all.y = F)

#combine taxids into one column
f$taxid3<-do.call(pmax, c(f[,c("taxid1","taxid2")], list(na.rm=TRUE)))

#get taxonomy path
write.table(f$taxid3,file = "taxids3.txt",row.names = F,col.names = F,quote = F)
g<-process$new(command = "taxonkit", args = c("lineage","taxids3.txt","--data-dir",ncbiTaxDir),
               echo_cmd = T,stdout = "lineage.txt")
g$wait()
g<-process$new(command = "taxonkit", args = c("reformat","lineage.txt","--data-dir",ncbiTaxDir),
               echo_cmd = T,stdout = "lineage.reformat.txt")
g$wait()

lineage.reformat<-read.table("lineage.reformat.txt",sep = "\t")

unlink("names.txt")
unlink("names2.txt")
unlink("taxids1.txt")
unlink("taxids2.txt")
unlink("taxids3.txt")
unlink("lineage.txt")
unlink("lineage.reformat.txt")

taxpaths<-stringr::str_split(lineage.reformat$V3,pattern = ";")
taxpaths2<-as.data.frame(t(as.data.frame(taxpaths)))
rownames(taxpaths2)<-NULL
taxpaths2$full<-paste0("kingdom=",taxpaths2$V1,"; phylum=",taxpaths2$V2,"; class=",taxpaths2$V3,"; order=",
                       taxpaths2$V4,"; family=",taxpaths2$V5,"; genus=",taxpaths2$V6,"; species=",taxpaths2$V7,";")
taxpaths2$full<-gsub("=;","=unknown;",taxpaths2$full)

#make new headers: "seqid taxid path definition"
seqid<-gsub(" .*","",formatted.exShort$seq.name)
definition<-gsub("^(.*) species=.*;","",formatted.exShort$seq.name)
newname<-paste0(seqid," taxid=",f$taxid3,"; ",taxpaths2$full,definition)
#make df for outputting as fasta
formatted.exShort.winfo<-as.data.frame(cbind(newname,as.character(formatted.exShort$seq.text)))
colnames(formatted.exShort.winfo)<-c("seq.name","seq.text")
phylotools::dat2fasta(formatted.exShort.winfo,outfile = out)
}

add.lineage.df<-function(df,ncbiTaxDir){
  if(is.null(df$taxids)) {stop("No column called taxids")}
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
  df[,(length(df)-6):length(df)] <- sapply(df[,(length(df)-6):length(df)],as.character)
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
