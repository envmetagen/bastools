ecopcr2refs<-function(ecopcrfile,outfile){
  message("Reminder: assumes -D option was used for ecopcr")
  #read results
  ecopcroutput<-data.table::fread(ecopcrfile,sep = "\t")
  #remove hits outside desired lengths
  ecopcroutput<-ecopcroutput[!ecopcroutput$amplicon_length<min_length,]
  ecopcroutput<-ecopcroutput[!ecopcroutput$amplicon_length>max_length,]
  #remove duplicates (i.e. pick one entry per AC, based on lowest mismatches)
  ecopcroutput$total_mismatches<-as.numeric(ecopcroutput$forward_mismatch)+as.numeric(ecopcroutput$reverse_mismatch)
  ecopcroutput <- ecopcroutput[order(ecopcroutput$AC,ecopcroutput$total_mismatches),]
  ecopcroutput<-ecopcroutput[!duplicated(ecopcroutput$AC),]
  #remove duplicate sequences...for this database
  ecopcroutput<-ecopcroutput[!duplicated(ecopcroutput$sequence),] 
  #output as fasta
  colnames(ecopcroutput)<-gsub("sequence", "seq.text",colnames(ecopcroutput))
  ecopcroutput$seq.text<-toupper(ecopcroutput$seq.text)
  ecopcroutput$seq.name<-paste0(ecopcroutput$AC," ","taxid=",ecopcroutput$taxid,"; family=",
                                ecopcroutput$family_name, "; ",ecopcroutput$definition)
  phylotools::dat2fasta(ecopcroutput[,c("seq.name","seq.text")],outfile)
}