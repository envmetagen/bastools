ecopcr2refs<-function(ecopcrfile,outfile){
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
  #select one per family
  ecopcroutput<-ecopcroutput[!duplicated(ecopcroutput$family_name),] 
  
  #output as fasta
  colnames(ecopcroutput)<-gsub("sequence", "seq.text",colnames(ecopcroutput))
  ecopcroutput$seq.text<-toupper(ecopcroutput$seq.text)
  ecopcroutput$seq.name<-paste0(ecopcroutput$AC," ","taxid=",ecopcroutput$taxid,"; family=",
                                ecopcroutput$family_name, "; ",ecopcroutput$definition)
  phylotools::dat2fasta(ecopcroutput[,c("seq.name","seq.text")],outfile)
}

ecopcr2refs2<-function(ecopcrfile,outfile,bufferecopcr){
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
  
  #the edges removal in next step takes a long time, so choose one seq per genus
  ecopcroutput <- ecopcroutput[order(ecopcroutput$family_name,ecopcroutput$amplicon_length,decreasing = T),]
  ecopcroutput <-ecopcroutput[!duplicated(ecopcroutput[,c("family_name")]),] 
  
  #remove seqs with edges less than buffer
  a<-strsplit(ecopcroutput$sequence,split = "")
  ecopcroutput$leftbuffer<-"0"
  ecopcroutput$rightbuffer<-"0"
  
  d<-list()
  pb = txtProgressBar(min = 0, max = length(ecopcroutput$leftbuffer), initial = 0,style = 3) 
  for(i in 1:length(ecopcroutput$leftbuffer)){
    ecopcroutput$leftbuffer[i]<-match(FALSE,a[[i]] %in% letters)
    d[[i]]<-a[[i]][ecopcroutput$leftbuffer[i]:length(a[[i]])]
    ecopcroutput$rightbuffer[i]<-length(d[[i]])-as.numeric(match(TRUE,d[[i]] %in% letters))
    setTxtProgressBar(pb,i)
  }
  #remove seqs with left buffer less than buffer
  ecopcroutput<-ecopcroutput[!ecopcroutput$leftbuffer<bufferecopcr,]
  #remove seqs with right buffer less than buffer
  ecopcroutput<-ecopcroutput[!ecopcroutput$rightbuffer<bufferecopcr,]
  
  #output as fasta
  colnames(ecopcroutput)<-gsub("sequence", "seq.text",colnames(ecopcroutput))
  ecopcroutput$seq.text<-toupper(ecopcroutput$seq.text)
  ecopcroutput$seq.name<-paste0(ecopcroutput$AC," taxid=",ecopcroutput$taxid,"; ",ecopcroutput$definition)
  phylotools::dat2fasta(ecopcroutput[,c("seq.name","seq.text")],outfile)
}
