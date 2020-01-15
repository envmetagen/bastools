make.blastdb.bas<-function(infasta,makeblastdb_exec="makeblastdb",addtaxidsfasta=T, ncbiTaxDir, dbversion=5){
  library(processx)
  
  if(!infasta %in% list.files()) stop("infasta not found in current directory")
  
  message("Reminder: Only works for blast 2.9.0: Current version:")
  system2(command = makeblastdb_exec,args = c("-version"))
  
  tempfileA<-paste0("temp",as.numeric(Sys.time()),".fasta")
  
  if(addtaxidsfasta==T){  
    message("assumes fasta with species=xxx;") 
    add.lineage.fasta.BAS(infasta,ncbiTaxDir,out = tempfileA)
    tempfasta<-phylotools::read.fasta(tempfileA)
    
    #remove seqs with no taxonomy 
    if(length(tempfasta$seq.name[grep("kingdom=unknown",tempfasta$seq.name)])>0){
      message(paste0("Discarding ",length(tempfasta$seq.name[grep("kingdom=unknown",tempfasta$seq.name)]),
                     " sequences for which taxonomy could not be found"))
      tempfasta<-tempfasta[-grep("kingdom=unknown",tempfasta$seq.name),]
    }
    phylotools::dat2fasta(tempfasta,tempfileA)
    
  } else{message("Assumes fasta with taxid=xxx;")
         file.copy(infasta,tempfileA)}
  
  #give unique ids
  message("giving unique ids")
  tempfileB<-paste0("temp",as.numeric(Sys.time()),".fasta")
  system2(command = "obiannotate", args=c("--uniq-id",tempfileA), stdout=tempfileB,stderr = "",wait = T)
  
  #Add "database name to header
  message("adding db name to headers")
  sedarg=paste0("s/taxid=/database=",gsub(".fasta","",infasta),"; taxid=/g")
  h<-process$new(command = "sed", args=c(sedarg, tempfileB), echo_cmd = T,
                 stdout=gsub(".fasta",".blastdbformatted.fasta",infasta))
  h$wait()
  
  #ensure ids are <50 characters
  tempfasta<-phylotools::read.fasta(gsub(".fasta",".blastdbformatted.fasta",infasta))
  defs<-suppressWarnings(do.call(rbind,strsplit(as.character(tempfasta$seq.name)," ")))
  defs[,1]<-stringr::str_trunc(as.character(defs[,1]),width = 49)
  #remove any quotes
  defs[,1]<-gsub('"',"",defs[,1])
  tempfasta$seq.name <- apply(as.data.frame(defs),1,paste,collapse = " " )
  phylotools::dat2fasta(tempfasta,gsub(".fasta",".blastdbformatted.fasta",infasta))
  
  #make mapping file
  taxids<-sub(x = stringr::str_match(string = tempfasta$seq.name,pattern = "taxid=(.*)"),pattern = ";.*",
              replacement = "")[,2]
  
  taxonmap<-as.data.frame(cbind(defs[,1],taxids))
  mappingfile<-paste0("mapping",as.numeric(Sys.time()),".txt")
  write.table(taxonmap,mappingfile,quote = F,sep = " ",row.names = F,col.names = F)
  
  system2(command = makeblastdb_exec, 
          args=c("-in", gsub(".fasta",".blastdbformatted.fasta",infasta), 
                 "-dbtype", "nucl", 
                 "-blastdb_version", dbversion,
                 "-parse_seqids","-taxid_map",mappingfile,"-out",
                 gsub(".fasta","",infasta)), stderr = "",wait = T)
  
  #testblast
  message("Running test blast")
  phylotools::dat2fasta(head(tempfasta,n=1),gsub(".fasta",".blastdbformatted.test.fasta",infasta))
  h<-blast.min.bas(gsub(".fasta",".blastdbformatted.test.fasta",infasta),refdb = gsub(".fasta","",infasta))
  check.blasts(gsub(".fasta",".blastdbformatted.test.fasta",infasta),h)
  message("Does new db have taxids?")
  print(data.table::fread(gsub(".fasta",".blastdbformatted.test.blast.txt",infasta)))
  system2(command = "blastdbcheck",args = c("-must_have_taxids","-db",gsub(".fasta","",infasta)))
  
  unlink(tempfileA)
  unlink(tempfileB)
  unlink(gsub(".fasta",".blastdbformatted.test.blast.txt",infasta))
  unlink(gsub(".fasta",".blastdbformatted.test.fasta",infasta))
}
