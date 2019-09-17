make.blastdb.bas<-function(infasta, makeblastdb_exec="makeblastdb"){
#GENERAL FORMATTING
message("assumes fasta with species=xxx;   should add option for taxid=xxx;")
  
  
tempfileA<-paste0("temp",as.numeric(Sys.time()),".fasta")
add.lineage.fasta.BAS(infasta,ncbiTaxDir,out = tempfileA)
tempfasta<-phylotools::read.fasta(tempfileA)

#remove seqs with no taxonomy 
if(length(tempfasta$seq.name[grep("kingdom=unknown",tempfasta$seq.name)])>0){
  message(paste0("Discarding ",length(tempfasta$seq.name[grep("kingdom=unknown",tempfasta$seq.name)]),
                 " sequences for which taxonomy could not be found"))
  tempfasta<-tempfasta[-grep("kingdom=unknown",tempfasta$seq.name),]
}

phylotools::dat2fasta(tempfasta,tempfileA)
#give unique ids
tempfileB<-paste0("temp",as.numeric(Sys.time()),".fasta")
system2(command = "obiannotate", args=c("--uniq-id",tempfileA), stdout=tempfileB,stderr = "",wait = T)
message("IDs made unique")

#Add "database=DB_MiFish12S_EUfwfish_04022019" to header
message("adding db name to headers")
sedarg=paste0("s/taxid=/database=",gsub(".fasta","",infasta),"; taxid=/g")
h<-process$new(command = "sed", args=c(sedarg, tempfileB), echo_cmd = T,
               stdout=gsub(".fasta",".blastdbformat.fasta",infasta))
h$wait()

system2(command = makeblastdb_exec, 
        args=c("-in", gsub(".fasta",".blastdbformat.fasta",infasta), 
               "-dbtype", "nucl", 
               "-blastdb_version", "5",
               "-parse_seqids","-out",
               gsub(".fasta","",infasta)), stderr = "",wait = T)
}
