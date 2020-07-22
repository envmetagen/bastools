threshold.bin.blast2<-function(df,qseqidcol="origseqid",qseqcol="origseq",TaxlevelTest="G",taxidcol="origtaxid",qpath="origpath"){
  #TaxlevelTest can be G or F
  if(!qseqidcol %in% colnames(df)) stop(qseqidcol," missing")
  if(!qseqcol %in% colnames(df)) stop(qseqcol," missing")
  if(!taxidcol %in% colnames(df)) stop(taxidcol," missing")
  if(!qpath %in% colnames(df)) stop(qpath," missing")
  
  if(TaxlevelTest!="F" & TaxlevelTest!="G") stop("TaxlevelTest must be F or G")

  colnames(df)<-gsub(qseqidcol,"qseqid",colnames(df))
  colnames(df)<-gsub(qseqcol,"seq.text",colnames(df))
  colnames(df)<-gsub(taxidcol,"taxids",colnames(df))
  colnames(df)<-gsub(qpath,"qpath",colnames(df))
  
  df$seq.name<-paste0(df$qseqid," taxid=",df$taxids,";")
  
  lineage<-as.data.frame(stringr::str_split(df$qpath,";",simplify = T))
  colnames(lineage)<-c("K","P","C","O","F","G","S")
  
  df<-cbind(df,lineage)
  
  if(TaxlevelTest=="G") {
    ex.seqid.group<-"S"
    out<-paste0("tempBLASTDB.",TaxlevelTest,".tsv")
    if(file.exists(out)) file.remove(out)
  }
  
  if(TaxlevelTest=="F") {
    ex.seqid.group<-"G"
    out<-paste0("tempBLASTDB.",TaxlevelTest,".tsv")
    file.remove(out)
    if(file.exists(out)) file.remove(out)
  }
  
  #write blastdb
  db<-df[!duplicated(df$qseqid),]
  phylotools::dat2fasta(db[,c("seq.name","seq.text")],"temp.db.fasta")
  #make blastdb using metabinkit
  system2("metabinkit_blastgendb",c("-f","temp.db.fasta","-o", "temp.db","-c"),wait = T)
a<-Sys.time()
  for(i in 1:length(unique(df$qseqid))){
    
    message("loop ",i)
    
    #make query fasta
    b<-df[df$qseqid==unique(df$qseqid)[i],]
    b<-b[1,]
    phylotools::dat2fasta(b[,c("seq.name","seq.text")],"temp.seq.fasta")
    
    #get seqids of query group
    ex.seqids<-unique(df[df[,ex.seqid.group]==b[,ex.seqid.group],"qseqid"])  
    write.table(ex.seqids,"ex.seqids.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = F)
    
    #blast (cant use metabin_blast because of -negative_seqidlist)
    # system2(command = "blastdb_aliastool",
    #         args=c("-seqid_file_in", "ex.seqids.txt","-seqid_file_out","ex.seqids.out.txt"), wait = T)
    # 
    system2(command = "blastn",
            args=c("-task", "megablast", "-query", "temp.seq.fasta", "-db","temp.db","-outfmt",
                   "'6 qseqid saccver ssciname evalue staxid pident qcovs'","-evalue",1,"-num_threads", 16, "-max_target_seqs", 
                   100, "-max_hsps",1,"-word_size", 20,"-perc_identity", 70,"-qcov_hsp_perc",98,
                   "-gapopen", 0, "-gapextend", 2, "-reward", 1, "-penalty", -1, "-dust","no", 
                 #  "-negative_seqidlist", "ex.seqids.out.txt", 
                   "-out",
                   "temp.seq.blast.txt"), wait = T)
    #could exclude evalue and qcovs, might save a bit of time
    
    #store blast results
    lblast<-system2("wc",c("-l","temp.seq.blast.txt"),wait=T)
    
    #if(lblast>0) {
    blastResults<-data.table::fread("temp.seq.blast.txt",sep = "\t",data.table = F)
    write.table(blastResults,file = out,append = T,quote = F,row.names = F,sep = "\t",col.names = F)
  }

b<-Sys.time()
  
  #hard coded for now, note the "taxids" rather than "staxid", which metabin does not accept
  headers<-paste0("'1i",paste(c("qseqid", "saccver", "ssciname","evalue", "taxids", "pident", "qcovs"),collapse = "\t"),"'")
  
  system2("sed",c("-i", headers, out),wait = T)
  
  message("Ouput save to ",out)
}
