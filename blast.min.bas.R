blast.min.bas<-function(infastas,refdb,blast_exec="blastn",wait=T,outpattern=NULL){
  if(is.null(outpattern)){
  for(i in 1:length(infastas)){
  system2(command = blast_exec, args=c("-query", infastas[i],
          "-task", "megablast","-db",refdb,"-outfmt", '"6 qseqid evalue staxid pident qcovs"', 
          "-num_threads", "16", "-max_target_seqs", "100", "-max_hsps","1", "-out",
          paste0(gsub(x = infastas[i],pattern = ".fasta",replacement = ".blast.txt"))),stdout = NULL,wait = wait)
  }}
  if(!is.null(outpattern)){
    for(i in 1:length(infastas)){
      system2(command = blast_exec, args=c("-query", infastas[i],
               "-task", "megablast","-db",refdb,"-outfmt", '"6 qseqid evalue staxid pident qcovs"', 
              "-num_threads", "16", "-max_target_seqs", "100", "-max_hsps","1", "-out",
              paste0(gsub(x = infastas[i],pattern = ".fasta",replacement = outpattern))),stdout = NULL,wait = wait)
    }}
}

split.blast.min.bas<-function(infastas,refdb,blast_exec="blastn",outpattern=NULL){
  
  origDir<-getwd()
  tempDir<-(paste0(origDir,"/",as.numeric(Sys.time())))
  dir.create(tempDir)
  setwd(tempDir)
  
  for(i in 1:length(infastas)){
    
    #copy query
    file.copy(paste0(origDir,"/",infastas[i]),paste0(tempDir,"/",infastas[i]))
    
    #split fasta
    system2(command = "pyfasta", args=c("split","-n","10",infastas[i]), wait = T,stderr = "")
  }
    
    #list files
    queries<-list.files(pattern = ".0[0-9].fasta")
    
    h<-list()
    for(i in 1:length(queries)){
      h[[i]]<-process$new(command = "blastn", args=c("-query", queries[i], "-task", "megablast","-db",refdb,"-outfmt",
            "6 qseqid evalue staxid pident qcovs",
           "-num_threads", "2", "-max_target_seqs", "100", "-max_hsps","1", "-out",
           paste0(gsub(x = queries[i],pattern = ".fasta",replacement = ".blast.txt"))),
           echo_cmd = T)
    }
    
    for(i in 1:length(h)){
      h[[i]]$wait()
    }
    

    #concatenate blast results
for(i in 1:length(infastas)){
  
    catblasts<-grep(".blast.txt",list.files(pattern = gsub(".fasta","",infastas[i])),value = T)
    system2(command = "cat", args=c(catblasts), 
            stdout=gsub(".fasta",".blast.txt",infastas[i]), wait=T,stderr = "")
    
    file.copy(paste0(tempDir,"/",gsub(".fasta",".blast.txt",infastas[i])),
              paste0(origDir,"/",gsub(".fasta",outpattern,infastas[i])))
    
}

unlink(tempDir)
}


  
  
  