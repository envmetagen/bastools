blast.min.bas<-function(infastas,refdb,blast_exec="blastn",wait=T,taxidlimit=NULL,taxidname=NULL,ncbiTaxDir=NULL){
  
  if(!is.null(taxidlimit)) if(is.null(ncbiTaxDir)) stop("to use taxidlimit, ncbiTaxDir must be supplied")
  if(!is.null(taxidlimit)) if(is.null(taxidname)) stop("to use taxidlimit, taxidname must be supplied")
  if(!is.null(taxidlimit)) message("Make sure infastas,taxidlimit & taxidname are in correct order")
  
  t1<-Sys.time()
  
  library(processx)
  
  if(length(infastas)==1 | length(infastas)==2 | length(infastas)==3) threads<-8
  if(length(infastas)==4 | length(infastas)==5 | length(infastas)==6) threads<-4
  if(length(infastas)==7 | length(infastas)==8 | length(infastas)==9) threads<-2
  if(length(infastas)>9) threads<-1
  
  continue<-data.frame("file"<-infastas,"response"="y")
  continue$response<-as.character(continue$response)
  for(i in 1:length(infastas)){
    if(paste0(gsub(x = infastas[i],pattern = ".fasta",replacement = ".blast.txt")) %in% list.files()){
      continue[i,2]<-readline(paste0("The following file already exists, Overwrite? (y/n):", "
                                     ",gsub(x = infastas[i],pattern = ".fasta",replacement = ".blast.txt")))
    }
    }
  
  if("n" %in% continue$response) stop("Abandoned blast due to overwrite conflict")
  
  if(!is.null(taxidlimit)){
    
    h<-list()
    
    for(i in 1:length(infastas)){
      if(length(list.files(pattern = paste0(taxidname[i],"_taxidlimit.txt")))==0){
      
      system2(command = "taxonkit",args = c("list", "--ids", taxidlimit[i], "--indent", '""',"--data-dir",ncbiTaxDir)
              ,wait=T,stdout = paste0(taxidname[i],"_taxidlimit.temp.txt"))
      
      #remove blank row
      taxidlist<-read.table(paste0(taxidname[i],"_taxidlimit.temp.txt"))
      write.table(taxidlist,paste0(taxidname[i],"_taxidlimit.txt"),row.names = F,quote = F,col.names = F)
      
      unlink(paste0(taxidname[i],"_taxidlimit.temp.txt"))
      message(paste("taxidlist saved to",paste0(taxidname[i],"_taxidlimit.txt")))
      }else{ message(paste("The file",paste0(taxidname[i],"_taxidlimit.txt"),"already exists. Using that file."))}
      
      h[[i]]<-process$new(command = blast_exec, 
                          args=c("-query", infastas[i], "-task", "megablast","-db",refdb,"-outfmt",
                                 "6 qseqid evalue staxid pident qcovs","-num_threads", threads, "-taxidlist", 
                                 paste0(taxidname[i],"_taxidlimit.txt"),"-max_target_seqs", "100", "-max_hsps","1", "-out",
                                 paste0(gsub(x = infastas[i],pattern = "\\.fasta",replacement = ".blast.txt"))),echo_cmd = T,
                                 stderr = paste0("blast.error.temp.processx.file",i))
    }
  }
  
  if(is.null(taxidlimit)){  
    
    h<-list()
    
    for(i in 1:length(infastas)){
      
      h[[i]]<-process$new(command = blast_exec,
                          args=c("-query", infastas[i], "-task", "megablast","-db",refdb,"-outfmt",
                                 "6 qseqid evalue staxid pident qcovs","-num_threads", threads, "-max_target_seqs", 
                                 "100", "-max_hsps","1", "-out",
                                 paste0(gsub(x = infastas[i],pattern = "\\.fasta",replacement = ".blast.txt"))),
                                 echo_cmd = T,stderr = paste0("blast.error.temp.processx.file",i))
    }
  }
  
  Sys.sleep(time = 2)
  
  exits<-list()
  for(i in 1:length(h)){
    exits[[i]]<-h[[i]]$get_exit_status()
  }
  
  if(1 %in% exits){
    message("
            ************
            There was a problem with ", infastas[match(1,exits)], ", aborting all blasts
            ************")
    for(i in 1:length(h)){
      h[[i]]$kill()
    }
  }
  
  if(wait==T){
    for(i in 1:length(h)){
      h[[i]]$wait()
      message(paste(infastas[i],"blast pid",h[[i]]$get_pid()))
      message(readLines(paste0("blast.error.temp.processx.file",i)))
      unlink(paste0("blast.error.temp.processx.file",i))
    }
  }
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  
  message(c("All blasts complete in ",t3," mins."))
  
  return(h)
}

check.blasts<-function(infastas,h){
  a<-list()
  for(i in 1:length(h)){
    message(infastas[i])
    a[[i]]<-tryCatch(h[[i]]$get_status(),error=function(e) print("not running"))
    if(nchar(a[[i]])!=11) message(gsub("sleeping","running",h[[i]]$get_status()))
    message(paste0("exit status: ",h[[i]]$get_exit_status()))
  }
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


  
  
  