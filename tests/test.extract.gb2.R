extract.gene.gb<-function(gbfile,gene,bastoolsDir){
  #split gb file by record
  system2(command = "cat", args=c(gbfile, "|", "sh",paste0(bastoolsDir,"split_gb.sh")),
          wait=T) 
  #remove last file cause its always empty
  files<-list.files(pattern = "^outTemp.*")
  a<-suppressWarnings(max(as.numeric(do.call(rbind,stringr::str_split(files,"outTemp"))[,2])))
  unlink(paste0("outTemp",a))
  
  script<-readLines(paste0(bastoolsDir,"parse-genbank-source.py"))
  
  if(gene=="18S") {
    
    features<-c("rRNA","misc")
    terms<-c("18S","small subunit ribosomal RNA")
    qualifiers<-c("product","note")
    
    scriptlist<-list()
    start=1
    for(i in 1:length(features)){
      for(j in 1:length(terms)){
        for(k in 1:length(qualifiers)){
          scriptlist[[start]]<-script
          scriptlist[[start]]<-gsub("PUT_FEAT_TYPE_HERE",features[i],scriptlist[[start]])
          scriptlist[[start]]<-gsub("PUT_TERM_HERE",terms[j],scriptlist[[start]])
          scriptlist[[start]]<-gsub("PUT_QUALIFIER_HERE",qualifiers[k],scriptlist[[start]])
          start<-start+1
        }
      }
    }
  }
   for(i in 1:length(scriptlist)){ 
    writeLines(scriptlist[[i]],paste0(bastoolsDir,"parse-genbank_",i))
   }
    
   
  #extract gene for each file
  for(i in 1:length(list.files(pattern = "^outTemp.*"))){
    a<-list.files(pattern = "^outTemp.*")[i]
    
    system2("python",args = c(paste0(bastoolsDir,"parse-genbank-Temp.py"), a),wait = T,
            stdout = gsub("outTemp","extract.outTemp",a),stderr = F)
    
    # #some files have /note instead of /product
    # count<-system2("wc",args = c("-l",gsub("outTemp","extract.outTemp",a)),wait = T,stdout = T)
    # if(as.numeric(do.call(rbind,stringr::str_split(count," "))[,1])==0){
    #   system2("python",args = c(gsub(".py","_note.py",script), a),wait = T,
    #           stdout = gsub("outTemp","extract.outTemp",a),stderr = F)
    # }
    # #some files have "small subunit ribosomal RNA" in /product (no "16S")
    # count<-system2("wc",args = c("-l",gsub("outTemp","extract.outTemp",a)),wait = T,stdout = T)
    # if(as.numeric(do.call(rbind,stringr::str_split(count," "))[,1])==0){
    #   system2("python",args = c(gsub(".py","_ssrrna.py",script), a),wait = T,
    #           stdout = gsub("outTemp","extract.outTemp",a),stderr = F)
    # }
    
    #some files fail for other reasons 
    count<-system2("wc",args = c("-l",gsub("outTemp","extract.outTemp",a)),wait = T,stdout = T)
    if(as.numeric(do.call(rbind,stringr::str_split(count," "))[,1])==0){
      message(paste("One record from",gbfile,"failed and was excluded"))
      unlink(gsub("outTemp","extract.outTemp",a))
    }
  }
  
  #cat files
  if(length(list.files(pattern = "^extract.outTemp.*"))!=0){
    
    system2("cat",args=c(list.files(pattern = "^extract.outTemp.*")),
            stdout = gsub(".gb",".extract.fasta",gbfile),wait = T)
    
    #remove extraneuos files
    unlink(list.files(pattern = "^extract.outTemp.*"))
    unlink(list.files(pattern = "^outTemp.*"))
  }
}