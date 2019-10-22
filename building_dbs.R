DL.nuccore.gb<-function(group.taxid,ncbiTaxDir,gene,rank_req="family"){
  
  children<-bas.get.children(group.taxid,ncbiTaxDir = ncbiTaxDir,rank_req = rank_req)
  
  message(paste("************",length(children$taxid), "children found at required rank ************"))
  
  pb = txtProgressBar(min = 0, max = length(children$taxid), initial = 0,style = 3)
  for(i in 1:length(children$taxid)){
    #could maybe add a trycatch to output remaining children df in case of failure?
    message(paste("Downloading data for",colnames(children)[3],":",children[i,3]))
    bas.get.nuccore(taxid = children[i,1],gene = gene,name = children[i,3])
    setTxtProgressBar(pb,i)
  }
  
  #remove files with empty result
  a<-system2("grep",args = c("'Empty result - nothing to do'",list.files(pattern = "*.gb")),wait = T,
             stderr = T,stdout = T)
  unlink(gsub(":\t.*","",a))
  a<-system2("grep",args = c("'Unable to obtain query'",list.files(pattern = "*.gb")),wait = T,
             stderr = T,stdout = T)
  unlink(gsub(":\t.*","",a))
}

extract.gene.gb<-function(gbfile,gene){
  #split gb file by record
  system2(command = "cat", args=c(gbfile, "|", "sh",paste0(bastoolsDir,"split_gb.sh")),
          wait=T) 
  #remove last file cause its always empty
  files<-list.files(pattern = "^outTemp.*")
  a<-suppressWarnings(max(as.numeric(do.call(rbind,stringr::str_split(files,"outTemp"))[,2])))
  unlink(paste0("outTemp",a))
  
  #extract gene for each file
  for(i in 1:length(list.files(pattern = "^outTemp.*"))){
    a<-list.files(pattern = "^outTemp.*")[i]
    if(gene=="18S") script<-paste0(bastoolsDir,"parse-genbank-18S.py")
    if(gene=="16S") script<-paste0(bastoolsDir,"parse-genbank-16S.py")
    if(gene=="COI") script<-paste0(bastoolsDir,"parse-genbank-COI.py")
    system2("python",args = c(script, a),wait = T,
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

bas.get.children<-function(group.taxid,ncbiTaxDir,rank_req=NULL){
  taxalist<-system2(command = "taxonkit",args = c("list","-r","--show-name", "--ids", group.taxid, "--indent", 
                                                  '""',"--data-dir",ncbiTaxDir),wait=T,stdout = T)
  taxalistdf<-data.table::fread(text = taxalist,sep = "\t",header = F)
  taxalistdf<-as.data.frame(do.call(rbind,stringr::str_split(taxalistdf$V1," \\[")))
  taxalistdf[,2:3]<-suppressWarnings(do.call(rbind,stringr::str_split(taxalistdf$V2,"\\] ")))
  
  if(!is.null(rank_req)){
  out<-taxalistdf[taxalistdf$V2==rank_req,]
  } else out <- taxalistdf
  
  colnames(out)<-c("taxid","rank",taxalistdf[1,3])
  return(out)
}


bas.get.nuccore<-function(taxid,name,gene){
  #######GENE CAN ONLY = 18S, 16S OR COI FOR NOW
  if(gene!="18S" & gene!="16S" & gene!="COI") stop("accepted values for gene are 16S, 18S or COI")
  
  if(gene=="18S") geneTerm<-
      "(18S ribosomal RNA[All Fields] OR 18S small subunit ribosomal RNA[All Fields] OR 18S*[Gene])"
  if(gene=="16S") geneTerm<-
      "(16S ribosomal RNA[All Fields] OR 16S small subunit ribosomal RNA[All Fields] OR 16S*[Gene])"
  if(gene=="COI") geneTerm<-
      "(COXi[Gene] OR CO1[Gene] OR COI[Gene] OR MTCO1[Gene] OR COX1[Gene]) AND (cytochrome oxidase subunit 1[All fields] OR cytochrome c oxidase subunit 1[All fields] OR cytochrome c oxidase subunit I[All fields] OR cytochrome oxidase subunit I[All fields])"
  if(gene=="12S") geneTerm<-
      "(12S ribosomal RNA[All Fields] OR 12S*[Gene])"
          #seems like 12S is a little different
  # rRNA            73..1024
  # /product="s-rRNA"
  # /note="12S ribosomal RNA"
  #OR
  # rRNA            <1..>170
  # /product="12S ribosomal RNA"
  out<-paste0(name,"_",taxid,"_",gene,".gb")
  
  if(!out %in% list.files()) {
  
  searchQ <- paste0("txid",taxid, "[Organism] AND ", geneTerm)
  search_results <- rentrez::entrez_search(db = "nuccore", term = searchQ, retmax = 9999999, use_history = T)
  
  if(length(search_results$ids)<600){
    DLseqs <- rentrez::entrez_fetch(db = "nuccore", web_history = search_results$web_history, rettype = "gb")
    
    writeLines(DLseqs,out)
  }
  
  if(length(search_results$ids)>=600){
    
    # #modified from primer miner:
    start <- 0
    chunks <- ceiling(length(search_results$ids)/100)
    
    for (i in 1:chunks) {
      DLseqs <- rentrez::entrez_fetch(db = "nuccore",
                                      web_history = search_results$web_history, rettype = "gb",
                                      retmax = 100, retstart = start)
      cat(DLseqs, file = out, sep = "",  append = T)
      start <- start + 100
      Sys.sleep(2.5)
    }
  }
  } else {message("Skipping download, file already exists")}
}

gb2fasta<-function(gbfile){
  #convert to fasta (including taxids)
  system2(command = "obiconvert", args=c("--genbank", gbfile, "--fasta-output"),wait = T,
          stdout=gsub(gbfile,pattern = ".gb",replacement = ".fasta"))
}



#' Download sequences from NCBI and BOLD based on taxonomy and gene search terms
#' @param groups Can be Families, Orders, higher?? As many groups as required can be entered
#' @param target_gene The gene to be downloaded. Make sure to include any variations of the gene (e.g. COI, COX1)
#'     as a string c("COI","COX1")
#' @param out name of output file
#' @param fix_names Make the sequence headers readable for \code{finalise_DB}
#' @return A fasta file of sequences found
#' @note there are no errors if a taxon name or target gene is not found, so please check that your taxa and genes are
#' spelled correctly for the databases
#' @examples
#' groups=c("Anura","vafs","faga")
#' subgroup1=c("Leiopelmatidae","Leiopelma")
#' subgroup2="Ursidae"
#' target_gene="COI"
#' out="results.fasta"
#' build.db(groups = groups,target_gene = target_gene,out=out)
#' @export
download.db<-function(groups,target_gene,BOLD=T,genbank=T,mito=F){

  origpath<-getwd()
  temppath<-paste0(getwd(),"/builddb-",as.numeric(Sys.time()))
  dir.create(path = temppath)

  #make folder for these downloads
  dlfolder<-paste0(temppath,"/downloaded_seqs")
  dir.create(path = dlfolder)
  setwd(dlfolder)

  ###################why is bold search limited to 2000bp?

#write taxa.csv
a1<-system.file("extdata", "config.txt", package = "bastools")
col2<-c(rep("NA",times=length(groups)))
df<-as.data.frame(cbind(groups,col2))
names(df)<-c("Order","Family")
write.csv(df,file = gsub(a1,pattern = "config.txt",replacement = "taxa.csv"),row.names = F)

#set config params
config.params<- readLines(a1)
config.params<-gsub(x=config.params, pattern = "\"TARGET_GENE\"",replacement =
  capture.output(cat(paste(shQuote(target_gene, type="cmd"), collapse=","))))

if(genbank==F){
config.params<-gsub(x=config.params,pattern = "download_GB = T",replacement = "download_GB = F")}
if(mito==F){
  config.params<-gsub(x=config.params,pattern = "download_mt = T",replacement = "download_mt = F")}

writeLines(config.params,con=gsub(a1,pattern = ".txt",replacement = "2.txt"))

#download sequences
#bold
if(BOLD==T){get_BOLD_BAS(groups)
}
#ncbi (short and long)
setwd(dlfolder)
batch_download.BAS(system.file("extdata", "taxa.csv", package = "bastools"),
               system.file("extdata", "config2.txt", package = "bastools"))
#back to orig folder
setwd(origpath)

#rm(add_mt,clipping_left_bold,clipping_left_GB,clipping_left_mt,clipping_rigth_bold,clipping_rigth_GB,
 #  clipping_rigth_mt,cmd,custom_query_GB,custom_query_mt,Download,download_bold,download_GB,download_mt,id,Marker,
  # maxlength_GB,maxlength_mt,Merge_and_Cluster_data,merge_bold,merge_GB,merge_mt,minlength_mt,no_marker,operating_system,
   #rm_dup,Skip_if_complete,Taxon_sep,Taxon_table,threshold,Version,vsearchpath)
}

format.downloads<-function(folder,out){

  ###grab downloaded file names
  origpath<-getwd()
  setwd(folder)
  a<-processx::run(command = "find", args=c("-type", "f"),echo=F,echo_cmd = F)
  b<-grep(pattern = "_BOLD",x = stringr::str_split(string = a$stdout,pattern = "\n")[[1]],value = T)
  d<-grep(pattern = "_GB",x = stringr::str_split(string = a$stdout,pattern = "\n")[[1]],value = T)
  g<-grep(pattern = "_mito",x = stringr::str_split(string = a$stdout,pattern = "\n")[[1]],value = T)

  cb <- function(line, proc) {cat(line, "\n")}

  #for mitochondrial downloads, convert from genbank to fasta and do cleanup
  if(length(g)>0){
  for(i in 1:length(g)){
    f<-process$new(command = "obiconvert", args=c("--genbank", g[i], "--fasta-output"),
                   echo_cmd = T,stdout=paste(gsub(g[i],pattern = ".gb",replacement = ".fasta")))
    f$wait()}

    a<-processx::run(command = "find", args=c("-type", "f"), echo=F,echo_cmd = F)
    l<-grep(pattern = "_mito.fasta",x = stringr::str_split(string = a$stdout,pattern = "\n")[[1]],value = T)
    for(i in 1:length(l)){
      f<-process$new(command = "sed", args=c("s/ /_/2;s/organism=/species=/", l[i]),
                     echo_cmd = T,stdout = gsub(".fasta",".clean.fasta",l[i]))
      f$wait()}
}

  #get filenames again
  a<-processx::run(command = "find", args=c("-type", "f"), echo=F,echo_cmd = F,stderr_line_callback = cb)
  l<-grep(pattern = "_mito.fasta",x = stringr::str_split(string = a$stdout,pattern = "\n")[[1]],value = T)

  #if bold downloaded, do cleanup on it
  if(length(b)>0){
    for(i in 1:length(b)){
      e<-process$new(command = "sed", args=c("s/_/ species=/;s/ /_/2;/^>/s/$/;/", b[i]),
                     echo_cmd = T,stdout = gsub(".fasta",".clean.fasta",b[i]))
      f$wait()}
  }

  #if NCBI downloaded, do cleanup on it
  if(length(d)>0){
    for(i in 1:length(d)){
      e<-process$new(command = "sed", args=c("s/ /_/2;s/ / species=/;s/ /; /2", d[i]),
                     echo_cmd = T,stdout = gsub(".fasta",".clean.fasta",d[i]))
      f$wait()}
  }

  #concatenate cleaned fastas
  a<-processx::run(command = "find", args=c("-type", "f"), echo=F,echo_cmd = F,stderr_line_callback = cb)
  b<-grep(pattern = "clean",x = stringr::str_split(string = a$stdout,pattern = "\n")[[1]],value = T)
  f<-process$new(command = "cat", args=c(b), echo_cmd = T,stdout="temp.fasta")
  f$wait()

  #remove extra newlines
  f<-process$new(command = "sed", args=c("/^[^>]/s/-//g", "temp.fasta"),echo_cmd = T,stdout="temp2.fasta")
  f$wait()
  #give unique ids
  f<-process$new(command = "obiannotate", args=c("--uniq-id","temp2.fasta"), echo_cmd = T,stdout=out)
  f$wait()
  #move
  f<-processx::run(command = "mv", args=c(out,origpath),echo=F,echo_cmd = T)
  #remove temp files, reset path
  unlink("temp.fasta")
  unlink("temp2.fasta")
  setwd(origpath)
}



