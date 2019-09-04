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



