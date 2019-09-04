#' Download NCBI taxonomy
#' @title Download NCBI taxonomy.
#' @param path Full path to download NCBI taxonomy to. If using current directory use \code{path="here"}
#' @return An unpacked NCBI taxonomy
#' @note Any previous versions of the taxonomy located the same directory are overwritten
#' @examples
#' getNCBItaxonomy("here")
#' @export
getNCBItaxonomy<-function(path){
  if(path!="here"){setwd(path)}
  list.files()
  file.remove("taxdump.tar.gz")
  cb <- function(line, proc) {cat(line, "\n")}
  processx::run(command = "wget",
                args=c('ftp://ftp.ncbi.nlm.nih.gov://pub/taxonomy/taxdump.tar.gz'),echo=F,stderr_line_callback = cb)
  taxdump<-"taxdump.tar.gz"
  processx::run(command = "tar",
                args=c("-xzf",taxdump),echo=F,stderr_line_callback = cb)
  b<-"SUCCESS!"
}

#' Convert NCBI taxonomy to obitools taxonomy
#' @title Convert NCBI taxonomy to obitools taxonomy.
#' @param path Full path to directory containing the (unzipped) NCBI taxonomy files. If using current directory use \code{path="here"}
#' @return An obitools taxonomy
#' @note This takes a long time.
#' @examples
#' NCBI2obitaxonomy("here","obitax_26-4-19")
#' @export
NCBI2obitaxonomy<-function(path,out){
  if(path!="here"){setwd(path)}
  a<-getwd()
  cb <- function(line, proc) {cat(line, "\n")}
  processx::run(command = "obitaxonomy", args=c("-t",a,"-d",out),echo=F,stderr_line_callback = cb,echo_cmd = T)
  b<-"SUCCESS!"
}

#' Download megan accession 2 taxonomy map
#' @title Download megan accession 2 taxonomy map.
#' @param path Full path to download acc2tax to. If using current directory use \code{path="here"}
#' @return An unpacked acc2tax map
#' @note This file is c. 1gAny previous versions of the taxonomy located the same directory are overwritten??????
#' @examples
#' get_acc2tax_map("/media/sf_Documents/WORK/CIBIO/STATS_AND_CODE/TAXONOMIES/")
#' @export
get_acc2tax_map<-function(path){
  if(path!="here"){setwd(path)}
  cb <- function(line, proc) {cat(line, "\n")}
  a<-processx::run(command = "wget",
                 args=c("-nc","-r","-nd","--no-parent","-A",'nucl_acc2tax*',"--spider",
                        'http://ab.inf.uni-tuebingen.de/data/software/megan6/download/'),echo=F,spinner = T)
  e<-stringr::str_extract_all(string = a$stderr,pattern =  "nucl_acc2tax(.*?)zip")[[1]][1]
  message(c("Found file from megan6/download/: ",e))
  f<-strsplit(e,split = ".",fixed = T)[[1]][1]
  print(list.files())
  if(length(grep(x = list.files(),pattern = f))>0){
    stop("Extracted file already exists in folder, not downloading",call. = F)}

  message("downloading...will take some time (>1Gb)")
  d<-processx::run(command = "wget",
                   args=c("-nc","-r","-nd","--no-parent","-A",'nucl_acc2tax*',
                          'http://ab.inf.uni-tuebingen.de/data/software/megan6/download/'),echo=F,spinner = T)
  message("unpacking...")
  processx::run(command = "unzip", args=c(e),echo=F,stderr_line_callback = cb)
  file.remove(e)
  b<-"SUCCESS!"
}


