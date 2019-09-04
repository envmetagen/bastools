build.db2.subsp<-function(out){
  message("not tested for BOLD")

  ###grab downloaded file names
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
      f<-process$new(command = "sed", args=c("-E", "s/;\\s/;/g;s/^(>.* organism=)([^ ;]+) /\\1\\2_/;s/^(>.* organism=)([^ ;]+) /\\1\\2_/;s/;/; /g;s/organism=/species=/", l[i]),
                     echo_cmd = T,stdout = gsub(x = l[i],pattern = "_mito",replacement = "_mito2"))
      f$wait()
      f<-process$new(command = "mv", args=c(gsub(x = l[i],pattern = "_mito",replacement = "_mito2"), l[i]),echo_cmd = T)
      f$wait()
      }
  }


  #get filenames again
  a<-processx::run(command = "find", args=c("-type", "f"), echo=F,echo_cmd = F,stderr_line_callback = cb)
  l<-grep(pattern = "_mito.fasta",x = stringr::str_split(string = a$stdout,pattern = "\n")[[1]],value = T)

  #if bold downloaded, do cleanup on it
 # if(length(b)>0){
  #  for(i in 1:length(b)){
   #   e<-processx::run(command = "sed", args=c("-i", "s/_/ species=/", b[i]),
  #                     echo=F,echo_cmd = T,stderr_line_callback = cb)
#    }
#    for(i in 1:length(b)){
#      e<-processx::run(command = "sed", args=c("-i", "s/ /_/2", b[i]),echo=F,echo_cmd = T,stderr_line_callback = cb)
#    }
#    for(i in 1:length(b)){
#      e<-processx::run(command = "sed", args=c("-i", "/^>/s/$/;/", b[i]),echo=F,echo_cmd = T)
#    }
#  }

  #if NCBI downloaded, do cleanup on it
  if(length(d)>0){
    for(i in 1:length(d)){
      e<-processx::run(command = "sed", args=c("-i", "s/ /_/2", d[i]),
                       echo=F,echo_cmd = T,stderr_line_callback = cb)
    }
    for(i in 1:length(d)){
      e<-processx::run(command = "sed", args=c("-i", "s/ / species=/", d[i]),
                       echo=F,echo_cmd = T,stderr_line_callback = cb)
    }
    for(i in 1:length(d)){
      e<-processx::run(command = "sed", args=c("-i", "s/ /; /2", d[i]),
                       echo=F,echo_cmd = T,stderr_line_callback = cb)
    }
  }

  #concatenate cleaned fastas
  if(length(b)>0 & length(d)>0 & length(l)>0) {f<-process$new(command = "cat", args=c(b,d,l),
                                                              echo_cmd = T,stdout=out)
  f$wait()}
  if(length(b)>0 & length(d)>0 & length(l)==0) {f<-process$new(command = "cat", args=c(b,d),
                                                               echo_cmd = T,stdout=out)
  f$wait()}
  if(length(b)==0 & length(d)>0 & length(l)>0) {f<-process$new(command = "cat", args=c(d,l),
                                                               echo_cmd = T,stdout=out)
  f$wait()}
  if(length(b)>0 & length(d)==0 & length(l)>0) {f<-process$new(command = "cat", args=c(b,l),
                                                               echo_cmd = T,stdout=out)
  f$wait()}
  if(length(b)>0 & length(d)==0 & length(l)==0) {f<-process$new(command = "cat", args=c(b),
                                                                echo_cmd = T,stdout=out)
  f$wait()}
  if(length(b)==0 & length(d)>0 & length(l)==0) {f<-process$new(command = "cat", args=c(d),
                                                                echo_cmd = T,stdout=out)
  f$wait()}
  if(length(b)==0 & length(d)==0 & length(l)>0) {f<-process$new(command = "cat", args=c(l),
                                                                echo_cmd = T,stdout=out)
  f$wait()}

  #remove extra newlines
  e<-processx::run(command = "sed", args=c("-i", "/^[^>]/s/-//g", out),
                   echo=F,echo_cmd = T)

}
