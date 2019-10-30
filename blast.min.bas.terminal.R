
args = commandArgs(trailingOnly=TRUE)

message("usage: command fasta refdb basToolsDir")

if(length(args)>0){
  infastas<-args[1]
  refdb<-args[2]
  basToolsDir<-args[3]
}

source(paste0(basToolsDir,"blast.min.bas.R"))

h<-blast.min.bas(infastas,refdb,blast_exec="blastn",wait=T)

check.blasts(infastas,h)