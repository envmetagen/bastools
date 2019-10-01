library(dplyr)
library(processx)

setwd("/home/bastian.egeter/git_bastools/bastools/")
source("bas.minion.pipeline1.R")
source("blast.min.bas.R")
source("googlesheet.foos.R")

googlesheets::gs_auth()

#run 8 epi2me
#https://epi2me.nanoporetech.com/workflow_instance/215579

mastersheet<-google.read.master("Mussels_datasheet")
expid<-"2019_September_001"
demuxed_fstq_dir<-"/media/sf_Documents/WORK/CIBIO/RUNS/MINION/run8/downloads/PASS/"
run="Run8"

bas.minion.pipeline1(emuxed_fstq_dir,mastersheet,expid,run)

setwd("/home/bastian.egeter/Minion_data/run8/")

infastas<-list.files(pattern = "test*")
refdb = "/mnt/Disk1/Tools/BLAST+/DBs/nt_v5/nt_v5"

######need to add exits status to this!

blast.min.bas(infastas,refdb,
              blast_exec = "/home/bastian.egeter/Tools/ncbi-blast-2.9.0+/bin/blastn",wait = F)  
  