library(dplyr)

#run 8 epi2me
#https://epi2me.nanoporetech.com/workflow_instance/215579

mastersheet<-google.read.master("Mussels_datasheet")
expid<-"2019_September_001"
demuxed_fstq_dir<-"/media/sf_Documents/WORK/CIBIO/RUNS/MINION/run8/downloads/PASS/"
run="Run8"

bas.minion.pipeline1<-function(demuxed_fstq_dir,mastersheet,expid,run)
  