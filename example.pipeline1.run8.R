library(dplyr)
library(processx)

setwd("/home/bastian.egeter/git_bastools/bastools/")
source("bas.minion.pipeline1.R")
source("blast.min.bas.R")
source("googlesheet.foos.R")
source("bin.blast.R")
source("add.taxids.fasta.BAS.R")
source("merge_MBC_otutab_with_bin_blast.R")
source("taxatab.filter.R")



#run 8 epi2me
#https://epi2me.nanoporetech.com/workflow_instance/215579
googlesheets::gs_auth()
mastersheet<-google.read.master("Mussels_datasheet")
write.table(mastersheet,sep = "\t",quote = F,row.names = F,file = "Mussels_datasheet_master.txt")
mastersheetfile<-"Mussels_datasheet_master.txt"
expid<-"2019_September_001"
demuxed_fstq_dir<-"/media/sf_Documents/WORK/CIBIO/RUNS/MINION/run8/downloads/PASS/"
run="Run8"
ncbiTaxDir<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/taxdump/"
obitaxdb<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/obitaxdump/obitaxdb"
blast_exec<-"/home/bastian.egeter/Tools/ncbi-blast-2.9.0+/bin/blastn"
setwd("/home/bastian.egeter/Minion_data/run8/")

bas.minion.pipeline1(demuxed_fstq_dir,mastersheet,expid,run)

infastas<-list.files(pattern = "test*")
refdb = "/mnt/Disk1/Tools/BLAST+/DBs/nt_v5/nt_v5"

blast.min.bas(infastas,refdb,
              blast_exec = "/home/bastian.egeter/Tools/ncbi-blast-2.9.0+/bin/blastn",wait = F)  

#############################################################################
#FILTER BLASTS
files<-list.files(pattern = ".blast.txt")
for(i in 1:length(files)){
  blastfile = files[i]
  out<-gsub(".blast.txt",".blast.filt.txt",files[i])
  filter.blast(blastfile = blastfile,ncbiTaxDir = ncbiTaxDir,out = out)
}

#############################################################################
#BIN READS
files<-list.files(pattern = ".blast.filt.txt")
for(i in 1:length(files)){
  filtered_blastfile<-files[i]
  binfile<-gsub(".blast.filt.txt",".bins.txt",files[i])
  bin.blast2(filtered_blastfile = filtered_blastfile,ncbiTaxDir = ncbiTaxDir,
             obitaxdb = obitaxdb,out = binfile,spident = 95,gpident = 92,fpident = 80,
             abspident = 70)
}

#############################################################################
#MERGING BINS AND OTU TABLES

files<-list.files(pattern = "wlen.tab")

for(i in 1:length(files)){
   obitabfile<-files[i]
   binfile<-gsub(".tab",".bins.txt",files[i])
   out<-gsub(".tab",".taxatable.txt",files[i])
   
   obitab_bin_blast_merge_minion(obitabfile = obitabfile,binfile = binfile,mastersheetfile = mastersheetfile,
                                experiment_id = expid,out=out)
}


#Apply taxa filter
files<-list.files(pattern = "*taxatable.txt")
taxon.filter.solo(files,filterpc=0.1)

#Output taxatables by project
files<-list.files(pattern = "*taxatable.tf.txt")
splice.taxatables(files,mastersheet)

#make contributor files
files<-list.files(pattern = "*spliced.txt")
for(i in 1:length(files)){
  check.low.res.df(
    filtered.taxatab = files[i],
    filtered_blastfile = gsub("taxatable.tf.spliced.txt","blast.filt.txt",
                                    strsplit(files[i],"-")[[1]][2]),
    binfile = gsub("taxatable.tf.spliced.txt","bins.txt",strsplit(files[i],"-")[[1]][2]))
}
  