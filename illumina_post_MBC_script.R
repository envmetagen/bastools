#put all final fastas (usually one per primer) in a folder called "final_fastas", with no other fastas.
#put all final otutabs (usually one per primer) in a folder called "final_otutabs", with no other tabs.
#the parent directory should be empty also, as the outputs will be put here

library(processx)
library(dplyr)
library(mgsub)

setwd("/home/bastian.egeter/git_bastools/bastools/")
file.sources<-c("add.taxids.fasta.BAS.R","bin.blast.R","merge_MBC_otutab_with_bin_blast.R","blast.min.bas.R",
                "taxatab.filter.R","do.count.R","googlesheet.foos.R")
sapply(file.sources,source)

setwd(rootdir)

setwd("./final_fastas")

#BLASTING NT
files<-list.files(pattern = "*.fasta")
blast.status<-blast.min.bas(infastas = files,refdb = refdb,blast_exec) 

check.blasts(infastas = files,h = blast.status)

files<-list.files(pattern = ".blast.txt")
setwd("../")
dir.create("blasts")
file.copy(list.files(full.names = T,path = "./final_fastas", pattern = ".blast.txt"),
          "./blasts")
file.remove(list.files(full.names = T,path = "./final_fastas", pattern = ".blast.txt"))
#############################################################################
#FILTER BLASTS
setwd("./blasts/")
files<-list.files(pattern = ".blast.txt")
for(i in 1:length(files)){
  message(paste("filtering blast results for",files[i]))
  blastfile = files[i]
  out<-gsub(".blast.txt",".blast.filt.txt",files[i])
  filter.blast(blastfile = blastfile,ncbiTaxDir = ncbiTaxDir,out = out)
}
#############################################################################
#BIN READS
files<-list.files(pattern = ".blast.filt.txt")
for(i in 1:length(files)){
  message(paste("binning filtered blast results for",files[i]))
  filtered_blastfile<-files[i]
  binfile<-gsub(".blast.filt.txt",".bins.txt",files[i])
  bin.blast2(filtered_blastfile = filtered_blastfile,ncbiTaxDir = ncbiTaxDir,
            obitaxdb = obitaxdb,out = binfile)
}
files<-list.files(pattern = ".bins.txt")
setwd("../")
dir.create("bins")
file.copy(list.files(full.names = T,path = "./blasts", pattern = ".bins.txt"),
          "./bins")
file.remove(list.files(full.names = T,path = "./blasts", pattern = ".bins.txt"))
#############################################################################
#MERGING BINS AND OTU TABLES
setwd("./final_otutabs/")
files<-list.files(pattern = ".tsv")

for(i in 1:length(files)){
  otutabfile<-files[i]
  binfile<-gsub(".tsv",".bins.txt",files[i])
  binfile<-list.files(pattern = binfile,path = "../bins",full.names = T)
  mergefile<-gsub(".tsv",".taxatable.txt",files[i])
  message(paste("mergeing bins and otutabs for",binfile, "and" , otutabfile))
  MBC_otu_bin_blast_merge(MBC_otutab = otutabfile,bin_blast_results = binfile,out = mergefile)
}
setwd("../")
dir.create("taxatables")
file.copy(list.files(full.names = T,path = "./final_otutabs", pattern = ".taxatable.txt"),
          "./taxatables")
file.remove(list.files(full.names = T,path = "./final_otutabs", pattern = ".taxatable.txt"))
#############################################################################
# Note: 
# If all hits for a particular OTU are removed due to filters, the results will be NA;NA;NA;NA;NA;NA;NA 
# If the lca for a particular OTU is above kingdom, the results will be unknown;unknown;unknown;unknown;NA;NA;NA  
# If there were no blast hits, the results will be No_hits;No_hits;No_hits;No_hits;No_hits;No_hits;No_hits
#############################################################################

#Apply taxa filter
setwd("taxatables/")
files<-list.files(pattern = "*taxatable.txt")
message(paste("applying taxon filters for all taxatables"))
taxon.filter.solo(files,filterpc=0.1)

#Output taxatables by project
files<-list.files(pattern = "*taxatable.tf.txt")
message(paste("splicing taxatables by project"))
splice.taxatables(files,mastersheet)

#make contributor files
files<-list.files(pattern = "*spliced.txt")
for(i in 1:length(files)){
  message(paste("making contributor file for",files[i]))
check.low.res.df(
  filtered.taxatab = files[i],
  filtered_blastfile = list.files(path = "../blasts",full.names = T,
   pattern = gsub("taxatable.tf.spliced.txt","blast.filt.txt",strsplit(files[i],"-")[[1]][2])),
  binfile = list.files(path = "../bins",full.names = T,
   pattern = gsub("taxatable.tf.spliced.txt","bins.txt",strsplit(files[i],"-")[[1]][2])))
}

#inspect contributor files by eye and make list of taxids to exclude...then build this into bin.blast.R

