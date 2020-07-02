message("settings:
        ")
print(ls.str())

library(processx)
library(dplyr)
source(paste0(bastoolsDir,"master_functions.R"))
setwd(outDir)

####################################################
#step 7 blast
if("step7" %in% stepstotake){
  
  message("STEP7-BLAST")
  
  message("Note to self, consider using different blastn options")
  
  #BLASTING NT
  if(class(startingfastas)=="data.frame") {
    stop ("not written this way for pipe3")
  } else{
    for(i in 1:length(startingfastas)){
    blast.status<-blast.min.bas2(infasta = startingfastas[i],refdb = refdb,blast_exec = blast_exec, wait = T,
                               taxidlimit = taxidlimit, ncbiTaxDir = ncbiTaxDir,opts = opts,overWrite = T)
    }
  }
  
  message("STEP7 complete")
  
}

####################################################
#step 8 filter 
if("step8" %in% stepstotake){  
  
  message("STEP8-filter blast")
  
  #startingfastas<-c("12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.fasta"
   #                 ,"16S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.fasta")
  #outDir="/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/2018_07/"
  
  
  #FILTER BLASTS
  files<-gsub(".fasta$",".blast.txt",startingfastas)
  
  for(i in 1:length(files)){
    message(paste("filtering blast results for",files[i]))
    blastfile = files[i]
    out<-gsub(".blast.txt",".blast.filt.txt",files[i])
    filter.blast3(blastfile = blastfile,ncbiTaxDir = ncbiTaxDir,out = out)
  }
  message("STEP8 complete")
  
}

####################################################
#step 9 BIN
if("step9" %in% stepstotake){ 
  
  message("STEP9 - bin")
  
  files<-gsub(".fasta$",".blast.filt.txt",startingfastas)
  
  for(i in 1:length(files)){
    message(paste("binning filtered blast results for",files[i]))
    filtered_blastfile<-files[i]
    binfile<-gsub(".blast.filt.txt",".bins.txt",files[i])
    bin.blast2(filtered_blastfile = filtered_blastfile,ncbiTaxDir = ncbiTaxDir,
               out = binfile,spident = spident,gpident = gpident,
               fpident = fpident,abspident = abspident)
  }
  message("STEP9 complete")
  
} 

####################################################
#step 10 merge with otutab

if("step10" %in% stepstotake){  
  
  message("STEP10 - merge with otutab")
  
  files<-gsub(".fasta$",".otutab.tsv",startingfastas)
  
  for(i in 1:length(files)){
    otutabfile<-files[i]
    binfile<-gsub(".otutab.tsv",".bins.txt",files[i])
    out<-gsub(".otutab.tsv",".taxatable.txt",files[i])
    
    MBC_otu_bin_blast_merge(MBC_otutab = otutabfile,bin_blast_results = binfile,out = out)
  }
  message("STEP10 complete")
  
}

####################################################
#step 11 apply taxon filter

if("step11" %in% stepstotake){  
  
  message("STEP11 - taxon filter")
  message("CHANGE SCRIPT TO ONLY TAKE RELEVANT FILES")
  
  files<-gsub(".fasta$",".taxatable.txt",startingfastas)
  
  taxon.filter.solo(files,filterpc)
  
  message("STEP11 complete")
  
}

####################################################
#step 12 splice taxa tables

if("step12" %in% stepstotake){  
  
  message("STEP12 - splice tables")
  message("CHANGE SCRIPT TO ONLY TAKE RELEVANT FILES")
  
  files<-gsub(".fasta$",".taxatable.tf.txt",startingfastas)
  
  splice.taxatables(files)
  
  message("STEP12 complete")
  
}

####################################################
#step 13 make contributor files

if("step13" %in% stepstotake){  
  
  message("STEP13 - make contributor files")
  
  message("CHANGE SCRIPT TO ONLY TAKE RELEVANT FILES")
  
  #files<-gsub(".fasta$",".taxatable.tf.spliced.txt",startingfastas)
  
  files<-list.files(pattern = ".taxatable.tf.spliced.txt$")
  
  #make contributor files
  for(i in 1:length(files)){
    message(paste("making contributor file for",files[i]))
    
    #first get blast file names (accounting for extra dashes)
    if(length(strsplit(files[i],"-")[[1]])==2) { 
      filtered_blastfile = list.files(pattern = gsub("taxatable.tf.spliced.txt","blast.filt.txt",
                                                     strsplit(files[i],"-")[[1]][2]))}
    
    if(length(strsplit(files[i],"-")[[1]])>2) { 
      filtered_blastfile = list.files(pattern = gsub("taxatable.tf.spliced.txt","blast.filt.txt",
                                                     paste(strsplit(files[i],"-")[[1]][2:length(strsplit(files[i],"-")[[1]])],
                                                           collapse = "-")))
    }
    
    check.low.res.df(
      filtered.taxatab = files[i],filtered_blastfile = filtered_blastfile,
      binfile = list.files(pattern = gsub("blast.filt.txt","bins.txt",filtered_blastfile))
      ,disabledTaxaFile = NULL,spident = spident,gpident = gpident,fpident = fpident,abspident = abspident)
  }
  message("STEP13 complete")
  
}

####################################################
#step 14 make krona plots

if("step14" %in% stepstotake){  
  
  message("STEP14 - krona plots")
  message("CHANGE SCRIPT TO ONLY TAKE RELEVANT FILES")
  
  files<-list.files(pattern = ".taxatable.tf.spliced.txt$")
  for(i in 1:length(files)){
    bas.krona.plot(files[i],KronaPath)
  }
  message("STEP14 complete")
  
}


