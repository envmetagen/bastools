
library(processx)
library(dplyr)

setwd(bastoolsDir)
if(use.google) googlesheets::gs_auth() 
source("blast.min.bas.R")
source("googlesheet.foos.R")
source("bin.blast.R")
source("add.taxids.fasta.BAS.R")
source("taxatab.filter.R")
source("merge_MBC_otutab_with_bin_blast.R")
source("plotting.R")

setwd(outDir)

####################################################
#step 0 master sheet creation

if("step0" %in% stepstotake){
  
  message("STEP0")
  
  experimentsheet<-google.make.exp.sheet.illumina(outDir,sheeturls,experiment_id) #in process, writes a sheet to file 
  
  message("STEP0 complete")
  
}


####################################################
#step 7 blast
if("step7" %in% stepstotake){
  
  message("STEP7")
  
  #BLASTING NT
  if(class(startingfastas)=="data.frame"){
    blast.status<-blast.min.bas(infastas = as.character(startingfastas[,1]),refdb = refdb,blast_exec = blast_exec,
                              wait = T,taxidlimit = startingfastas[,2],taxidname = as.character(startingfastas[,3]),
                              ncbiTaxDir = ncbiTaxDir)
  } else{blast.status<-blast.min.bas(infastas = startingfastas,refdb = refdb,blast_exec = blast_exec,
                                     wait = T, ncbiTaxDir = ncbiTaxDir)
  }
  
  if(class(startingfastas)=="data.frame"){
    check.blasts(infastas = as.character(startingfastas[,1]),h = blast.status)
  } else{check.blasts(infastas = startingfastas,h = blast.status)
  }
  
  message("STEP7 complete")
  
}

####################################################
#step 8 filter 
if("step8" %in% stepstotake){  
  
  message("STEP8")
  
  #FILTER BLASTS
  files<-list.files(pattern = ".blast.txt$")
                            
  
  for(i in 1:length(files)){
    message(paste("filtering blast results for",files[i]))
    blastfile = files[i]
    out<-gsub(".blast.txt",".blast.filt.txt",files[i])
    filter.blast(blastfile = blastfile,ncbiTaxDir = ncbiTaxDir,out = out, top = top)
  }
  message("STEP8 complete")
  
}

####################################################
#step 9 BIN
if("step9" %in% stepstotake){ 
  
  message("STEP9")
  
  files<-list.files(pattern = ".blast.filt.txt")
  
  for(i in 1:length(files)){
    message(paste("binning filtered blast results for",files[i]))
    filtered_blastfile<-files[i]
    binfile<-gsub(".blast.filt.txt",".bins.txt",files[i])
    bin.blast2(filtered_blastfile = filtered_blastfile,ncbiTaxDir = ncbiTaxDir,
               obitaxdb = obitaxdb,out = binfile,spident = spident,gpident = gpident,fpident = fpident,abspident = abspident)
  }
  message("STEP9 complete")
  
} 

####################################################
#step 10 merge with otutab

if("step10" %in% stepstotake){  
  
  message("STEP10")
  
  files<-list.files(pattern = ".otutab.tsv$")
  
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
  
  message("STEP11")
  
  files<-list.files(pattern = ".taxatable.txt$")
  
  taxon.filter.solo(files,filterpc)
  
  message("STEP11 complete")
  
}

####################################################
#step 12 splice taxa tables

if("step12" %in% stepstotake){  
  
  message("STEP12")
  
  files<-list.files(pattern = ".taxatable.tf.txt$")
  splice.taxatables(files,mastersheet = paste0(outDir,paste(experiment_id,collapse = "_"),"_experiment_sheet.txt"))
  
  message("STEP12 complete")
  
}

####################################################
#step 13 make contributor files

if("step13" %in% stepstotake){  
  
  message("STEP13")
  
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
  
  message("STEP14")
  
  files<-list.files(pattern = ".taxatable.tf.spliced.txt$")
  for(i in 1:length(files)){
    bas.krona.plot(files[i],KronaPath = "/home/bastian.egeter/Tools/Krona.install/bin/ktImportText")
  }
  message("STEP14 complete")
  
}


