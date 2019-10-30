
library(processx)
library(dplyr)

setwd(bastoolsDir)
googlesheets::gs_auth() 
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
  
  experimentsheet<-google.make.experiment.sheet(outDir,sheeturls,experiment_id) #in process, writes a sheet to file 
  
  message("STEP0 complete")
  
}


####################################################
#step 7 blast
if("step7" %in% stepstotake){
  
  message("STEP7")
  
  #BLASTING NT
  blast.status<-blast.min.bas(infastas = as.character(startingfastas[,1]),refdb = refdb,blast_exec = blast_exec,
                              wait = T,taxidlimit = startingfastas[,2],taxidname = as.character(startingfastas[,3]),
                              ncbiTaxDir = ncbiTaxDir)
  
  check.blasts(infastas = files,h = blast.status)
  
  message("STEP7 complete")
  
}

####################################################
#step 8 filter 
if("step8" %in% stepstotake){  
  
  message("STEP8")
  
  #FILTER BLASTS
  files<-paste0(outDir,list.files(path = outDir,pattern = ".blast.txt"))
                            
  
  for(i in 1:length(files)){
    message(paste("filtering blast results for",files[i]))
    blastfile = files[i]
    out<-gsub(".blast.txt",".blast.filt.txt",files[i])
    filter.blast(blastfile = blastfile,ncbiTaxDir = ncbiTaxDir,out = out)
  }
  message("STEP8 complete")
  
}

####################################################
#step 9 BIN
if("step9" %in% stepstotake){ 
  
  message("STEP9")
  
  files<-paste0(outDir,grep(experiment_id,list.files(path = outDir,pattern = ".blast.filt.txt"),value = T))
  
  for(i in 1:length(files)){
    message(paste("binning filtered blast results for",files[i]))
    filtered_blastfile<-files[i]
    binfile<-gsub(".blast.filt.txt",".bins.txt",files[i])
    bin.blast2(filtered_blastfile = filtered_blastfile,ncbiTaxDir = ncbiTaxDir,
               obitaxdb = obitaxdb,out = binfile)
  }
  message("STEP9 complete")
  
} 

####################################################
#step 10 merge with otutab

if("step10" %in% stepstotake){  
  
  message("STEP10")
  
  files<-paste0(outDir,grep(experiment_id,list.files(path = outDir,pattern = ".tab$"),value = T))
  binfiles<-paste0(outDir,grep(experiment_id,list.files(path = outDir,pattern = ".bins.txt$"),value = T))
  
  for(i in 1:length(files)){
    otutabfile<-files[i]
    binfile<-gsub(".tab",".bins.txt",files[i])
    out<-gsub(".tab",".taxatable.txt",files[i])
    
    otutab_bin_blast_merge_minion(otutabfile = otutabfile,binfile = binfile,experimentsheetfile = 
                                    paste0(outDir,experiment_id,"_experiment_sheet.txt"),
                                  hascount = hascount,  experiment_id = experiment_id,out=out)
  }
  message("STEP10 complete")
  
}

####################################################
#step 11 apply taxon filter

if("step11" %in% stepstotake){  
  
  message("STEP11")
  
  files<-paste0(outDir,grep(experiment_id,list.files(path = outDir,pattern = ".taxatable.txt$"),value = T))
  
  taxon.filter.solo(files,filterpc)
  
  message("STEP11 complete")
  
}

####################################################
#step 12 splice taxa tables

if("step12" %in% stepstotake){  
  
  message("STEP12")
  
  files<-grep(experiment_id,list.files(path = outDir,pattern = ".taxatable.tf.txt$"),value = T)
  splice.taxatables(files,mastersheet = paste0(outDir,experiment_id,"_experiment_sheet.txt"))
  
  message("STEP12 complete")
  
}

####################################################
#step 13 make contributor files

if("step13" %in% stepstotake){  
  
  message("STEP13")
  
  files<-grep(experiment_id,list.files(path = outDir,pattern = ".taxatable.tf.spliced.txt$"),value = T)
  
  #make contributor files
  for(i in 1:length(files)){
    message(paste("making contributor file for",files[i]))
    
    #first get blast file names (accounting for extra dashes)
    if(length(strsplit(files[i],"-")[[1]])==2) { 
      filtered_blastfile = list.files(pattern = gsub("taxatable.tf.spliced.txt","blast.filt.txt",
                                                     strsplit(files[i],"-")[[1]][2]))}
    
    if(length(strsplit(files[i],"-")[[1]])>2) { 
      filtered_blastfile = list.files(pattern = gsub("taxatable.tf.spliced.txt","blast.filt.txt",
                                                     paste(strsplit(files[i],"-")[[1]][2:length(strsplit(files[i],"-")[[1]])],collapse = "-")))
    }
    
    #then get bin file names (accounting for extra dashes)
    
    check.low.res.df(
      filtered.taxatab = files[i],filtered_blastfile,
      binfile = list.files(pattern = gsub("blast.filt.txt","bins.txt",filtered_blastfile)))
  }
  message("STEP13 complete")
  
}

####################################################
#step 14 make krona plots

if("step14" %in% stepstotake){  
  
  message("STEP14")
  
  files<-grep(experiment_id,list.files(path = outDir,pattern = ".taxatable.tf.spliced.txt$"),value = T)
  for(i in 1:length(files)){
    bas.krona.plot(files[i],KronaPath = "/home/bastian.egeter/Tools/Krona.install/bin/ktImportText")
  }
  message("STEP14 complete")
  
}


