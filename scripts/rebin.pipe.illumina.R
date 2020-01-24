message("settings:
        ")
ls.str()

library(processx)
library(dplyr)
setwd(bastoolsDir)
googlesheets4::sheets_auth(email = email) 
source("master_functions.R")
source("bin.blast.R")
setwd(outDir)

####################################################
#step9 rebin

if("step9" %in% stepstotake){  
  
  message("STEP9")
  
  for(i in 1:length(files)){
    message(paste("binning filtered blast results for",files[i]))
    filtered_blastfile<-files[i]
    binfile<-paste0(experiment_id[i],"_",gsub(".blast.filt.txt",".rebins.txt",
                                stringr::str_split(files[i],"/")[[1]][length(stringr::str_split(files[i],"/")[[1]])]))
    bin.blast2(filtered_blastfile = filtered_blastfile,ncbiTaxDir = ncbiTaxDir,
               obitaxdb = obitaxdb,out = binfile,spident = spident,gpident = gpident,fpident = fpident,
               abspident = abspident,disabledTaxaFiles = disabledTaxaFiles,disabledTaxaOut = disabledTaxaOut)
  }
  
  message("STEP9 complete")
  
}

####################################################
#step 10 merge with otutab

if("step10" %in% stepstotake){  
  
  message("STEP10")
  
  for(i in 1:length(files)){
    otutabfile<-gsub(".blast.filt.txt",".otutab.tsv",files[i])
    binfile<-paste0(experiment_id[i],"_",gsub(".blast.filt.txt",".rebins.txt",
                                stringr::str_split(files[i],"/")[[1]][length(stringr::str_split(files[i],"/")[[1]])]))
    out<-gsub(".rebins.txt",".rebins.taxatable.txt",binfile)
    MBC_otu_bin_blast_merge(MBC_otutab = otutabfile,bin_blast_results = binfile,out = out)
  }
  message("STEP10 complete")
  
}

####################################################
#step 11 apply taxon filter

if("step11" %in% stepstotake){  
  
  message("STEP11")
  
  taxatabs<-list.files(pattern = ".taxatable.txt$")
  
  taxon.filter.solo(taxatabs,filterpc)
  
  message("STEP11 complete")
  
}

####################################################
#step 12 cut out project samples 

if("step12" %in% stepstotake){  
  
  message("STEP12")
  
  message("Note:sample names must contain project name with dash")
  
  taxatabs<-list.files(pattern = ".taxatable.tf.txt$")
  for(i in 1:length(taxatabs)){
    taxatable<-data.table::fread(taxatabs[i],sep = "\t",header = T,data.table = F)
    taxatable<-taxatable[,c(1,grep(project_name,colnames(taxatable)))]
    taxatable<-taxatable[rowSums(taxatable[,2:length(colnames(taxatable))])!=0,]
    write.table(taxatable,
                paste0(project_name,gsub("taxatable.tf.txt","taxatable.tf.spliced.txt",taxatabs[i])),
                row.names = F,quote = F,sep = "\t")
  }
  
  message("STEP12 complete")
  
}

####################################################
#step 13 make contributor files

if("step13" %in% stepstotake){  
  
  message("STEP13")
  
  #make contributor files
  for(i in 1:length(files)){
    message(paste("making contributor file for",files[i]))
    
    #first get taxatab file names
    filtered.taxatab = paste0(project_name,experiment_id[i],"_" ,gsub(".blast.filt.txt",".rebins.taxatable.tf.spliced.txt",
       stringr::str_split(files[i],"/")[[1]][length(stringr::str_split(files[i],"/")[[1]])]))
    
    check.low.res.df(
      filtered.taxatab = filtered.taxatab,filtered_blastfile = files[i],
      binfile<-paste0(experiment_id[i],"_",gsub(".blast.filt.txt",".rebins.txt",
                    stringr::str_split(files[i],"/")[[1]][length(stringr::str_split(files[i],"/")[[1]])]))
      ,disabledTaxaFile = disabledTaxaOut,spident = spident,gpident = gpident,fpident = fpident,abspident = abspident)
  }
  message("STEP13 complete")
  
}

####################################################
#step 14 make krona plots

if("step14" %in% stepstotake){  
  
  message("STEP14")
  
  files<-list.files(pattern = ".taxatable.tf.spliced.txt$")
  for(i in 1:length(files)){
    bas.krona.plot(files[i],KronaPath)
  }
  message("STEP14 complete")
  
}
