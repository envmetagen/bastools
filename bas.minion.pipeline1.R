
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

#get barcodes used  
barcodes.used<-unique(experimentsheet$barcode_id)
barcodes.used <- barcodes.used[!is.na(barcodes.used)]
barcodes.used<-gsub("BC","barcode",barcodes.used)

#size select, for each fragment, I checked and seqs appear to have primers plus one base (at each end)
experimentsheet$primer_combo<-paste0(experimentsheet$Primer_F,experimentsheet$Primer_R)
primer_combo<-unique(experimentsheet$primer_combo)

#barcodes in each primer combo
primer_combo.bcs<-list()
for(i in 1:length(primer_combo)){
  primer_combo.bcs[[i]]<-unique(experimentsheet[experimentsheet$primer_combo==primer_combo[i],"barcode_id"])
  primer_combo.bcs[[i]]<-gsub("BC","barcode",primer_combo.bcs[[i]])
  names(primer_combo.bcs[[i]])<-gsub(" ","",
                                     experimentsheet[experimentsheet$primer_combo==primer_combo[i],"Primer_set"][1])
}

#size calc
message("size calc includes primer lengths at the moment")
minlength<-list()
for(i in 1:length(primer_combo)){
  minlength[[i]]<-sum(experimentsheet[experimentsheet$primer_combo==primer_combo[i],
                                      "Min_length"][1],nchar(primer_combo[i]))
  names(minlength[[i]])<-names(primer_combo.bcs[[i]][1])
}

maxlength<-list()
for(i in 1:length(primer_combo)){
  maxlength[[i]]<-sum(experimentsheet[experimentsheet$primer_combo==primer_combo[i],
                                      "Max_length"][1],nchar(primer_combo[i]))
  names(maxlength[[i]])<-names(primer_combo.bcs[[i]][1])
}

message("STEP0 complete")

}

####################################################
#step 1 extract gz files and convert to fasta
if("step1" %in% stepstotake){
  message("STEP1")
  
#extract gz to new files
if(in_folders) {
  for(i in 1:length(barcodes.used)){
    folder<-paste0(demuxed_fstq_dir,barcodes.used[i],"/")
    system2(command = "zcat",args = c(paste0(folder,list.files(folder))),
            stdout = paste0(outDir,barcodes.used[i],".fastq"),wait = T)
  }
} else stop ("Not written yet")

#these files should be the ones for DMP
#file.rename(list.files(pattern="water_*.img"), paste0("water_", 1:700))

  #convert to fasta
  files<-list.files(path = outDir,pattern = "*.fastq")
  for(i in 1:length(files)){
  system2(command = "sed",args = c("-n", "'1~4s/^@/>/p;2~4p'",paste0(outDir,files[i])),
          wait=T,stdout= gsub(".fastq",".fasta",paste0(outDir,files[i])))
  }
  
  #reformat for obitools
  files<-list.files(path = outDir,pattern = "*.fasta")
  for(i in 1:length(files)){
    system2(command = "sed",
            args = c("'s/ /; /g;s/; / /1;s/   / /g;s/barcode=barcode\\([0-9][0-9]\\)/barcode=barcode\\1;/g'",
                     paste0(outDir,files[i])),
            wait=T,stdout= gsub(".fasta",".obi.fasta",paste0(outDir,files[i])))
  }
  
  #add ss_sample_id titles to headers
  for(i in 1:length(barcodes.used)){
    file<-paste0(outDir,barcodes.used[i],".obi.fasta")
    ss_sample_id<-experimentsheet[experimentsheet$barcode_id==gsub("barcode","BC",barcodes.used[i]),"ss_sample_id"]
    system2(command = "sed",
    args = c(paste0("'s/\\(barcode=barcode[0-9][0-9];\\)/\\1 ss_sample_id=",ss_sample_id,";/g'"),
                     paste0(file)),
            wait=T,stdout= gsub(".obi.fasta",".ss.obi.fasta",paste0(file)))
  }
  
  #remove intermediates
  unlink(paste0(outDir,barcodes.used,".fastq"))
  unlink(paste0(outDir,barcodes.used,".fasta"))
  unlink(paste0(outDir,list.files(path = outDir,pattern = "barcode[0-9][0-9].obi.fasta")))
  
  message("STEP1 complete")
}
####################################################
#step 2 cat by frag
if("step2" %in% stepstotake){  
  
  message("STEP2")
  
    #cat each frag
    files<-list.files(path = outDir,pattern = "*.ss.obi.fasta")
    for(i in 1:length(primer_combo.bcs)){
      primercombofiles<-files[files %in% paste0(primer_combo.bcs[[i]],".ss.obi.fasta")]
      system2(command = "cat",args = c(paste0(outDir,primercombofiles)),
              wait = T,stdout = paste0(outDir,experiment_id,"_",names(primer_combo.bcs[[i]])[1],
                                       ".fasta")) 
    }
    
    message("STEP2 complete")
}    
####################################################
#step 3 size select
if("step3" %in% stepstotake){
  
  message("STEP3")
  
  #add length to fasta
    
    #add seq lengths
    files<-grep(experiment_id,list.files(path = outDir,pattern = ".fasta"),value = T)
    for(i in 1:length(files)){
    system2(command = "obiannotate",args = c("--length", paste0(outDir,files[i])),wait = T,
            stdout= gsub(".fasta",".wlen.obi.fasta",paste0(outDir,files[i])))
    }
    
  #obigrep lengths
  files<-grep(experiment_id,list.files(path = outDir,pattern = ".wlen.obi.fasta"),value = T)
  for(i in 1:length(minlength)){
    file<-paste0(outDir,grep(names(minlength[[i]]),files,value = T))
    system2(command = "obigrep",args = c("-l",minlength[[i]],"-L",maxlength[[i]], file),wait = T,
              stdout =  gsub(".wlen.obi.fasta",".filtlen.wlen.obi.fasta",file ))
  }
  
  message("STEP3 complete")
  
} 


####################################################
#step 4 demultiplex
if("step4" %in% stepstotake){
  
  message("STEP4")

#do obiuniq
files<-list.files(path = outDir,pattern="*.filtlen.wlen.obi.fasta")
for(i in 1:length(files)){
  system2(command = "obiuniq",args = c("-c", "ss_sample_id",paste0(outDir,files[i])),wait = T,
          stdout =  gsub(".filtlen.wlen.obi.fasta",".uniq.filtlen.wlen.obi.fasta",paste0(outDir,files[i])))
}
message("STEP4 complete")

}
####################################################
#step 5 make otu tabs
if("step5" %in% stepstotake){
  
  message("STEP5")
  
  if("step5" %in% stepstotake) {
    files<-list.files(path = outDir,pattern="*.uniq.filtlen.wlen.obi.fasta")} else {
      files<-list.files(path = outDir,pattern="*.filtlen.wlen.obi.fasta")
    }
  
  #make obitabs
  
  for(i in 1:length(files)){
    system2(command = "obitab",args = c(paste0(outDir,files[i])),wait = T,
              stdout =  gsub(".fasta",".tab",paste0(outDir,files[i])))
  }
  message("STEP5 complete")
  
}
  
####################################################
# #step 6 create and tidy folders
# if("step6" %in% stepstotake){
# 
#   dir.create(paste0(outDir,"post.minion.pipe"))
#   dir.create(paste0(outDir,"post.minion.pipe/final_fastas"))
#   dir.create(paste0(outDir,"post.minion.pipe/final_otutabs"))
#   dir.create(paste0(outDir,"post.minion.pipe/blasts"))
#   
#   #move otutabs 
#   files<-list.files(path = outDir,pattern="*uniq.filtlen.wlen.obi.tab")
#   for(i in 1:length(files)){
#     file.copy(paste0(outDir,files[i]),paste0(outDir,"post.minion.pipe/final_otutabs"))
#     unlink(paste0(outDir,files[i]))
#   }
#   
#   #move fastas 
#   files<-list.files(path = outDir,pattern="*uniq.filtlen.wlen.obi.fasta")
#   for(i in 1:length(files)){
#     file.copy(paste0(outDir,files[i]),paste0(outDir,"post.minion.pipe/final_fastas"))
#     unlink(paste0(outDir,files[i]))
#   }
# } 
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
  
  check.blasts(infastas = as.character(startingfastas[,1]),h = blast.status)
  
  message("STEP7 complete")
}
 
####################################################
#step 8 filter 
if("step8" %in% stepstotake){  
  
  message("STEP8")
  
  #FILTER BLASTS
  files<-paste0(outDir,grep(experiment_id,list.files(path = outDir,pattern = ".filtlen.wlen.obi.blast.txt"),
                            value = T))
  
  for(i in 1:length(files)){
    message(paste("filtering blast results for",files[i]))
    blastfile = files[i]
    out<-gsub(".blast.txt",".blast.filt.txt",files[i])
    filter.blast(blastfile = blastfile,ncbiTaxDir = ncbiTaxDir,out = out,top = top)
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
               obitaxdb = obitaxdb,out = binfile,spident = spident,gpident = gpident,fpident = fpident,abspident = abspident)
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
    
    
