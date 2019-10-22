
library(processx)
library(dplyr)

####################################################
#step 0 master sheet creation

if("step0" %in% stepstotake){
  
#read master sheets
master<-list()
headers<-c("barcode_id","Primer_set","Primer_F","Primer_R","Min_length","Max_length","ss_sample_id","experiment_id")
for(i in 1:length(sheeturls)){
master[[i]]<-google.read.master.url(sheeturls[i])
if(length(headers)!=sum(headers %in% colnames(master[[i]]))){
  stop (c("one of the following headers missing: ", paste(headers)))}
master[[i]]<-master[[i]][,headers]
}

#make a processing sheet
experimentsheet<-as.data.frame(data.table::rbindlist(master))
experimentsheet<-experimentsheet[experimentsheet$experiment_id==experiment_id,]

#get barcodes used  
barcodes.used<-unique(experimentsheet$barcode_id)
barcodes.used <- barcodes.used[!is.na(barcodes.used)]
barcodes.used<-gsub("BC","barcode",barcodes.used)

#size select, for each fragment, I checked an seqs appear to have primers plus one base (at each end)
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

}

####################################################
#step 1 extract gz files and convert to fasta
if("step1" %in% stepstotake){
  
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
}
####################################################
#step 2 cat by frag
if("step2" %in% stepstotake){  
    #cat each frag
    files<-list.files(path = outDir,pattern = "*.ss.obi.fasta")
    for(i in 1:length(primer_combo.bcs)){
      primercombofiles<-files[files %in% paste0(primer_combo.bcs[[i]],".ss.obi.fasta")]
      system2(command = "cat",args = c(paste0(outDir,primercombofiles)),
              wait = T,stdout = paste0(outDir,experiment_id,"_",names(primer_combo.bcs[[i]])[1],
                                       ".fasta")) 
    }
    
}    
####################################################
#step 3 size select
if("step3" %in% stepstotake){
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
} 


####################################################
#step 4 demultiplex
if("step4" %in% stepstotake){

#do obiuniq
files<-list.files(path = outDir,pattern="*.filtlen.wlen.obi.fasta")
for(i in 1:length(files)){
  system2(command = "obiuniq",args = c("-c", "ss_sample_id",paste0(outDir,files[i])),wait = T,
          stdout =  gsub(".filtlen.wlen.obi.fasta",".uniq.filtlen.wlen.obi.fasta",paste0(outDir,files[i])))
}

}
####################################################
#step 5 make otu tabs
if("step5" %in% stepstotake){
  if("step4" %in% stepstotake) {
    files<-list.files(path = outDir,pattern="*.uniq.filtlen.wlen.obi.fasta")} else {
      files<-list.files(path = outDir,pattern="*.filtlen.wlen.obi.fasta")
    }
  
  #make obitabs
  
  for(i in 1:length(files)){
    system2(command = "obitab",args = c(paste0(outDir,files[i])),wait = T,
              stdout =  gsub(".fasta",".tab",paste0(outDir,files[i])))
  }

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

  #BLASTING NT
  files<-paste0(outDir,grep(experiment_id,list.files(path = outDir,pattern = ".filtlen.wlen.obi.fasta"),value = T))
  blast.status<-blast.min.bas(infastas = files,refdb = refdb,blast_exec = blast_exec) 
  check.blasts(infastas = files,h = blast.status)
  
}
 
####################################################
#step 8 filter and bin 
if("step8" %in% stepstotake){  

  #############################################################################
  #FILTER BLASTS
  files<-paste0(outDir,grep(experiment_id,list.files(path = outDir,pattern = ".filtlen.wlen.obi.blast.txt"),
                            value = T))
  
  for(i in 1:length(files)){
    message(paste("filtering blast results for",files[i]))
    blastfile = files[i]
    out<-gsub(".blast.txt",".blast.filt.txt",files[i])
    filter.blast(blastfile = blastfile,ncbiTaxDir = ncbiTaxDir,out = out)
  }
  #############################################################################
  #BIN READS
  files<-grep(experiment_id,list.files(path = outDir,pattern = ".blast.filt.txt"),value = T)
  
  for(i in 1:length(files)){
    message(paste("binning filtered blast results for",files[i]))
    filtered_blastfile<-files[i]
    binfile<-gsub(".blast.filt.txt",".bins.txt",files[i])
    bin.blast2(filtered_blastfile = filtered_blastfile,ncbiTaxDir = ncbiTaxDir,
               obitaxdb = obitaxdb,out = binfile)
  }
} 

####################################################
#step 9 merge with obitab

if("step8" %in% stepstotake){  
  message("not done yet")
  files<-grep(experiment_id,list.files(path = outDir,pattern = ".tab"),value = T)
  
  binfiles<-grep(experiment_id,list.files(path = outDir,pattern = ".tab"),value = T)
  
  for(i in 1:length(files)){
    obitabfile<-files[i]
    binfile<-gsub(".tab",".bins.txt",files[i])
    out<-gsub(".tab",".taxatable.txt",files[i])
    
    obitab_bin_blast_merge_minion(obitabfile = obitabfile,binfile = binfile,mastersheetfile = mastersheetfile,
                                  experiment_id = expid,out=out)
  }
  
  obitab_bin_blast_merge_minion(paste0(outDir,"post.minion.pipe/final_otutabs"),
                                binfile,mastersheetfile=NA,experiment_id,used.obiuniq=F,out)
}


    
    
