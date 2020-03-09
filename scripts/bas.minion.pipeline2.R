message("settings:
        ")
print(ls.str())

library(processx)
library(dplyr)
source(paste0(bastoolsDir,"master_functions.R"))
source(paste0(bastoolsDir,"bin.blast.R"))
setwd(filesDir)

####################################################

#step 0 master sheet creation

if("step0" %in% stepstotake){
  
  message("STEP0")
 
  mastersheet_df<-data.table::fread(mastersheet,data.table = F,sep = "\t")
  
  ##subset master sheet for this study 
  message("Subsetting datasheet")
  ms_ss<-subset_mastersheet(mastersheet_df, subsetlist)
  #check again to see the subset made sense
  master_xtabs(master_sheet = ms_ss,columns=c("sample_type",names(subsetlist)))
  
  message("STEP0 complete")

}

####################################################
#step 1 Cutadapt
if("step1" %in% stepstotake){
  message("STEP1 - Cutadapt")
  
  files<-paste0(ms_ss$barcode_name,".fastq.gz")
  
  cutadaptlist<-list()
  
  for(i in 1:length(files)){
    file<-files[i]
    message("Running cutadapt on ",file)
    if(!file %in% list.files(path = filesDir)) stop("File not found")
    barcode<-ms_ss$barcode_name[i]
    Pf<-ms_ss[match(ms_ss$barcode_name[i],ms_ss$barcode_name),"Primer_F"]
    Pr<-ms_ss[match(ms_ss$barcode_name[i],ms_ss$barcode_name),"Primer_R"]
    
    message("Pf=",Pf)
    message("Pr=",Pr)
    
    cutadaptlist[[i]]<-system2(command = "cutadapt",args = c("-g", paste0("'Or1=",Pf,";error_rate=0.2...",insect::rc(Pr),";max_error_rate=0.2'"), "-g",
                                          paste0("'Or2=",Pr,";error_rate=0.2...",insect::rc(Pf),";max_error_rate=0.2'")
                                          ,"-o", gsub(".fastq.gz",".trimmed.fastq",file),file),wait = T,stderr = T,stdout = T)
    
    print(grep("Total reads processed",cutadaptlist[[i]],value = T))
    print(grep("Reads with adapters",cutadaptlist[[i]],value = T))
  }
  
  message("STEP1 complete")
}

####################################################
#step 2 print length data
if("step2" %in% stepstotake){
  
  message("STEP2 - print length data")
  
  files<-paste0(ms_ss$barcode_name,".trimmed.fastq")
  
  #print histograms first
  for(i in 1:length(files)){
    message(files[i])
    system2("seqkit", args=c("watch","--dump",files[i]),wait = T)
  }
  
  message("STEP2 complete")
}

####################################################
#step 3 size select
if("step3" %in% stepstotake){
  
  message("STEP3 - size select")
  
  files<-paste0(ms_ss$barcode_name,".trimmed.fastq")
  
  #convert to fasta first
  for(i in 1:length(files)){
    file=files[i]
    message("Converting fastq to fasta for ",file)
    system2("seqkit", args=c("fq2fa",file),wait = F,stdout = gsub(".fastq",".fasta",file))
  }
  
  files<-paste0(ms_ss$barcode_name,".trimmed.fasta")
  
  #size select
  for(i in 1:length(files)){
    file=files[i]
    message("Running size selection on ",file)
    
    minlength<-ms_ss[match(ms_ss$barcode_name[i],ms_ss$barcode_name),"min_length"]
    maxlength<-ms_ss[match(ms_ss$barcode_name[i],ms_ss$barcode_name),"max_length"]
    
    message("minlength=",minlength)
    message("maxlength=",maxlength)
    
    xx<-system2("seqkit", args=c("seq","-M",maxlength,"-m",minlength,file),wait = T, stdout = gsub(".trimmed.fasta",".lenFilt.trimmed.fasta",file)
                ,stderr = F)
    
    #count seqs original
    a<-system2("seqkit", args=c("stats","-T",file),wait = T,stderr = T,stdout = T)
    a<-read.table(text = a,header = T)$num_seqs
    #count seqs after
    b<-system2("seqkit", args=c("stats","-T",gsub(".trimmed.fasta",".lenFilt.trimmed.fasta",file)),wait = T,stderr = T,stdout = T)
    b<-read.table(text = b,header = T)$num_seqs
    message(b," reads removed, from ",a, " (",round(b/a*100,digits = 2)," %)")
  }
  message("STEP3 complete")
}
  
####################################################
#step 4 
if("step4" %in% stepstotake){
  
  message("STEP4 - blast")
  
  files<-paste0(ms_ss$barcode_name,".lenFilt.trimmed.fasta")  
  
  for(i in 1:length(infastas)){
    if(!is.null(taxidlimit)){

    blast.status<-blast.min.bas(infastas = files[i],refdb = refdb,blast_exec = blast_exec,
                                wait = T,taxidlimit = taxidlimit,taxidname = taxidname, ncbiTaxDir = ncbiTaxDir,)
  } else{blast.status<-blast.min.bas(infastas = startingfastas,refdb = refdb,blast_exec = blast_exec, wait = T, ncbiTaxDir = ncbiTaxDir)
  }
  
  if(class(startingfastas)=="data.frame"){
    check.blasts(infastas = as.character(startingfastas[,1]),h = blast.status)
  } else{check.blasts(infastas = startingfastas,h = blast.status)
  }
  
  
  message("STEP4 complete")
    
    
    
            
  
  
  
    
  


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
#step 4 dereplicate
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
  
  if("step4" %in% stepstotake) {
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
  
  if(class(startingfastas)=="data.frame"){
  check.blasts(infastas = as.character(startingfastas[,1]),h = blast.status)
  } else {check.blasts(infastas = as.character(startingfastas),h = blast.status)}
  
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
  
  taxon.filter.solo(files = files,filterpc = filterpc)
  
  message("STEP11 complete")
  
}

####################################################
#step 12 splice taxa tables

if("step12" %in% stepstotake){  
  
  message("STEP12")
  
  files<-grep(experiment_id,list.files(path = outDir,pattern = ".taxatable.tf.txt$"),value = T)
  splice.taxatables(files = files,mastersheet = paste0(outDir,experiment_id,"_experiment_sheet.txt"))
  
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
      binfile = list.files(pattern = gsub("blast.filt.txt","bins.txt",filtered_blastfile)),
      spident = spident, gpident = gpident,fpident = fpident,abspident = abspident)
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
    
    
