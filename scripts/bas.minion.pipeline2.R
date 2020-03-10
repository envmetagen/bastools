message("settings:
        ")
print(ls.str())

library(processx)
library(dplyr)
library(ggplot2)
source(paste0(bastoolsDir,"master_functions.R"))
source(paste0(bastoolsDir,"bin.blast.R"))
setwd(filesDir)

####################################################

#step 0 subset master sheet 

  message("STEP0 - subset mastersheet")
 
  mastersheet_df<-data.table::fread(mastersheet,data.table = F,sep = "\t")
  
  ##subset master sheet for this study 
  message("Subsetting datasheet")
  ms_ss<-subset_mastersheet(mastersheet_df, subsetlist)
  #check again to see the subset made sense
  master_xtabs(master_sheet = ms_ss,columns=c("sample_type",names(subsetlist)))
  
  message("STEP0 complete")

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
    
    cutadaptlist[[i]]<-system2(command = "cutadapt",args = c("-g", paste0("'Or1=",Pf,";max_error_rate=",ca.error.rate,"...",insect::rc(Pr),
                                                                          ";max_error_rate=",ca.error.rate,"'"), "-g",
                                          paste0("'Or2=",Pr,";max_error_rate=",ca.error.rate,"...",insect::rc(Pf),";max_error_rate=",ca.error.rate,"'")
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
  #get lengths of seqs in files
  datalist<-list()
  for(i in 1:length(files)){
    message(files[i])
    datalist[[i]]<-system2("seqkit", args=c("fx2tab","-l","-i","-n",files[i]),wait = T,stdout = T)
    datalist[[i]]<-read.table(text = datalist[[i]],header = T,sep = "\t")
    datalist[[i]][,5]<-files[i]
    datalist[[i]]<-datalist[[i]][,c(4:5)]
    colnames(datalist[[i]])<-c("length","file")
  }

  datalistdf<-as.data.frame(do.call(rbind,datalist))
  
  minlength<-ms_ss[match(ms_ss$barcode_name[i],ms_ss$barcode_name),"min_length"]
  maxlength<-ms_ss[match(ms_ss$barcode_name[i],ms_ss$barcode_name),"max_length"]
  
  datalistdf$InsideRange<-!(datalistdf$length<minlength | datalistdf$length>maxlength)
  
  lengthplot<-ggplot(datalistdf,mapping = aes(x = length,fill=InsideRange))+geom_histogram(binwidth = 10) +
    theme(axis.text.x=element_text(size=8,angle=45, hjust=1)) +
    scale_x_continuous(breaks = round(seq(0, max(datalistdf$length), by = 100),1)) +
    ggtitle("Range",subtitle = paste("Minlength=",minlength,"Maxlength=",maxlength))
  print(lengthplot)
  print(lengthplot + facet_wrap("file",scales="free_y"))
  
  plotfile1<-paste(gsub(", ",".",toString(unlist(subsetlist),)),"lengthplot.pdf",sep = ".")
  plotfile2<-paste(gsub(", ",".",toString(unlist(subsetlist),)),"lengthplot.facet.pdf",sep = ".")
  
  ggsave(filename = plotfile1,plot = lengthplot,device = "pdf", width = 15,height = 10)
  ggsave(filename = plotfile2, plot = lengthplot + facet_wrap("file",scales="free_y") ,device = "pdf", width = 15,height = 10)
  
  message(paste("Saving length plots to", plotfile1, "and", plotfile2))
  
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
    message(a-b," reads removed, from ",a, " (",round((a-b)/a*100,digits = 2)," %)")
  }
  message("STEP3 complete")
}

####################################################
#step 4 - cat files
if("step4" %in% stepstotake){  
  
  message("STEP4 - cat files")
  
  files<-paste0(ms_ss$barcode_name,".lenFilt.trimmed.fasta") 
  
  #put ss_sample_ids in headers
  for(i in 1:length(files)){
    file=files[i]
    ss_sample_id<-ms_ss[match(ms_ss$barcode_name[i],ms_ss$barcode_name),"ss_sample_id"]
    
    system2("seqkit", args=c("replace","-p", "'sampleid='", "-r", paste0("'ss_sample_id=",ss_sample_id," sampleid='"),file),wait = T,
            stdout = gsub(".fasta",".ids.fasta",file))
  }
  
  files<-paste0(ms_ss$barcode_name,".lenFilt.trimmed.ids.fasta") 
  
  #cat files
  catted_file<-paste(gsub(", ",".",toString(unlist(subsetlist),)),"lenFilt.trimmed.ids.fasta",sep = ".")
  system2(command = "cat",args = c(files),wait = T,stdout = catted_file)
}
    
####################################################
#step 4 
if("step4" %in% stepstotake){
  
  message("STEP4 - blast")
  
  if(!is.null(taxidlimit)){
     blast.status<-blast.min.bas(infastas = catted_file,refdb = refdb,blast_exec = blast_exec,
                                wait = T,taxidlimit = taxidlimit,taxidname = taxidname, ncbiTaxDir = ncbiTaxDir,task="blastn")
    } else blast.status<-blast.min.bas(infastas = startingfastas,refdb = refdb,blast_exec = blast_exec, wait = T, ncbiTaxDir = ncbiTaxDir,task="blastn")
  
    check.blasts(infastas = files,h = blast.status)
  
  message("STEP4 complete")
}
  
####################################################
#step 5 filter 
if("step5" %in% stepstotake){  
  
  message("STEP5")
  
  files<-paste0(ms_ss$barcode_name,".lenFilt.trimmed.blast.txt")  
  
  for(i in 1:length(files)){
    message(paste("filtering blast results for",files[i]))
    blastfile = files[i]
    out<-gsub(".blast.txt",".blast.filt.txt",files[i])
    filter.blast(blastfile = blastfile,ncbiTaxDir = ncbiTaxDir,out = out, top = top)
  }
  
  message("STEP5 complete")
  
}

####################################################
#step 6 BIN
if("step6" %in% stepstotake){ 
  
  message("STEP6 - Bin")
  
  files<-paste0(ms_ss$barcode_name,".lenFilt.trimmed.blast.filt.txt")  
  
  for(i in 1:length(files)){
    message(paste("binning filtered blast results for",files[i]))
    filtered_blastfile<-files[i]
    binfile<-gsub(".blast.filt.txt",".bins.txt",files[i])
    bin.blast2(filtered_blastfile = filtered_blastfile,ncbiTaxDir = ncbiTaxDir,
               obitaxdb = obitaxdb,out = binfile,spident = spident,gpident = gpident,fpident = fpident,abspident = abspident)
  }
  
  message("STEP6 complete")
  
} 

####################################################

message("Need to make otutabs")

  
####################################################
#step 7 make contributor files

if("step7" %in% stepstotake){  
  
  message("STEP7 - Contributor files")
  
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
  message("CHANGE SCRIPT TO ONLY TAKE RELEVANT FILES")
  
  files<-list.files(pattern = ".taxatable.tf.spliced.txt$")
  for(i in 1:length(files)){
    bas.krona.plot(files[i],KronaPath)
  }
  message("STEP14 complete")
  
}



  

  
