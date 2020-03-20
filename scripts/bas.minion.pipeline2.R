message("settings:
        ")
print(ls.str())

library(rlang,preffered_rlibsfolder)
library(tidyr,preffered_rlibsfolder)
library(processx)
library(dplyr)
library(ggplot2)
source(paste0(bastoolsDir,"master_functions.R"))
setwd(filesDir)
if(!is.null(catted_suffix)) {
  catted_file<-paste(gsub(", ",".",toString(unlist(subsetlist))),"lenFilt.trimmed.ids",catted_suffix,"fasta",sep = ".")
} else catted_file<-paste(gsub(", ",".",toString(unlist(subsetlist))),"lenFilt.trimmed.ids.fasta",sep = ".")


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
  
  t1<-Sys.time()
  
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
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP1 COMPLETE in ", t3, " min")
}

####################################################
#step 2 print length data
if("step2" %in% stepstotake){
  
  message("STEP2 - print length data")
  
  t1<-Sys.time()
  
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
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP2 COMPLETE in ", t3, " min")
}

####################################################
#step 3 size select
if("step3" %in% stepstotake){
  
  message("STEP3 - size select")
  
  t1<-Sys.time()
  
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
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP3 COMPLETE in ", t3, " min")
}

####################################################
#step 4 - cat files
if("step4" %in% stepstotake){  
  
  message("STEP4 - cat files")
  
  t1<-Sys.time()
  
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
  system2(command = "cat",args = c(files),wait = T,stdout = catted_file)
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP4 COMPLETE in ", t3, " min")
}
    
####################################################
#step 5 blast
if("step5" %in% stepstotake){
  
  message("STEP5 - blast")
  
  t1<-Sys.time()
  
  options=c("-word_size", 6, "-perc_identity", 50, "-qcov_hsp_perc", min_qcovs, "-gapopen", 0, "-gapextend", 2, "-reward", 1, "-penalty", -1)
  
  message("Using blastn for Minion data. Defaulting to ", options)
  
  if(!is.null(taxidlimit)){
     blast.status<-blast.min.bas(infastas = catted_file,refdb = refdb,blast_exec = blast_exec,
                                wait = T,taxidlimit = taxidlimit,taxidname = taxidname, ncbiTaxDir = ncbiTaxDir,task="blastn",
                                max_target_seqs = max_target_seqs,more = options)
    } else blast.status<-blast.min.bas(infastas = catted_file,refdb = refdb,blast_exec = blast_exec, wait = T, ncbiTaxDir = ncbiTaxDir,
                                       task="blastn",overWrite = T,max_target_seqs = max_target_seqs,more = options)
  
  check.blasts(infastas = catted_file,h = blast.status)
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP5 COMPLETE in ", t3, " min")
}
  
####################################################
#step 6 filter blast results 
if("step6" %in% stepstotake){  
  
  message("STEP6 - filter blast results")
  
  t1<-Sys.time()
  
  blastfile<-gsub(".fasta",".blast.txt",catted_file)
  message(paste("filtering blast results for",blastfile))
  out<-gsub(".blast.txt",".blast.filt.txt",blastfile)
  filter.blast2(blastfile = blastfile,ncbiTaxDir = ncbiTaxDir,out = out, min_qcovs = min_qcovs,max_evalue = max_evalue)
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP6 COMPLETE in ", t3, " min")
}

####################################################
#step 7 BIN
if("step7" %in% stepstotake){ 
  
  message("STEP7 - Bin")
  
  t1<-Sys.time()
  
  filtered_blastfile<-gsub(".fasta",".blast.filt.txt",catted_file)

  message(paste("binning filtered blast results for",filtered_blastfile))
  binfile<-gsub(".blast.filt.txt",".bins.txt",filtered_blastfile)
  
  bin.blast3(filtered_blastfile = filtered_blastfile,ncbiTaxDir = ncbiTaxDir,
               obitaxdb = obitaxdb,out = binfile,spident = spident,gpident = gpident,
             fpident = fpident,abspident = abspident,topS=topS,topG=topG,topF=topF,topAbs=topAbs)
  
  message("STEP7 complete")
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP7 COMPLETE in ", t3, " min")
  
} 

####################################################
#step 8 make otutab
  
if("step8" %in% stepstotake){  
    
  message("STEP8 - Making otutab")
  
  t1<-Sys.time()

  otutab<-data.table::fread(cmd = paste("grep", "'>'", catted_file),header = F,select = c(1,3),col.names = c("otu","ss_sample_id"),
                            data.table = F)
  otutab$ss_sample_id<-gsub("ss_sample_id=","",otutab$ss_sample_id)
  otutab$otu<-gsub(">","",otutab$otu)
  otutab$count<-1
  otutab<-as.data.frame(tidyr::pivot_wider(otutab,names_from = ss_sample_id,values_from=count,values_fill = list(count = 0)))
  
  write.table(otutab,gsub(".fasta",".otutab.txt",catted_file),append = F,quote = F,sep = "\t",row.names = F)
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP8 COMPLETE in ", t3, " min")
}

####################################################
#step 9 merge with otutab
  
if("step9" %in% stepstotake){  
  
  message("STEP9 - Merging bins and otutab into taxatab")
  
  t1<-Sys.time()
    
  otutab<-data.table::fread(gsub(".fasta",".otutab.txt",catted_file),data.table = F)
  taxon_input<-data.table::fread(gsub(".fasta",".bins.txt",catted_file) , sep = "\t",data.table = F)
  taxon_input$path<-paste(taxon_input$K,taxon_input$P,taxon_input$C,taxon_input$O,taxon_input$F,taxon_input$G,taxon_input$S,sep = ";")
  merged.table<-merge(otutab,taxon_input[,c("qseqid","path")],by.x = "otu",by.y = "qseqid",all = TRUE)
  merged.table$otu<-NULL
  merged.table<-merged.table %>% select(path,everything())
  taxatab<-aggregate(merged.table[,-1],by = list(merged.table$path),FUN=sum)
  colnames(taxatab)[1]<-"taxon"
  
  write.table(taxatab,gsub(".fasta",".taxatab.txt",catted_file),append = F,quote = F,sep = "\t",row.names = F)
    
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP9 COMPLETE in ", t3, " min")
}    
  
####################################################
#step10  make contributor files

if("step10" %in% stepstotake){  
  
  message("STEP10 - Contributor files")
  
  message(paste("making contributor file for",gsub(".fasta",".otutab.txt",catted_file)))
  
  t1<-Sys.time()
  
  filtered_blastfile<-gsub(".fasta",".blast.filt.txt",catted_file)
  binfile<-gsub(".blast.filt.txt",".bins.txt",filtered_blastfile)
  taxatab<-gsub(".fasta",".taxatab.txt",catted_file)
  
  check.low.res.df(
      filtered.taxatab = taxatab,filtered_blastfile = filtered_blastfile,
      binfile = binfile,disabledTaxaFile = NULL,spident = spident,gpident = gpident,fpident = fpident,abspident = abspident)
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP10 COMPLETE in ", t3, " min")
  
}

####################################################
#step 11 make krona plots

if("step11" %in% stepstotake){  
  
  t1<-Sys.time()
  
  message("STEP11 - make krona plots")

  bas.krona.plot(gsub(".fasta",".taxatab.txt",catted_file),KronaPath)

  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP11 COMPLETE in ", t3, " min")
  
}



  

  
