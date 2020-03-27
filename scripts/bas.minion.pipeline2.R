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

catted_file<-paste(gsub(", ",".",toString(unlist(subsetlist))),"lenFilt.trimmed.ids.fasta",sep = ".")

if(!is.null(catted_suffix)) catted_file<-gsub(".fasta",paste0(".",catted_suffix,".fasta"),catted_file)
  
if(polished) catted_file<-gsub(".fasta",".pol.fasta",catted_file)

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
  
  if(polished) stop("If using polished reads should start from step 2")
  
  t1<-Sys.time()
  
  files<-paste0(origfilesDir,ms_ss$barcode_name,".fastq.gz")
  
  cutadaptlist<-list()
  
  for(i in 1:length(files)){
    file<-files[i]
    message("Running cutadapt on ",file)
    if(!file %in% gsub("//","/",list.files(path = origfilesDir,full.names = T))) stop("File not found")
    barcode<-ms_ss$barcode_name[i]
    Pf<-ms_ss[match(ms_ss$barcode_name[i],ms_ss$barcode_name),"Primer_F"]
    Pr<-ms_ss[match(ms_ss$barcode_name[i],ms_ss$barcode_name),"Primer_R"]
    
    message("Pf=",Pf)
    message("Pr=",Pr)
    
    out<-gsub(".fastq.gz",".trimmed.fastq",paste0("barcode",do.call(rbind,stringr::str_split(file,"barcode"))[,2]))
    
    cutadaptlist[[i]]<-system2(command = "cutadapt",args = c("--discard-untrimmed","-g", paste0("'Or1=",Pf,";max_error_rate=",ca.error.rate,"...",insect::rc(Pr),
                                                                          ";max_error_rate=",ca.error.rate,"'"), "-g",
                                          paste0("'Or2=",Pr,";max_error_rate=",ca.error.rate,"...",insect::rc(Pf),";max_error_rate=",ca.error.rate,"'")
                                          ,"-o", out,file),wait = T,stderr = T,stdout = T)
    
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
  
  #format polished results
  if(polished){
    message("assumes polished reads are in results.fasta.gz")
    ##remove seqs with "adapter=no_adapter"
    system2("seqkit", args=c("grep","-r","-n","-v","-p","'adapter=no_adapter'","results.fasta.gz"),wait = T,stdout = "results.fasta.temp.gz")
    
    #split polished by barcode
    for(i in 1:length(ms_ss$barcode_name)){
      system2("seqkit", args=c("grep","-r","-n","-p",ms_ss$barcode_name[i],"results.fasta.temp.gz"),wait = T,
              stdout = paste0(ms_ss$barcode_name[i],".trimmed.fastq"))
    }
  }

  
  files<-paste0(ms_ss$barcode_name,".trimmed.fastq")
  #get lengths of seqs in files
  datalist<-list()
  for(i in 1:length(files)){
    message(files[i])
    datalist[[i]]<-system2("seqkit", args=c("fx2tab","-l","-i","-n",files[i]),wait = T,stdout = T)
    datalist[[i]]<-read.table(text = datalist[[i]],header = T,sep = "\t")
    datalist[[i]][,5]<-files[i]
    
    if(polished){
      datalist[[i]][,6]<-as.numeric(do.call(rbind,stringr::str_split(do.call(rbind,stringr::str_split(datalist[[i]][,1],"size="))[,2],":"))[,1])
      
      rep.row<-function(x,n){
        matrix(rep(x,each=n),nrow=n)
      }
      expanded_lengths<-list()
      for(j in 1:nrow(datalist[[i]])){
        expanded_lengths[[j]]<-rep.row(datalist[[i]][j,4],n = datalist[[i]][j,6])
      }
      expanded_lengthsdf<-as.data.frame(unlist(expanded_lengths))
      expanded_lengthsdf$file<-files[i]
      datalist[[i]]<-expanded_lengthsdf
    } else datalist[[i]]<-datalist[[i]][,c(4:5)]
    
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
  
  plotfile1<-paste(gsub(", ",".",toString(unlist(subsetlist))),"lengthplot.pdf",sep = ".")
  plotfile2<-paste(gsub(", ",".",toString(unlist(subsetlist))),"lengthplot.facet.pdf",sep = ".")
  
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
    
    system2("seqkit", args=c("seq","-M",maxlength,"-m",minlength,file),wait = T, stdout = gsub(".trimmed.fasta",".lenFilt.trimmed.fasta",file)
                ,stderr = F)
    
    #count seqs original
    a<-system2("seqkit", args=c("stats","-T",file),wait = T,stderr = T,stdout = T)
    a<-read.table(text = a,header = T)$num_seqs
    #count seqs after
    b<-system2("seqkit", args=c("stats","-T",gsub(".trimmed.fasta",".lenFilt.trimmed.fasta",file)),wait = T,stderr = T,stdout = T)
    if(length(grep("DNA",b))>0) b<-read.table(text = b,header = T)$num_seqs else b=0
    message(a-b," sequences removed, from ",a, " (",round((a-b)/a*100,digits = 2)," %)")
  }
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP3 COMPLETE in ", t3, " min")
}
  
####################################################
#step 4 plot step counts
  
if("step4" %in% stepstotake){  
  
    t1<-Sys.time()
    
    message("STEP4 - calculate step counts")
    
    #count seqs in starting files
    files<-paste0(origfilesDir,ms_ss$barcode_name,".fastq.gz")
    
    #table for storing counts
    count.step.table<-data.frame(file=files)

    datalist<-list()
    for(i in 1:length(files)){
      datalist[[i]]<-system2("seqkit", args=c("stats","-T",files[i]),wait = T,stdout = T)
      datalist[[i]]<-read.table(text = datalist[[i]],header = T,sep = "\t")
    }
    
    starting.counts<-do.call(rbind,datalist)
    
    count.step.table<-merge(count.step.table,starting.counts[,c("file","num_seqs")], by = "file",all.x = T)
    colnames(count.step.table)[2]<-"starting.counts"
    
    starting.counts<-sum(starting.counts$num_seqs)  
    
    #count after cutadapt
    files<-paste0(ms_ss$barcode_name,".trimmed.fastq")
    datalist<-list()
    
    if(polished==F){
      for(i in 1:length(files)){
        datalist[[i]]<-system2("seqkit", args=c("stats","-T",files[i]),wait = T,stdout = T)
        datalist[[i]]<-read.table(text = datalist[[i]],header = T,sep = "\t")
      }
      
      after.cutadapt<-do.call(rbind,datalist)
      
      count.step.table$file_map<-gsub(".fastq.gz",".trimmed.fastq",paste0("barcode",do.call(rbind,stringr::str_split(count.step.table$file,"barcode"))[,2]))
      count.step.table<-merge(count.step.table,after.cutadapt[,c("file","num_seqs")], by.x = "file_map",by.y="file",all.x = T)
      colnames(count.step.table)<-gsub("num_seqs","after.cutadapt",colnames(count.step.table))
      
      after.cutadapt<-sum(after.cutadapt$num_seqs) 
      
      } else {
          count.polished<-function(fastx.file){
            counts<-system2("seqkit", args=c("fx2tab","-l","-i","-n",fastx.file),wait = T,stdout = T)
            if(length(counts)!=0) {
              counts<-read.table(text = counts,header = F,sep = "\t")
              counts$size<-as.numeric(do.call(rbind,stringr::str_split(do.call(rbind,stringr::str_split(counts[,1],"size="))[,2],":"))[,1])
              sum(counts$size)
            } else 0
        }
        
        for(i in 1:length(files)){
          datalist[[i]]<-as.data.frame(count.polished(files[i]))
          datalist[[i]]$file<-files[i]
          colnames(datalist[[i]])[1]<-"num_seqs"
        }
          after.polishing<-do.call(rbind,datalist)
          
          count.step.table$file_map<-gsub(".fastq.gz",".trimmed.fastq",paste0("barcode",do.call(rbind,stringr::str_split(count.step.table$file,"barcode"))[,2]))
          count.step.table<-merge(count.step.table,after.polishing[,c("file","num_seqs")], by.x = "file_map",by.y="file",all.x = T)
          colnames(count.step.table)<-gsub("num_seqs","after.polishing",colnames(count.step.table))
          
          after.polishing<-sum(after.polishing$num_seqs)
      }
   
    #after length filtering
    datalist<-list()
    files<-paste0(ms_ss$barcode_name,".lenFilt.trimmed.fasta")
    
    if(polished==F){
      for(i in 1:length(files)){
          datalist[[i]]<-system2("seqkit", args=c("stats","-T",files[i]),wait = T,stdout = T)
          datalist[[i]]<-read.table(text = datalist[[i]],header = T,sep = "\t")
      } 
      
      after.length.filt<-do.call(rbind,datalist)
      
      count.step.table$file_map<-gsub(".trimmed.fastq",".lenFilt.trimmed.fasta",count.step.table$file_map)
      count.step.table<-merge(count.step.table,after.length.filt[,c("file","num_seqs")], by.x = "file_map",by.y="file",all.x = T)
      colnames(count.step.table)<-gsub("num_seqs","after.length.filt",colnames(count.step.table))
      
      after.length.filt<-sum(after.length.filt$num_seqs) 
      
    } else {
        for(i in 1:length(files)){
          datalist[[i]]<-as.data.frame(count.polished(files[i]))
          datalist[[i]]$file<-files[i]
          colnames(datalist[[i]])[1]<-"num_seqs"
        }
      
        after.length.filt<-do.call(rbind,datalist)
        
        count.step.table$file_map<-gsub(".trimmed.fastq",".lenFilt.trimmed.fasta",count.step.table$file_map)
        count.step.table<-merge(count.step.table,after.length.filt[,c("file","num_seqs")], by.x = "file_map",by.y="file",all.x = T)
        colnames(count.step.table)<-gsub("num_seqs","after.length.filt",colnames(count.step.table))
        
        after.length.filt<-sum(after.length.filt$num_seqs)
      }
     
    if(polished){
      counts<-data.frame(step=c("Starting fastq files","After polishing", "After length filtering"),
                         reads=c(starting.counts,after.polishing,after.length.filt))
    } else {
      counts<-data.frame(step=c("Starting fastq files","After primer trimming", "After length filtering"),
                       reads=c(starting.counts,after.cutadapt,after.length.filt))
    }
    
    counts$step<-factor(counts$step,levels = counts$step)
    
    countplot<-ggplot(counts,aes(x = step,y=reads))+
      theme(axis.text.x=element_text(size=8,angle=45, hjust=1))+
      ggtitle(gsub(".fasta","",catted_file))+
      geom_bar(stat = "identity")+ 
      scale_y_continuous(labels = scales::comma)
    
    ggsave(filename = gsub(".fasta",".stepCounts.pdf",catted_file),plot = countplot,device = "pdf", width = 15,height = 10)
    
    count.step.table<-count.step.table[,-match("file_map",colnames(count.step.table))]
    write.table(count.step.table,file = gsub(".fasta",".stepCounts.by.sample.txt",catted_file),append = F,row.names = F,quote = F,sep = "\t")
    
    t2<-Sys.time()
    t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
    message("STEP4 COMPLETE in ", t3, " min")
  }
 
####################################################
#step 5 - cat files
if("step5" %in% stepstotake){  
  
  message("STEP5 - cat files")
  
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
  message("STEP5 COMPLETE in ", t3, " min")
}
  
####################################################
#step 6 blast
if("step6" %in% stepstotake){
  
  message("STEP6 - blast")
  
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
  message("STEP6 COMPLETE in ", t3, " min")
}
  
####################################################
#step 7 plot hit pidents blast results 
if("step7" %in% stepstotake){  
    
 message("STEP7 - plot pidents blast results")
    
 t1<-Sys.time()
 
 blastfile<-gsub(".fasta",".blast.txt",catted_file)
 
 a<-report.blast.maxmin(blastfile)
 
 write.table(a,file = gsub(".txt",".minmax.hits.txt",blastfile),append = F,quote = F,row.names = F,sep = "\t")
 
 b<-ggplot(a,aes(x=minmax,y=pident)) + geom_violin() 
  
 ggsave(filename = gsub(".txt",".minmax.hits.pdf",blastfile),b,device = "pdf", width = 15,height = 10)
 
 rm(b)
 
 t2<-Sys.time()
 t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
 message("STEP7 COMPLETE in ", t3, " min")
 
}
  
####################################################
#step 8 filter blast results 
if("step8" %in% stepstotake){  
  
  message("STEP8 - filter blast results")
  
  t1<-Sys.time()
  
  blastfile<-gsub(".fasta",".blast.txt",catted_file)
  message(paste("filtering blast results for",blastfile))
  out<-gsub(".blast.txt",".blast.filt.txt",blastfile)
  filter.blast2(blastfile = blastfile,ncbiTaxDir = ncbiTaxDir,out = out, min_qcovs = min_qcovs,max_evalue = max_evalue)
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP8 COMPLETE in ", t3, " min")
}

####################################################
#step 9 BIN
if("step9" %in% stepstotake){ 
  
  message("STEP9 - Bin")
  
  t1<-Sys.time()
  
  filtered_blastfile<-gsub(".fasta",".blast.filt.txt",catted_file)

  message(paste("binning filtered blast results for",filtered_blastfile))
  binfile<-gsub(".blast.filt.txt",".bins.txt",filtered_blastfile)
  
  bin.blast3(filtered_blastfile = filtered_blastfile,ncbiTaxDir = ncbiTaxDir,
               out = binfile,spident = spident,gpident = gpident,
             fpident = fpident,abspident = abspident,topS=topS,topG=topG,topF=topF,topAbs=topAbs)
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP9 COMPLETE in ", t3, " min")
  
} 

####################################################
#step 10 make otutab
if("step10" %in% stepstotake){  
    
  message("STEP10 - Making otutab")
  
  t1<-Sys.time()
  
  if(polished==F){

    otutab<-data.table::fread(cmd = paste("grep", "'>'", catted_file),header = F,select = c(1,3),col.names = c("otu","ss_sample_id"),
                              data.table = F)
    otutab$ss_sample_id<-gsub("ss_sample_id=","",otutab$ss_sample_id)
    otutab$otu<-gsub(">","",otutab$otu)
    otutab$count<-1
    otutab<-as.data.frame(tidyr::pivot_wider(otutab,names_from = ss_sample_id,values_from=count,values_fill = list(count = 0)))
  }  else {
      otutab<-data.table::fread(cmd = paste("grep", "'>'", catted_file),header = F, data.table = F)
      otutab$otu<-otutab[,1]
      otutab$otu<-gsub(">","",otutab$otu)
      otutab$ss_sample_id<-do.call(rbind,stringr::str_split(otutab$otu,"ss_sample_id="))[,2]
      
      otutab$count<-as.numeric(do.call(rbind,stringr::str_split(do.call(rbind,stringr::str_split(otutab[,2],"size="))[,2],":"))[,1])
      otutab<-as.data.frame(tidyr::pivot_wider(otutab[,c(-1,-2)],names_from = ss_sample_id,values_from=count,values_fill = list(count = 0)))
  }
  
    write.table(otutab,gsub(".fasta",".otutab.txt",catted_file),append = F,quote = F,sep = "\t",row.names = F)
    
    t2<-Sys.time()
    t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
    message("STEP10 COMPLETE in ", t3, " min")
}

####################################################
#step 11 merge with otutab
  
if("step11" %in% stepstotake){  
  
  message("STEP11 - Merging bins and otutab into taxatab")
  
  t1<-Sys.time()
    
  otutab<-data.table::fread(gsub(".fasta",".otutab.txt",catted_file),data.table = F)
  taxon_input<-data.table::fread(gsub(".fasta",".bins.txt",catted_file) , sep = "\t",data.table = F)
  taxon_input$path<-paste(taxon_input$K,taxon_input$P,taxon_input$C,taxon_input$O,taxon_input$F,taxon_input$G,taxon_input$S,sep = ";")
  merged.table<-merge(otutab,taxon_input[,c("qseqid","path")],by.x = "otu",by.y = "qseqid",all = TRUE)
  merged.table$otu<-NULL
  merged.table<-merged.table %>% select(path,everything())
  merged.table$path[is.na(merged.table$path)]<-"no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits"
  taxatab<-aggregate(merged.table[,-1],by = list(merged.table$path),FUN=sum)
  colnames(taxatab)[1]<-"taxon"
  
  write.table(taxatab,gsub(".fasta",".taxatab.txt",catted_file),append = F,quote = F,sep = "\t",row.names = F)
    
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP11 COMPLETE in ", t3, " min")
}    
  
####################################################
#step 12 apply taxon filter
  
if("step12" %in% stepstotake){  
    
    message("STEP12 - taxon filter")

    taxon.filter.solo(gsub(".fasta",".taxatab.txt",catted_file),taxonpc)
    
    message("STEP12 complete")
    
}
  
####################################################
#step13  make contributor files

if("step13" %in% stepstotake){  
  
  message("STEP13 - Contributor files")
  
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
  message("STEP13 COMPLETE in ", t3, " min")
  
}

####################################################
#step 14 make krona plots

if("step14" %in% stepstotake){  
  
  t1<-Sys.time()
  
  message("STEP14 - make krona and blast.hit plots")

  bas.krona.plot(gsub(".fasta",".taxatab.txt",catted_file),KronaPath)

  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  message("STEP14 COMPLETE in ", t3, " min")
  
}

