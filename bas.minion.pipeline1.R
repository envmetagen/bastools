
bas.minion.pipeline1<-function(demuxed_fstq_dir,mastersheet,expid,run){

#for unseperated files...
setwd(demuxed_fstq_dir)

  #cat 
  system2(command = "cat",args = c("*"),stdout = paste0(run,"_all.fastq"),wait = T)   
  
  #convert to fasta
  system2(command = "sed",args = c("-n", "'1~4s/^@/>/p;2~4p'",paste0(run,"_all.fastq")),wait=T,stdout= paste0(run,"_all.fasta"))
  
  #reformat for obitools
  system2(command = "sed",args = c("'s/ /; /g;s/; / /1'",paste0(run,"_all.fasta")),wait=T,stdout= paste0(run,"_all_rf.fasta"))
  
  #split by barcode
  system2(command = "obisplit", args = c("-t", "barcode", paste0(run,"_all_rf.fasta")),wait=T)
  #these files should be the ones for DMP
  message("should replace barcodes with sample names here and change filenames for DMP")
          
  message("not removing primers in this pipeline")
  
  #size select, for each fragment, I checked an seqs appear to have primers plus one base (at each end)
    #list barcodes in each frag in run
    #this run
    mastersheetrun<-mastersheet[mastersheet$experiment_id==expid,]
    #used bcs each frag  
    primersets<-unique(mastersheetrun$Primer_set)
    a<-list()
    for(i in 1:length(primersets)){
      a[[i]]<-mastersheetrun[mastersheetrun$Primer_set==primersets[i],"barcode_id"]
    }
    names(a)<-primersets
    
    #cat each frag
    for(i in 1:length(a)){
      system2(command = "cat",args = c(paste0(a[[i]]$barcode_id,".fasta")),wait = T,stdout = paste0(run,"_",names(a)[i],".fasta")) 
    }
    
    #size calc
    fraglengths<-data.frame("primerset"=primersets)
    for(i in 1:length(primersets)){
    thisprimer<-match(primersets[i],mastersheetrun$Primer_set)
    minL<-mastersheetrun$Min_length[thisprimer]+nchar(mastersheetrun$Primer_F[thisprimer])+
      nchar(mastersheetrun$Primer_R[thisprimer])
    maxL<-mastersheetrun$Max_length[thisprimer]+nchar(mastersheetrun$Primer_F[thisprimer])+
      nchar(mastersheetrun$Primer_R[thisprimer])
    fraglengths[match(primersets[i],fraglengths$primerset),2]<-minL
    fraglengths[match(primersets[i],fraglengths$primerset),3]<-maxL
    }
    
    #add length to fasta
    for(i in 1:length(names(a))){
    system2(command = "obiannotate",args = c("--length", paste0(run,"_",names(a)[i],".fasta")),wait = T,
            stdout =  paste0(run,"_",names(a)[i],"_wlen.fasta"))
    }
    
    #obigrep lengths
    for(i in 1:length(fraglengths$primerset)){
      system2(command = "obigrep",args = c("-l",fraglengths[i,2],"-L",fraglengths[i,3],
                                           paste0(run,"_",names(a)[i],"_wlen.fasta")),wait = T,
              stdout =  paste0(run,"_",names(a)[i],"_l",fraglengths[i,2],"L",fraglengths[i,3],"_wlen.fasta"))
    }
    
    #obitab all
    for(i in 1:length(fraglengths$primerset)){
      system2(command = "obitab",args = c(paste0(run,"_",names(a)[i],
                                                 "_l",fraglengths[i,2],"L",fraglengths[i,3],"_wlen.fasta")),wait = T,
              stdout =  paste0(run,"_",names(a)[i],"_l",fraglengths[i,2],"L",fraglengths[i,3],"_wlen.tab"))
    }
  
  # #dereplicate - not worth it, doesnt reduce anything!
 
    
}
    
    
