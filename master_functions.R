taxon.filter.solo<-function(files,filterpc=0.1){
  #remove detections where less than x% of total reads for taxon
  taxatables<-list()
  for(i in 1:length(files)){
    taxatables[[i]]<-data.table::fread(files[i],data.table = F)
    rownames(taxatables[[i]])<-taxatables[[i]]$taxon
    taxatables[[i]]$taxon=NULL
    taxatables[[i]]<-taxatables[[i]][rowSums(taxatables[[i]])!=0,]
    taxatab.PCS<-sweep(taxatables[[i]], MARGIN = 1, STATS = rowSums(taxatables[[i]]), FUN = "/")*100
    taxatab.PCS[taxatab.PCS<filterpc]<-0 
    taxatables[[i]][taxatab.PCS==0]<-0
    taxatables[[i]]<-taxatables[[i]][rowSums(taxatables[[i]])!=0,]
    taxatables[[i]]$taxon<-rownames(taxatables[[i]])
    taxatables[[i]]<-taxatables[[i]][,c(length(colnames(taxatables[[i]])),1:(length(colnames(taxatables[[i]]))-1))]
    write.table(taxatables[[i]],gsub("taxatable.txt","taxatable.tf.txt",files[i]),row.names = F,quote = F,sep = "\t")
  }
}

taxon.filter.solo.df<-function(taxatab,taxonpc=0.1){
  message("Applying taxon_pc filter.")
  taxatab2<-taxatab[,-1]
  a<-sum(taxatab2)
  taxatab.PCS<-sweep(taxatab2, MARGIN = 1, STATS = rowSums(taxatab2), FUN = "/")*100
  taxatab.PCS[taxatab.PCS<taxonpc]<-0 
  taxatab2[taxatab.PCS==0]<-0
  b<-sum(taxatab2)
  taxatab3<-cbind(taxon=taxatab$taxon,taxatab2)
  message(paste(a-b,"reads removed"))
  taxatab4<-rm.0readtaxSam(taxatab3)
}

sample.filter.solo<-function(taxatab,samplepc=0.1){
  message("Applying sample_pc filter. Note: this removes samples with no reads")
  taxatab2<-taxatab[,2:length(colnames(taxatab))]
  taxon<-taxatab$taxon
  taxatab<-taxatab2[,!colSums(taxatab2)==0]
  taxatab.PCS<-sweep(taxatab, MARGIN = 2, STATS = colSums(taxatab), FUN = "/")*100
  taxatab.PCS[taxatab.PCS<samplepc]<-0
  
  a<-sum(taxatab)
  taxatab[taxatab.PCS==0]<-0
  b<-sum(taxatab)
  message(paste(a-b,"reads removed"))
  
  taxatab.out<-cbind(taxon,taxatab)
  taxatab.out<-rm.0readtaxSam(taxatab.out)
}

rm.0readtaxSam<-function(taxatab){
  taxatab2<-taxatab[rowSums(taxatab[,-1])!=0,]
  taxatab3<-cbind(taxon=taxatab2[,1],taxatab2[,-1][,colSums(taxatab2[,-1])!=0])
}

rm.0readOTUSam<-function(taxatab){
  taxatab2<-taxatab[rowSums(taxatab[,-1])!=0,]
  taxatab3<-cbind(OTU=taxatab2$OTU,taxatab2[,-1][,colSums(taxatab2[,-1])!=0])
}

bas.merge.taxatabs<-function(taxatabs){
  require(tidyverse)
  
  taxatabs.list<-list()
  counts<-data.frame(file=taxatabs,reads=0,taxa=0,samples=0)
  
  for(i in 1:length(taxatabs)){
    taxatabs.list[[i]]<-data.table::fread(taxatabs[i],sep = "\t",data.table = F)
    counts[i,2]<-sum(taxatabs.list[[i]][,-1])
    counts[i,3]<-length(taxatabs.list[[i]][,1])
    counts[i,4]<-length(colnames(taxatabs.list[[i]][,-1]))
    
    message(counts[i,1])
    message(paste0("reads: ",counts[i,2],", taxa: ",counts[i,3],", samples: ",counts[i,4]))
  }
  
  all.taxatabs<-taxatabs.list %>% purrr::reduce(full_join, by = "taxon")
  #remove NAs
  all.taxatabs[is.na(all.taxatabs)]<-0
  
  #check sums
  a<-sum(all.taxatabs[,-1])
  if(a==sum(counts$reads)) message("Read counts all good") else stop("Read counts do not match")
  
  message("merged taxatable")
  message(paste0("reads: ",a,", taxa: ",length(all.taxatabs[,1]),", samples: ",length(colnames(all.taxatabs[,-1]))))
  
  return(all.taxatabs)
}


filter.dxns<-function(taxatab,filter_dxn=50){
  taxatab2<-taxatab[,-1]
  reads1<-sum(taxatab2)
  dxns1<-sum(taxatab2>0)
  taxatab2[taxatab2<filter_dxn] <- 0
  reads2<-sum(taxatab2)
  dxns2<-sum(taxatab2>0)
  taxatab2<-cbind(taxon=taxatab$taxon,taxatab2)  
  taxatab2<-rm.0readtaxSam(taxatab2)
  message(paste("reads removed:",reads1-reads2,"from",reads1, "; detections removed:",dxns1-dxns2,"from",dxns1))
  return(taxatab2)
}

count.dxns.by.taxon<-function(taxatab){
  taxatab2<-taxatab[,-1]
  return(data.frame(taxon=taxatab$taxon,n.samples=as.numeric(apply(taxatab2,1,function(x) sum(x>0)))))
}

sum.reads.by.taxon<-function(taxatab){
  taxatab2<-taxatab[,-1]
  return(data.frame(taxon=taxatab$taxon,total.reads=as.numeric(apply(taxatab2,1,function(x) sum(x)))))
}


range.dxns.by.taxon<-function(taxatab){
  taxatab2<-taxatab[,-1]
  taxatab2[taxatab2==0]<-NA
  rangetab<-apply(taxatab2,1,function(x) range(x,na.rm = T))
  return(data.frame(taxon=taxatab$taxon,low=rangetab[1,],high=rangetab[2,]))
}

summary.dxns.by.taxon<-function(taxatab){
  taxatab2<-count.dxns.by.taxon(taxatab)
  taxatab3<-range.dxns.by.taxon(taxatab)
  taxatab4<-sum.reads.by.taxon(taxatab)
  taxatab5<-merge(taxatab2,taxatab4,by="taxon")
  taxatab6<-merge(taxatab5,taxatab3,by="taxon")
  return(taxatab6)
}

write.taxatab<-function(taxatab,out){
  write.table(taxatab,file = out,append = F,quote = F,sep = "\t",row.names = F)
}

#group taxa
bas.group.taxa<-function(taxatab,taxon, jointo){
  
  taxontable1<-taxatab[taxatab$taxon==taxon,]
  taxontable2<-taxatab[taxatab$taxon==jointo,]
  taxontable3<-rbind(taxontable1,taxontable2)
  taxontable4<-cbind(taxon=jointo,as.data.frame(t(colSums(taxontable3[,-1]))))
  
  taxatab<-taxatab[!taxatab$taxon==jointo,]
  taxatab<-taxatab[!taxatab$taxon==taxon,]
  
  taxatab<-rbind(taxatab,taxontable4)
}

negs.stats<-function(taxatab,ms_ss,real,ex_hominidae=T,printnegs=T){
  message("Ignoring the following taxa: NA;NA;NA;NA;NA;NA;NA & no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits")
  if(ex_hominidae) message(" & Hominidae")
  
  if(is.null(ms_ss$sample_type)) stop("No column called sample_type")
  
  taxatab2<-taxatab[taxatab$taxon!="NA;NA;NA;NA;NA;NA;NA",]
  taxatab2<-taxatab2[taxatab2$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]
  if(ex_hominidae) taxatab2<-taxatab2[-grep("Hominidae",taxatab2$taxon),]
  
  #find negatives with reads
  negs<-ms_ss[!ms_ss$sample_type %in% real,c("ss_sample_id","sample_type")]
  reads.in.negs<-as.data.frame(colSums(taxatab2[colnames(taxatab2) %in% negs$ss_sample_id]))
  colnames(reads.in.negs)<-"reads"
  reads.in.negs$ss_sample_id<-rownames(reads.in.negs)
  read.in.negs<-reads.in.negs[reads.in.negs$reads!=0,]
  
  if(length(reads.in.negs$ss_sample_id)>0){
    
    #taxatable for negs
    taxatab.negs<-cbind(taxon=taxatab2$taxon,taxatab2[colnames(taxatab2) %in% read.in.negs$ss_sample_id])
    taxatab.negs<-taxatab.negs[rowSums(taxatab.negs[,-1,drop=FALSE])!=0,]
    
    #taxatab by sample type
    taxatab.negs.list<-list()
    for(i in 1:length(unique(negs$sample_type))){
      taxatab.negs.list[[i]]<-cbind(taxon=taxatab.negs$taxon,taxatab.negs[colnames(taxatab.negs) %in% negs$ss_sample_id[negs$sample_type==unique(negs$sample_type)[i]]])
      taxatab.negs.list[[i]]<-taxatab.negs.list[[i]][rowSums(taxatab.negs.list[[i]][,-1,drop=FALSE])!=0,]
      names(taxatab.negs.list)[i]<-unique(negs$sample_type)[i]
    }
    #summary sentence
    for(i in 1:length(unique(negs$sample_type))){
      message(gsub("-1","0",paste("from",length(negs$ss_sample_id[negs$sample_type==unique(negs$sample_type)[i]]), unique(negs$sample_type)[i],"samples", 
                                  length(colnames(taxatab.negs.list[grep(unique(negs$sample_type)[i],names(taxatab.negs.list))][[1]]))-1, "contained reads")))
      
      if(printnegs==T){
      if(length(colnames(taxatab.negs.list[grep(unique(negs$sample_type)[i],names(taxatab.negs.list))][[1]]))>0){
        print(taxatab.negs.list[grep(unique(negs$sample_type)[i],names(taxatab.negs.list))][[1]])
        }
      }
    }
    
    return(taxatab.negs.list)
  } else message("No negatives found with reads")
  
}

#remove contaminant taxa from samples, based on occurrences in negatives
remove.contaminant.taxa<-function(master_sheet,taxatab,negatives,group.code,printcontaminations=T){
  
  totalstart<-sum(taxatab[,-1,drop=F])
  
  for(i in 1:length(negatives)){
    
    negdf<-negatives[[i]]
    
    negative_type<-names(negatives)[i]
    
    if(length(negdf)!=0){
      
      for(i in 2:length(colnames(negdf))){
        neg<-colnames(negdf)[i]
        negtaxa<-negdf[,c("taxon",neg)]
        sumnegtaxa<-sum(negtaxa[,-1,drop=F])
        negtaxa<-as.character(negtaxa[rowSums(negtaxa[,-1,drop=F])>0,"taxon"])
        
        
        #get site of a sample
        group.id<-master_sheet[grep(neg,master_sheet[,"ss_sample_id"]),group.code]
        #get other samples in site
        group.samples<-master_sheet[master_sheet[,group.code]==group.id,"ss_sample_id"]
        group.samples<-group.samples[!is.na(group.samples)]
        
        #temp
        #negtaxa<-as.character(taxatab$taxon)
        
        message(paste("Based on",negative_type,neg, ", Removing detections of"))
        print(negtaxa)
        message("from")
        print(group.samples)
        message(paste("which all belong to",group.code,":",group.id))
        
        sumb4<-sum(taxatab[,-1])
        
        #put taxon counts to 0
        contaminations<-cbind(taxon=taxatab$taxon[taxatab$taxon %in% negtaxa],taxatab[taxatab$taxon %in% negtaxa,colnames(taxatab) %in% group.samples])
        contaminations<-cbind(taxon=contaminations[,1],contaminations[,-1,drop=F][,colSums(contaminations[,-1,drop=F])>0,drop=F])
        taxatab[taxatab$taxon %in% negtaxa,colnames(taxatab) %in% group.samples]<-0
        
        sumafter<-sum(taxatab[,-1])
        
        message(paste(sumb4-sumafter,"reads removed from",sum(contaminations[,-1,drop=F]>0),"detection(s) in",
                      length(colnames(contaminations[,-1,drop=F])),"sample(s), of which",sumnegtaxa,"were in the negative. See table below for details:"))
        if(printcontaminations==T) {
          print(contaminations)
          }else(message("Not printing contamination table"))
        }
    }
  }
  taxatab<-rm.0readtaxSam(taxatab)
  totalend<-sum(taxatab[,-1,drop=F])
  message(paste("***********A total of", totalstart-totalend, "reads removed"))
  return(taxatab)
}

#keep only species-level assignments
keep.only.spL.assigns<-function(taxatab){
  #if(level=="species"){
  taxatab<-taxatab[-grep(";NA;NA$",taxatab$taxon),]
  taxatab<-taxatab[-grep(";NA$",taxatab$taxon),]
  taxatab<-rm.0readtaxSam(taxatab)
  return(taxatab)
}

adonis.bas<-function(taxatab,master_sheet,factor1,samLevel="ss_sample_id",stratum=NULL){
  taxatab<-rm.0readtaxSam(taxatab)
  taxatab2<-binarise.taxatab(taxatab)
  distance_matrix<-taxatab2bray(taxatab2)
  
  if(samLevel=="ss_sample_id") {
    if(!"ss_sample_id" %in% colnames(master_sheet)) stop("No column called ss_sample_id")
    master_sheet2<-master_sheet[master_sheet$ss_sample_id  %in% colnames(taxatab),]
  }
  
  if(samLevel=="Sample_Name") {
    if(!"Sample_Name" %in% colnames(master_sheet)) stop("No column called ss_sample_id")
    master_sheet2<-master_sheet[master_sheet$Sample_Name  %in% colnames(taxatab),]
    master_sheet2<-master_sheet2[!duplicated(master_sheet2$Sample_Name),]
  }
  
  
  if(!is.null(stratum)){
    #model using adonis
    vegan::adonis(distance_matrix~master_sheet2[,factor1],by="term",method="bray", data = master_sheet2,
                  strata = master_sheet2[,stratum],permutations = 10000)
  } else {
    #model using adonis
    vegan::adonis(distance_matrix~master_sheet2[,factor1],by="term",method="bray", data = master_sheet2,permutations = 10000)
  }
  
}

add.lineage.df<-function(df,ncbiTaxDir){
  if(is.null(df$taxids)) {stop("No column called taxids")}
  df$taxids<-as.integer(as.character(df$taxids)) 
  #write taxids to file
  taxids_fileA<-paste0("taxids",as.numeric(Sys.time()),".txt")
  write.table(unique(df$taxids),file = taxids_fileA,row.names = F,col.names = F,quote = F)
  
  #get taxonomy from taxids and format in 7 levels
  taxids_fileB<-paste0("taxids",as.numeric(Sys.time()),".txt")
  system2(command = "taxonkit",args =  c("lineage","-r",taxids_fileA,"-c","--data-dir",ncbiTaxDir)
          ,stdout = taxids_fileB,stderr = "",wait = T)
  taxids_fileC<-paste0("taxids",as.numeric(Sys.time()),".txt")
  system2(command = "taxonkit",args =  c("reformat",taxids_fileB,"-i",3,"--data-dir",ncbiTaxDir)
          ,stdout = taxids_fileC,stderr = "",wait = T)
  lineage<-as.data.frame(data.table::fread(file = taxids_fileC,sep = "\t"))
  colnames(lineage)<-gsub("V1","taxids",colnames(lineage))
  colnames(lineage)<-gsub("V2","new_taxids",colnames(lineage))
  colnames(lineage)<-gsub("V5","path",colnames(lineage))
  
  #merge with df
  message("replacing taxids with updated taxids. Saving old taxids in old_taxids.")
  df<-merge(df,lineage[,c(1,2,5)],by = "taxids")
  df$old_taxids<-df$taxids
  df$taxids<-df$new_taxids
  df$new_taxids=NULL
  df<-cbind(df,do.call(rbind, stringr::str_split(df$path,";")))
  colnames(df)[(length(df)-6):length(df)]<-c("K","P","C","O","F","G","S")
  df$K<-as.character(df$K)
  df$P<-as.character(df$P)
  df$C<-as.character(df$C)
  df$O<-as.character(df$O)
  df$F<-as.character(df$F)
  df$G<-as.character(df$G)
  df$S<-as.character(df$S)
  #not sure why the following command wasnt working
  #df[,(length(df)-6):length(df)] <- sapply(df[,(length(df)-6):length(df)],as.character)
  df[,(length(df)-6):length(df)][df[,(length(df)-6):length(df)]==""]<- "unknown"
  df$path=NULL
  unlink(taxids_fileA)
  unlink(taxids_fileB)
  unlink(taxids_fileC)
  return(df)
}

add.lineage.fasta.BAS<-function(infasta,ncbiTaxDir,out,taxids=T){
  
  if(taxids) {
    message("Reminders: headers must contain \"taxid=taxid;\"")
    fasta.table<-phylotools::read.fasta(infasta)
    fasta.table$taxids<-stringr::str_match(fasta.table$seq.name, "taxid=(.*?);")[,2]
    
  } else {
    message("Reminders: headers must contain \"species=species_name;\"")
    fasta.table<-phylotools::read.fasta(infasta)
    fasta.table$name<-gsub("_"," ",stringr::str_match(fasta.table$seq.name, "species=(.*?);")[,2])
    fasta.table$taxids<-names2taxids(vector = fasta.table$name,ncbiTaxDir)
  }
  
  new_lineage<-add.lineage.df(fasta.table,ncbiTaxDir)
  
  new_lineage$full<-paste0("kingdom=",new_lineage$K,"; phylum=",new_lineage$P,"; class=",new_lineage$C,
                           "; order=",new_lineage$O,"; family=",new_lineage$F,"; genus=",new_lineage$G,
                           "; species=",new_lineage$S,";")
  new_lineage$full<-gsub("=;","=unknown;",new_lineage$full)
  
  #make new headers: "seqid taxid path definition"
  seqid<-gsub(" .*","",new_lineage$seq.name)
  definition<-gsub("^(.*) species=.*;","",new_lineage$seq.name)
  newname<-paste0(seqid," taxid=",new_lineage$taxids,"; ",new_lineage$full,definition)
  #make df for outputting as fasta
  fasta.table.winfo<-as.data.frame(cbind(newname,as.character(new_lineage$seq.text)))
  colnames(fasta.table.winfo)<-c("seq.name","seq.text")
  phylotools::dat2fasta(fasta.table.winfo,outfile = out)
  
}

check.blasts<-function(infastas,h){
  a<-list()
  for(i in 1:length(h)){
    message(infastas[i])
    a[[i]]<-tryCatch(h[[i]]$get_status(),error=function(e) print("not running"))
    if(nchar(a[[i]])!=11) message(gsub("sleeping","running",h[[i]]$get_status()))
    message(paste0("exit status: ",h[[i]]$get_exit_status()))
  }
}

blast.min.bas<-function(infastas,refdb,blast_exec="blastn",wait=T,taxidlimit=NULL,taxidname=NULL,ncbiTaxDir=NULL){
  
  if(!is.null(taxidlimit)) if(is.null(ncbiTaxDir)) stop("to use taxidlimit, ncbiTaxDir must be supplied")
  if(!is.null(taxidlimit)) if(is.null(taxidname)) stop("to use taxidlimit, taxidname must be supplied")
  if(!is.null(taxidlimit)) message("Make sure infastas,taxidlimit & taxidname are in correct order")
  
  t1<-Sys.time()
  
  library(processx)
  
  if(length(infastas)==1 | length(infastas)==2 | length(infastas)==3) threads<-8
  if(length(infastas)==4 | length(infastas)==5 | length(infastas)==6) threads<-4
  if(length(infastas)==7 | length(infastas)==8 | length(infastas)==9) threads<-2
  if(length(infastas)>9) threads<-1
  
  continue<-data.frame("file"<-infastas,"response"="y")
  continue$response<-as.character(continue$response)
  for(i in 1:length(infastas)){
    if(paste0(gsub(x = infastas[i],pattern = ".fasta",replacement = ".blast.txt")) %in% list.files()){
      continue[i,2]<-readline(paste0("The following file already exists, Overwrite? (y/n):", "
                                     ",gsub(x = infastas[i],pattern = ".fasta",replacement = ".blast.txt")))
    }
  }
  
  if("n" %in% continue$response) stop("Abandoned blast due to overwrite conflict")
  
  if(!is.null(taxidlimit)){
    
    h<-list()
    
    for(i in 1:length(infastas)){
      if(length(list.files(pattern = paste0(taxidname[i],"_taxidlimit.txt")))==0){
        
        system2(command = "taxonkit",args = c("list", "--ids", taxidlimit[i], "--indent", '""',"--data-dir",ncbiTaxDir)
                ,wait=T,stdout = paste0(taxidname[i],"_taxidlimit.temp.txt"))
        
        #remove blank row
        taxidlist<-read.table(paste0(taxidname[i],"_taxidlimit.temp.txt"))
        write.table(taxidlist,paste0(taxidname[i],"_taxidlimit.txt"),row.names = F,quote = F,col.names = F)
        
        unlink(paste0(taxidname[i],"_taxidlimit.temp.txt"))
        message(paste("taxidlist saved to",paste0(taxidname[i],"_taxidlimit.txt")))
      }else{ message(paste("The file",paste0(taxidname[i],"_taxidlimit.txt"),"already exists. Using that file."))}
      
      h[[i]]<-process$new(command = blast_exec, 
                          args=c("-query", infastas[i], "-task", "megablast","-db",refdb,"-outfmt",
                                 "6 qseqid evalue staxid pident qcovs","-num_threads", threads, "-taxidlist", 
                                 paste0(taxidname[i],"_taxidlimit.txt"),"-max_target_seqs", "100", "-max_hsps","1", "-out",
                                 paste0(gsub(x = infastas[i],pattern = "\\.fasta",replacement = ".blast.txt"))),echo_cmd = T,
                          stderr = paste0("blast.error.temp.processx.file",i))
    }
  }
  
  if(is.null(taxidlimit)){  
    
    h<-list()
    
    for(i in 1:length(infastas)){
      
      h[[i]]<-process$new(command = blast_exec,
                          args=c("-query", infastas[i], "-task", "megablast","-db",refdb,"-outfmt",
                                 "6 qseqid evalue staxid pident qcovs","-num_threads", threads, "-max_target_seqs", 
                                 "100", "-max_hsps","1", "-out",
                                 paste0(gsub(x = infastas[i],pattern = "\\.fasta",replacement = ".blast.txt"))),
                          echo_cmd = T,stderr = paste0("blast.error.temp.processx.file",i))
    }
  }
  
  Sys.sleep(time = 2)
  
  exits<-list()
  for(i in 1:length(h)){
    exits[[i]]<-h[[i]]$get_exit_status()
  }
  
  if(1 %in% exits){
    message("
            ************
            There was a problem with ", infastas[match(1,exits)], ", aborting all blasts
            ************")
    for(i in 1:length(h)){
      h[[i]]$kill()
    }
  }
  
  if(wait==T){
    for(i in 1:length(h)){
      h[[i]]$wait()
      message(paste(infastas[i],"blast pid",h[[i]]$get_pid()))
      message(readLines(paste0("blast.error.temp.processx.file",i)))
      unlink(paste0("blast.error.temp.processx.file",i))
    }
  }
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  
  message(c("All blasts complete in ",t3," mins."))
  
  return(h)
}

count.MBC<-function(MBCtsvDir,ms_ss,otutab,illumina_script_taxatab,illumina_script_taxatab_tf){
  ##########################################
  message("Currently set up for one run, one primer, modify later if needed")
  if(!"ss_sample_id" %in% colnames(ms_ss)) stop("No column in ms_ss called ss_sample_id")
  #import counts step1
  step1<-data.table::fread(paste0(MBCtsvDir,"step1_stats.tsv"),data.table = F)
  step1$ss_sample_id<-gsub(".none.*$","",step1$Stats)
  demuliplexed_files<-sum(step1[step1$ss_sample_id %in% ms_ss$ss_sample_id,"Number of reads"])/2
  
  #import counts step2
  step2<-data.table::fread(paste0(MBCtsvDir,"step2_stats.tsv"),data.table = F)
  step2$ss_sample_id<-gsub(".none.*$","",step2$Stats)
  after_paired_end<-sum(step2[step2$ss_sample_id %in% ms_ss$ss_sample_id,"Number of reads"])
  
  #import counts step3
  step3<-data.table::fread(paste0(MBCtsvDir,"step3_stats.tsv"),data.table = F)
  step3$ss_sample_id<-gsub(".none.*$","",step3$Stats)
  after_cutadapt<-sum(step3[step3$ss_sample_id %in% ms_ss$ss_sample_id,"nseqs"])
  
  # #import counts step4 - THIS IS UNIQ, SO I GUESS IT DOESNT AFFECT READ COUNT?
  # step4<-data.table::fread(paste0(MBCtsvDir,"step4_stats.tsv"),data.table = F)
  # step4$ss_sample_id<-gsub(".none.*$","",step4$Stats)
  # OTUs_after_cutadapt<-sum(step4[step4$ss_sample_id %in% ms_ss$ss_sample_id,"nseqs"])
  # 
  # sum(step4[,"nseqs"])
  
  #import counts step5
  #command used to make "step5_stats_BAS.tsv"
  
  # find . -name \*none.flash2_merged.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.fasta.gz -print0 
  # | xargs -0 zgrep ">" | sed 's/\.\///;s/.*\///;s/.none.*size=/ /' > step5_stats_BAS.tsv         
  # 
  # step5<-data.table::fread(paste0(MBCtsvDir,"step5_stats_BAS.tsv"),data.table = F)
  # step5_A<-step5[step5$V1 %in% ms_ss$ss_sample_id,]
  # sum(step5_A$V2)
  
  #huh, well after all that, this is the same read count as in otutab, so can skip it
  
  
  #import otu tab
  otutab_A<-data.table::fread(otutab,data.table = F)
  #subset 
  otutab_B<-cbind(OTU=otutab_A$`#OTU ID`,otutab_A[,colnames(otutab_A) %in% ms_ss$ss_sample_id])
  #remove 0-read OTUs and samples
  otutab_C<-rm.0readtaxSam(taxatab = otutab_B)
  #OTUs.in.otutab<-length(otutab_C$OTU) #NOT REALLY NECESSARY
  after_size_select<-sum(otutab_C[,-1])
  
  #import first taxa tab
  taxatab_A<-data.table::fread(illumina_script_taxatab,data.table = F)
  #subset
  taxatab_B<-cbind(taxon=taxatab_A$taxon,taxatab_A[,colnames(taxatab_A) %in% ms_ss$ss_sample_id])
  #remove 0-read OTUs and samples
  taxatab_C<-rm.0readtaxSam(taxatab = taxatab_B)
  ##this is the same as otutab
  
  after_blast<-sum(taxatab_C[taxatab_C$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",-1])
  after_blast_filt<-after_blast-sum(taxatab_C[taxatab_C$taxon=="NA;NA;NA;NA;NA;NA;NA",-1])
  
  #import filtered taxa tab
  taxatab_A<-data.table::fread(illumina_script_taxatab_tf,data.table = F)
  #subset 
  taxatab_B<-cbind(taxon=taxatab_A$taxon,taxatab_A[,colnames(taxatab_A) %in% ms_ss$ss_sample_id])
  #remove 0-read OTUs and samples
  taxatab_C<-rm.0readtaxSam(taxatab = taxatab_B)
  #remove no hits and NA
  taxatab_D<-taxatab_C[taxatab_C$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]
  taxatab_D<-taxatab_D[taxatab_D$taxon!="NA;NA;NA;NA;NA;NA;NA",]
  after.taxon.filter<-sum(taxatab_D[,-1])
  
  
  out<-data.frame("Demuliplexed files"=demuliplexed_files,
                  "After paired end alignment"=after_paired_end,
                  "After primer trimming"=after_cutadapt,
                  "After size selection and singleton removal"=after_size_select,
                  "After blast"=after_blast,
                  "After blast filters"=after_blast_filt,
                  "After initial taxon filter"=after.taxon.filter)
  
  colnames(out)<-gsub("\\."," ",colnames(out))
  
  return(out)
}

google.read.master.url<-function(sheeturl,out=NULL,ws="Master_Samplesheet"){
  url2<-stringr::str_split(sheeturl,"/d/")[[1]][2]
  url2<-stringr::str_split(url2,"/")[[1]][1]
  ss_info<-googlesheets4::sheets_get(ss = url2)
  if(ws == "ENA_sample_data") {
    ss_data<-googlesheets4::read_sheet(ss = url2,sheet = ws,col_types = "c") 
    colnames(ss_data)<-ss_data[2,]
    ss_data<-ss_data[3:length(ss_data$sample_alias),]
    ss_data<-as.data.frame(ss_data[!is.na(ss_data$sample_alias),])
  } else{
    if(ws == "ENA_library_data") { 
      ss_data<-googlesheets4::read_sheet(ss = url2,sheet = ws,col_types = "c") 
      ss_data<-as.data.frame(ss_data[!is.na(ss_data$sample_alias),])
    } else{
      if(ws == "Library_data") { 
        ss_data<-googlesheets4::read_sheet(ss = url2,sheet = ws,col_types = "c") 
        ss_data<-as.data.frame(ss_data[!is.na(ss_data$run_alias),])
      } else{
        ss_data<-googlesheets4::read_sheet(ss = url2,sheet = ws,col_types = "c") 
        ss_data<-as.data.frame(ss_data[!is.na(ss_data$Sample_Name),])
      }
    }
  }
  
  
  if(!is.null(out)){
    write.table(ss_data,file = paste0(gsub(" ","_",ss_info$name),"_",gsub(" ","_",ws),".txt"),quote = F,row.names = F,sep = "\t")
    message(paste("file saved as",paste0(gsub(" ","_",ss_info$name),"_",gsub(" ","_",ws),".txt")))
  }
  return(ss_data)
}

#' Interpret glm(er) with a binary response
#' @title Interpret glm(er) with a binary response
#' @param model A model produced by \code{\link[stats]{glm}} or \code{\link[lmer4]{glmer}}, where the response variable is binary.
#' @return Prints one sentence for each predictor variable in the form: "for each increase of 1 unit \code{predictor},
#'     \code{response} increases \code{x} times (p=\code{p value})"
#' @note The response variable must be binary to be accurate, as the function reports the exp of the log odds.
#'
#' @examples
#' interp.glm.bin(model)
#'
#' @export
interp.glm.bin<- function(model){
  for (i in 2:length(summary(model)$coefficients[,1])){
    print(paste0("for each increase of 1 unit ", rownames(summary(model)$coefficients)[i], ", ",
                 strsplit(as.character(summary(model)$call[2])," ~")[[1]][1]," increases ",
                 round(exp(summary(model)$coefficients[i,"Estimate"]),digits = 4),
                 " times (p=",round(summary(model)$coefficients[i,"Pr(>|z|)"],digits=4),")"))
  }}

#' Interpret lm(er)
#' @title Interpret lm(er)
#' @param model A model produced by \code{\link[stats]{lm}} or \code{\link[lmer4]{lmer}}.
#' @return Prints one sentence for each predictor variable in the form: "for each increase of 1 unit (predictor),
#'     response increases \code{x} units (p=(p value))"
#'
#' @examples
#' interp.lm(model)
#'
#' @export
interp.lm <- function(model){
  for (i in 2:length(summary(model)$coefficients[,1])){
    print(paste0("for each increase of 1 unit ", rownames(summary(model)$coefficients)[i],", ",
                 strsplit(as.character(summary(model)$call[2])," ~")[[1]][1]," increases ",
                 round(summary(model)$coefficients[i,"Estimate"],digits = 4),
                 " units (p=",round(summary(model)$coefficients[i,"Pr(>|t|)"],digits=4),")"))
  }}


make.blastdb.bas<-function(infasta,makeblastdb_exec="makeblastdb",addtaxidsfasta=T, ncbiTaxDir, dbversion=5){
  library(processx)
  
  if(!infasta %in% list.files()) stop("infasta not found in current directory")
  
  message("Reminder: Only works for blast 2.9.0: Current version:")
  system2(command = makeblastdb_exec,args = c("-version"))
  
  tempfileA<-paste0("temp",as.numeric(Sys.time()),".fasta")
  
  if(addtaxidsfasta==T){  
    message("assumes fasta with species=xxx;") 
    add.lineage.fasta.BAS(infasta,ncbiTaxDir,out = tempfileA)
    tempfasta<-phylotools::read.fasta(tempfileA)
    
    #remove seqs with no taxonomy 
    if(length(tempfasta$seq.name[grep("kingdom=unknown",tempfasta$seq.name)])>0){
      message(paste0("Discarding ",length(tempfasta$seq.name[grep("kingdom=unknown",tempfasta$seq.name)]),
                     " sequences for which taxonomy could not be found"))
      tempfasta<-tempfasta[-grep("kingdom=unknown",tempfasta$seq.name),]
    }
    phylotools::dat2fasta(tempfasta,tempfileA)
    
  } else{message("Assumes fasta with taxid=xxx;")
    file.copy(infasta,tempfileA)}
  
  #give unique ids
  message("giving unique ids")
  tempfileB<-paste0("temp",as.numeric(Sys.time()),".fasta")
  system2(command = "obiannotate", args=c("--uniq-id",tempfileA), stdout=tempfileB,stderr = "",wait = T)
  
  #Add "database name to header
  message("adding db name to headers")
  sedarg=paste0("s/taxid=/database=",gsub(".fasta","",infasta),"; taxid=/g")
  h<-process$new(command = "sed", args=c(sedarg, tempfileB), echo_cmd = T,
                 stdout=gsub(".fasta",".blastdbformatted.fasta",infasta))
  h$wait()
  
  #ensure ids are <50 characters
  tempfasta<-phylotools::read.fasta(gsub(".fasta",".blastdbformatted.fasta",infasta))
  defs<-suppressWarnings(do.call(rbind,strsplit(as.character(tempfasta$seq.name)," ")))
  defs[,1]<-stringr::str_trunc(as.character(defs[,1]),width = 49)
  #remove any quotes
  defs[,1]<-gsub('"',"",defs[,1])
  tempfasta$seq.name <- apply(as.data.frame(defs),1,paste,collapse = " " )
  phylotools::dat2fasta(tempfasta,gsub(".fasta",".blastdbformatted.fasta",infasta))
  
  #make mapping file
  taxids<-sub(x = stringr::str_match(string = tempfasta$seq.name,pattern = "taxid=(.*)"),pattern = ";.*",
              replacement = "")[,2]
  
  taxonmap<-as.data.frame(cbind(defs[,1],taxids))
  mappingfile<-paste0("mapping",as.numeric(Sys.time()),".txt")
  write.table(taxonmap,mappingfile,quote = F,sep = " ",row.names = F,col.names = F)
  
  system2(command = makeblastdb_exec, 
          args=c("-in", gsub(".fasta",".blastdbformatted.fasta",infasta), 
                 "-dbtype", "nucl", 
                 "-blastdb_version", dbversion,
                 "-parse_seqids","-taxid_map",mappingfile,"-out",
                 gsub(".fasta","",infasta)), stderr = "",wait = T)
  
  #testblast
  message("Running test blast")
  phylotools::dat2fasta(head(tempfasta,n=1),gsub(".fasta",".blastdbformatted.test.fasta",infasta))
  h<-blast.min.bas(gsub(".fasta",".blastdbformatted.test.fasta",infasta),refdb = gsub(".fasta","",infasta))
  check.blasts(gsub(".fasta",".blastdbformatted.test.fasta",infasta),h)
  message("Does new db have taxids?")
  print(data.table::fread(gsub(".fasta",".blastdbformatted.test.blast.txt",infasta)))
  system2(command = "blastdbcheck",args = c("-must_have_taxids","-db",gsub(".fasta","",infasta)))
  
  unlink(tempfileA)
  unlink(tempfileB)
  unlink(gsub(".fasta",".blastdbformatted.test.blast.txt",infasta))
  unlink(gsub(".fasta",".blastdbformatted.test.fasta",infasta))
}

#' Convert sequence files between various formats
#' @param infile A file of type genbank, embl (not tested), fasta, fastq, ecopcrdb.
#' @param in_type Accepted values: "genbank", "embl", "fasta", "fastq", "ecopcrdb".
#' @param taxo An obitools-formatted taxonomy database as an R object, generated by \code{\link[ROBITaxonomy]{read.taxonomy}}
#'     or \code{\link[bastools]{NCBI2obitaxonomy}}. This is only used when creating an ecopcrdb - taxids must be in header in format
#'     > seqX taxids=8276
#'     see \code{obiaddtaxids.Bas}.
#' @param skip_errors logical. skips entries where taxonomy could not be found. Only used when creating an ecopcrdb.
#' ########################Hangs without this...make mandatory. also taxids are mandatory
#' @param out_type Accepted values: "--fasta-output", "--ecopcrdb-output"
#' @note Creating an ecopcrdb takes much longer here than running in a terminal, not sure why! (e.g. 20 seqs took).
#'     Looks like an issue creating the .adx file of new taxonomy.
#' @return A file of type fasta or ecopcrdb
#' @examples
#' fastq->fasta
#' a<-system.file("extdata", "head.test.minion.fastq", package = "bastools")
#' obitaxonomydb<-"/media/sf_Documents/WORK/CIBIO/STATS_AND_CODE/TAXONOMIES/obitax_26-4-19"
#' obiconvert.Bas2(a,taxo = obitaxonomydb, out_type = "--fasta-output",out="head.test.minion.fasta",in_type="fastq")
#' ok
#'
#' fastq->ecopcrdb, error
#' ok
#'
#' fasta without taxids->ecopcrdb
#' a<-system.file("extdata", "head.test.fasta", package = "bastools")
#' obitaxonomydb<-"/media/sf_Documents/WORK/CIBIO/STATS_AND_CODE/TAXONOMIES/obitax_26-4-19"
#' obiconvert.Bas(a,taxo = obitaxonomydb, out_type = "--ecopcrdb-output",out="head.test.ecopcrdb",in_type="fasta")
#'#############hangs
#'
#' fasta with taxids->ecopcrdb
#' a<-system.file("extdata", "wTaxids.head20.test.fasta", package = "bastools")
#' obitaxonomydb<-"/media/sf_Documents/WORK/CIBIO/STATS_AND_CODE/TAXONOMIES/obitax_26-4-19"
#' obiconvert.Bas(a,taxo = obitaxonomydb, out_type = "--ecopcrdb-output",out="wTaxids.head.test.ecopcrdb",in_type="fasta")
#' #################NOT WORKING PROPERLY, SUPER SLOW.....
#' #14:27 start (20 seqs with taxids)
#'
#' genbank->fasta
#' a<-system.file("extdata", "head_ALL_VERTS_REFSEQ_MTDNA.gb", package = "bastools")
#' obitaxonomydb<-"/media/sf_Documents/WORK/CIBIO/STATS_AND_CODE/TAXONOMIES/obitax_26-4-19"
#' obiconvert.Bas(a,taxo = obitaxonomydb, out_type = "--ecopcrdb-output",out="head.test.ecopcrdb")
#'
#' genbank to ecopcrdb
#'
#' @export
obiconvert.Bas<-function(infile,in_type,taxo,out_type,out){
  if(in_type!="fasta" & out_type=="--ecopcrdb-output") stop("Only fasta files can be converted into an ecopcrdb")
  
  cb <- function(line, proc) {cat(line, "\n")}
  if(out_type=="--fasta-output"){
    b<-processx::run(command = "obiconvert", args=c(infile,"-d",taxo,out_type),
                     echo=F,stderr_line_callback = cb,echo_cmd = T)
    writeLines(b$stdout,con = out)
  }
  if(out_type=="--ecopcrdb-output"){
    message("Make sure input fasta file has taxids in sequence headers in the format \" taxid=XXXX;\"")
    processx::run(command = "obiconvert", args=c(infile,"-d",taxo,paste0(out_type,"=",out),"--skip-on-error"),
                  echo=F,stderr_line_callback = cb,echo_cmd = T)
    d<-"SUCCESS!"
  }
  h<-processx::run(command = "obicount", args=c(infile),echo=F,echo_cmd = F,stderr_line_callback = cb)
  j<-stringr::str_locate_all(string = h$stdout,pattern = " " )[[1]][,1]
  k<-substr(h$stdout,1,j-1)
  l<-paste0("input_read_count=",k)
  h<-processx::run(command = "obicount", args=c(out),echo=F,echo_cmd = F,stderr_line_callback = cb)
  j<-stringr::str_locate_all(string = h$stdout,pattern = " " )[[1]][,1]
  k<-substr(h$stdout,1,j-1)
  m<-paste0("output_read_count=",k)
  print("SUCCESS!")
  print(l)
  print(m)
}

#' @export
bascount.fasta<-function(infile){
  f<-processx::run(command = "grep",args = c("-c",'>',infile),echo_cmd = F,echo = T)
  as.numeric(f$stdout)
}

#' @export
bascount.fastas.recur<-function(folder){
  cb <- function(line, proc) {cat(line, "\n")}
  origpath<-getwd()
  setwd(folder)
  a<-list.files(pattern = "*.fasta",recursive = T)
  b<-processx::run(command = "grep",args = c("-cH",'>',a),echo_cmd = F,echo = F,stderr_line_callback = cb)
  d<-read.table(text = b$stdout,sep = "\n")
  print(d)
  e<-stringr::str_split(string = d$V1,pattern = ":")
  setwd(origpath)
  sum(sapply(e, function(x){as.numeric(x[2])}))
}

bascount.gb.recur<-function(folder){
  cb <- function(line, proc) {cat(line, "\n")}
  origpath<-getwd()
  setwd(folder)
  a<-list.files(pattern = "*.gb",recursive = T)
  b<-processx::run(command = "grep",args = c("-cH",'//',a),echo_cmd = F,echo = F,stderr_line_callback = cb)
  d<-read.table(text = b$stdout,sep = "\n")
  print(d)
  e<-stringr::str_split(string = d$V1,pattern = ":")
  setwd(origpath)
  sum(sapply(e, function(x){as.numeric(x[2])}))
}

bascount.gb<-function(gbfile){
  b<-processx::run(command = "grep",args = c('-c','//',gbfile),echo_cmd = F,echo = F)
  d<-read.table(text = b$stdout,sep = "\n")
  return(d$V1)
}

remove.dashes.fasta<-function(infasta,outfasta){
  f<-process$new(command = "sed",args = c("-e","/^>/!s/-//g;/^$/d",infasta),echo_cmd = T, stdout = outfasta)
  f$wait()
}

#' Merge MBC OTU table and bin_blast_results
#'
#' @param MBC_otutab MBC_otutab
#' @param bin_blast_results bin_blast_results
#' @return Merged tab-separated taxa table 
#' @export
MBC_otu_bin_blast_merge<-function(MBC_otutab, bin_blast_results,out){
  otutab<-data.table::fread(file=MBC_otutab,sep = "\t",header = T)
  bins<-data.table::fread(file=bin_blast_results,sep = "\t",header = T)
  
  colnames(otutab)[1]<-"OTU_name"
  
  bins$OTU_name<-do.call(rbind,stringr::str_split(bins$qseqid,";"))[,1]
  
  otutab.bins<-as.data.frame(merge(x = otutab,y = bins[,2:9], by = "OTU_name"))
  otutab.bins.all<-as.data.frame(merge(x = otutab,y = bins[,2:9], by = "OTU_name",all.x=T))
  no.hit.otutab<-otutab.bins.all[!otutab.bins.all$OTU_name %in% otutab.bins$OTU_name,]
  if(length(no.hit.otutab$OTU_name)!=0){
    no.hit.otutab[,c((length(colnames(otutab.bins))-6):length(colnames(otutab.bins)))]<-"no_hits"}
  otutab.bins<-rbind(otutab.bins,no.hit.otutab)
  
  #Reorder columns
  otutab.bins<-otutab.bins[,c(1,(length(colnames(otutab.bins))-6):length(colnames(otutab.bins)),
                              2:(length(colnames(otutab.bins))-7))]
  #sep paths from read counts
  path<-paste0(otutab.bins$K,";",otutab.bins$P,";",otutab.bins$C,";", otutab.bins$O,";",otutab.bins$F,";",
               otutab.bins$G,";",otutab.bins$S)
  otutab.bins<-otutab.bins[,c(9:length(colnames(otutab.bins)))]
  taxa.table<-aggregate(x = otutab.bins,by = list(path),FUN = sum)
  colnames(taxa.table)<-gsub("Group.1","taxon",colnames(taxa.table))
  write.table(x = taxa.table,file = out,sep="\t",quote = F,row.names = F)
}

#' Combine an OTU table with a table of taxon names. Most commonly where \code{obitab} was used to convert a fasta file to a table,
#'    and \code{MEGAN} was used to assign taxonomy from BLAST results of the same fasta file.
#' @title Merge reads and assigned taxonomy
#' @param obitab.txt Any tab-delineated text file with a column "id" containing read names,
#'     such as the file created by \code{obitab}.
#' @param megan.taxa Can be:
#'    \itemize{
#'     \item A simple, headerless text file where the first column consist of the read names
#'     and the second column consists of taxa, as manually output using the \code{MEGAN::readName_to_taxonName} option.
#'     \item A dataframe consisting of a taxonomy table with read names in a column named "id", as ouput by \code{rma2info.BAS}}
#' @return A dataframe which is equal to \code{obitab.txt} but, depending on input, either has one new column "taxon",
#'     or multiple columns corresponding to taxon path (SK,P,C,O,F,G,S) and lowest taxonomic level reached.
#'
#' @examples
#'   \itemize{
#'     \item test<-merge.tab.taxon(c20.uniq.l85L105.PRIMER_16S.tab", "Mblast.c20.uniq.l85L105.PRIMER_16S-taxon.txt")
#'     \item test<-merge.tab.taxon(c20.uniq.l85L105.PRIMER_16S.tab", taxon.table)
#'     \item blast2rma.BAS("primer_16S.uniq.l75L120.c20.xml",outfile = "primer_16S.uniq.l75L120.c20.TEST.blast2rma.rma6",
#'               a2t = "nucl_acc2tax-Nov2018.abin")
#'              inspect file and disable taxa as necessary
#'              primer_16S.uniq.l75L120.c20.TEST.taxon.table<-rma2info.BAS("primer_16S.uniq.l75L120.c20.TEST.blast2rma.rma6")
#'              final.table<-merge.tab.taxon(obitab.txt = "primer_16S.uniq.l75L120.c20.tab",primer_16S.uniq.l75L120.c20.TEST.taxon.table)}
#' @export
otutab_bin_blast_merge_minion<-function(otutabfile,binfile,experimentsheetfile=NULL,experiment_id,hascount,out){
  
  otutab_input<-data.table::fread(otutabfile, sep = "\t",header = T)
  
  if(hascount) if(!"count" %in% colnames(otutab_input)) stop("no column named count in otutab")
  if(hascount==F) otutab_input$count<-1
  
  taxon_input<-data.table::fread(file = binfile, sep = "\t")
  
  taxon_input$path<-paste0(taxon_input$K,";",taxon_input$P,";",taxon_input$C,";",taxon_input$O,";",
                           taxon_input$F,";",taxon_input$G,";",taxon_input$S)
  if("barcode" %in% colnames(otutab_input)){
    merged.table<-merge(taxon_input[,c("qseqid","path")],otutab_input[,c("id","barcode","count")],
                        by.x = "qseqid",by.y = "id",all = TRUE)
    
    taxatable<-reshape2::dcast(merged.table[,c("path","barcode","count")],path~barcode,value.var = "count",
                               fun.aggregate = sum)
    
    if(!is.null(experimentsheetfile)){
      #read sheet
      experimentsheet<-as.data.frame(read.table(paste0(outDir,experiment_id,"_experiment_sheet.txt"),header = T))
      #replace barcodes with sample_names
      final.barcodes<-as.data.frame(colnames(taxatable[,2:length(colnames(taxatable))]))
      colnames(final.barcodes)<-"barcode_id"
      mapping.samples<-experimentsheet[grep(experiment_id, experimentsheet$experiment_id),c("barcode_id","ss_sample_id")]
      mapping.samples$barcode_id<-gsub("BC","barcode",mapping.samples$barcode_id)
      final.samples<-merge(final.barcodes,mapping.samples,by = "barcode_id",all.y = F,all.x = T)
      final.samples$barcode_id<-as.character(final.samples$barcode_id)
      final.samples$ss_sample_id<-as.character(final.samples$ss_sample_id)
      #final.samples$ss_sample_id2<-do.call(pmax, c(final.samples, na.rm = TRUE))
      colnames(taxatable)<-c("taxon",as.character(final.samples$ss_sample_id))
    }
    
  } else {
    merged.table<-merge(taxon_input[,c("qseqid","path")],otutab_input[,c("id","ss_sample_id","count")],
                        by.x = "qseqid",by.y = "id",all = TRUE)
    
    taxatable<-reshape2::dcast(merged.table[,c("path","ss_sample_id","count")],path~ss_sample_id,value.var = "count",
                               fun.aggregate = sum)
    colnames(taxatable)<-gsub("path","taxon",colnames(taxatable))
  }
  
  
  taxatable$taxon[is.na(taxatable$taxon)]<-"No_hits"
  
  write.table(x = taxatable,file = out,sep="\t",quote = F,row.names = F)
}


######################################################################
#SPLIT TAXATABLES
splice.taxatables<-function(files,mastersheet){
  
  message("Note:sample names must contain project name with dash")
  
  #read files
  taxatables<-list()
  for(i in 1:length(files)){
    taxatables[[i]]<-data.table::fread(files[i],data.table = F,sep = "\t")
  }
  
  #split by project
  ssdf<-data.table::fread(file = mastersheet,sep = "\t")
  projectnames<-suppressWarnings(unique(do.call(rbind,stringr::str_split(string = ssdf$ss_sample_id,pattern = "-"))[,1]))
  projectnames<-paste0(projectnames,"-")
  
  taxatablesplit<-list()
  taxatablesplit2<-list()
  
  for(i in 1:length(taxatables)){
    for(j in 1:length(projectnames)){
      taxatablesplit[[j]]<-as.data.frame(taxatables[[i]][,grep(projectnames[j],colnames(taxatables[[i]]))])
      colnames(taxatablesplit[[j]])<-grep(projectnames[j],colnames(taxatables[[i]]),value = T) #to ovecome dfs with one sample only
      rownames(taxatablesplit[[j]])<-taxatables[[i]]$taxon
    }
    taxatablesplit2[[i]]<-taxatablesplit
    names(taxatablesplit2[[i]])<-projectnames
  }
  
  #write
  for(i in 1:length(taxatablesplit2)){
    for(j in 1:length(projectnames)){
      if(length(taxatablesplit2[[i]][[j]])>1) {
        #remove 0 read taxa
        taxatablesplit2[[i]][[j]]<-taxatablesplit2[[i]][[j]][rowSums(taxatablesplit2[[i]][[j]])!=0,]
        taxatablesplit2[[i]][[j]]$taxon<-rownames(taxatablesplit2[[i]][[j]])
        taxatablesplit2[[i]][[j]]<-taxatablesplit2[[i]][[j]][,c(length(colnames(taxatablesplit2[[i]][[j]])),1:(length(colnames(taxatablesplit2[[i]][[j]]))-1))]
        if(length(taxatablesplit2[[i]][[j]]$taxon)>0){
          write.table(taxatablesplit2[[i]][[j]],
                      paste0(projectnames[j],gsub("taxatable.tf.txt","taxatable.tf.spliced.txt",files[i])),
                      row.names = F,quote = F,sep = "\t")
        }
      }
    }
  }
}

#taxatable to krona format
bas.krona.plot<-function(taxatable,KronaPath=NULL){
  
  a<-read.table(taxatable,header = T,sep = "\t")
  
  b<-as.data.frame(do.call(rbind, stringr::str_split(a[,1],";")))
  colnames(b)<-c("K","P","C","O","F","G","S")
  
  a$all<-rowSums(a[,2:length(colnames(a))])
  d<-colnames(a[,2:length(colnames(a))])
  
  for(i in 1:length(d)){
    sample<-cbind(a[,d[i]],b)
    colnames(sample)[1]<-d[i]
    write.table(sample,row.names = F,file = paste0(d[i],".krona.txt"),quote = F,sep = "\t",col.names = F)
  }
  
  if(!is.null(KronaPath)){
    command<-KronaPath} else {command<- "ktImportText"}
  
  system2(command = command,args = c(list.files(pattern = "*krona.txt"),"-o", paste0(gsub(".txt",".krona.html",taxatable)))
          ,stdout = F,stderr = "",wait = T)
  
  unlink(list.files(pattern = "*krona.txt"))
}

taxatab.stackplot<-function(taxatab,master_sheet=NULL,column=NULL,as.percent=T,as.dxns=F,facetcol=NULL,hidelegend=F){
  taxa<-do.call(rbind,str_split(taxatab$taxon,";"))
  taxa<-cbind(taxa,do.call(rbind,str_split(taxa[,7]," ")))
  taxa2<-as.data.frame(substr(taxa,start = 1,stop = 3))
  taxatab$taxon<-apply(taxa2,MARGIN = 1,FUN = function(x) paste0(x[1],".",x[2],".",x[3],".",x[4],".",x[5],".",x[6],".",x[8],"_",x[9]))
  
  if(as.dxns==T) taxatab<-binarise.taxatab(taxatab,t=T)
  
  long<-reshape2::melt(taxatab)
  long<-long[long$value>0,]
  
  if(!is.null(master_sheet)) {
    long<-merge(x = long,y = master_sheet[,c(grep("ss_sample_id",colnames(master_sheet)),
                                             grep(column,colnames(master_sheet)))],
                by.x="variable",by.y="ss_sample_id")
    
    if(!is.null(facetcol)){
      long<-merge(x = long,y = master_sheet[,c(grep("ss_sample_id",colnames(master_sheet)),
                                               grep(facetcol,colnames(master_sheet)))],
                  by.x="variable",by.y="ss_sample_id")
    }
    
    long$variable<-long[,grep(column,colnames(long))]
    
  } else message("No master_sheet or column specified, making default plot")
  
  #long<-long[order(long$taxon),]
  
  if(as.dxns) long<-long[!duplicated(long),]
  
  a<-ggplot2::ggplot(data=long , aes(y=value, x=variable, fill=taxon))+
    theme(legend.title = element_text(size=10), legend.text=element_text(size=10),
          axis.text.x=element_text(size=8,angle=45, hjust=1),legend.position="right",legend.direction="vertical")
  
  if(length(unique(long$taxon))<29)  a<-a+scale_fill_manual(values = MyCols) 
  
  if(as.percent) {a<-a+geom_bar(position="fill", stat="identity")} else {a<-a+geom_bar(stat = "identity")}
  
  if(!is.null(facetcol))  a<-a+facet_wrap(as.formula(paste0("~",facetcol)))
  
  if(hidelegend) {
    message("Outputting as a list where first element is plot and second is legend")
    plotlist<-list()
    plotlist[[1]]<-a+theme(legend.position="none")
    plotlist[[2]]<-ggpubr::as_ggplot(cowplot::get_legend(a))
    return(plotlist)
  }else return(a)
}


#Plotting bray distance matrix PCA
taxatab.pca.plot<-function(taxatab,master_sheet,MS_colID="ss_sample_id",factor1,lines=F,longnames=F,shortnames=F,ellipse=T){
  taxatab<-rm.0readtaxSam(taxatab)
  taxatab2<-binarise.taxatab(taxatab)
  distance_matrix<-taxatab2bray(taxatab2)
  
  cmds<-cmdscale(distance_matrix,k=4, list. = T, eig = T)
  cor.cmds<-cor(taxatab2,cmds$points)
  VarExplainedPC1<-round(cor(vegan::vegdist(cmds$points[,1],method = "euclidean"),distance_matrix)^2,digits = 2)
  VarExplainedPC2<-round(cor(vegan::vegdist(cmds$points[,2],method = "euclidean"),distance_matrix)^2,digits = 2)
  
  #create loadings for plotting
  loadings<-as.data.frame(cor.cmds)
  cmdspoints<-as.data.frame(cmds$points)
  cmdspoints$ss_sample_id<-rownames(cmdspoints)
  cmdspoints<-merge(cmdspoints,master_sheet,by=MS_colID,all.x = T)
  
  p<-ggplot(cmdspoints,aes(x=V1,y=V2))+
    geom_point(aes(size=1,colour=cmdspoints[,factor1]))+
    #geom_jitter(position = )
    xlab(bquote("Variance explained =" ~ .(VarExplainedPC1)))+
    ylab(bquote("Variance explained =" ~ .(VarExplainedPC2))) +
    theme_bw()+
    guides(size = FALSE)+
    theme(legend.title = element_blank())+
    theme(legend.text=element_text(size=12))+
    theme(legend.spacing.x = unit(0.2, 'cm'))+
    theme(axis.title = element_text(size = 12))
  
  message("Principal Component Analysis plot of community similarity using Bray-Curtis distances")
  
  if(lines){
    p<- p +geom_segment(data = loadings, aes(x=0,y=0,xend=V1,yend=V2),arrow=arrow(length=unit(0.1,"cm")))
    if(longnames) if(shortnames) stop("Can only use EITHER long OR short names")
    if(!longnames) if(!shortnames) message("No names added")
    if(longnames) if(!shortnames) p<- p + geom_text(data = loadings, aes(x=V1, y=V2, label=colnames(taxatab2)))
    if(shortnames) if(!longnames){
      taxa<-do.call(rbind,str_split(colnames(taxatab2),";"))
      taxa<-cbind(taxa,do.call(rbind,str_split(taxa[,7]," ")))
      taxa<-as.data.frame(substr(taxa,start = 1,stop = 3))
      taxa<-apply(taxa,MARGIN = 1,FUN = function(x) paste0(x[1],".",x[2],".",x[3],".",x[4],".",x[5],".",x[6],".",x[8],".",x[9]))
      p<- p + geom_text(data = loadings, aes(x=V1, y=V2, label=taxa))
    }
  }
  
  if(ellipse){
    p<- p +stat_ellipse(aes(linetype=cmdspoints[,factor1]),type = "norm", level=0.90) 
    message("ellipses are drawn with a confidence level of 0.90")
  }
  
  p
}

#Plotting bray distance matrix PCA
taxatab.pca.plot.col<-function(taxatab,master_sheet,MS_colID="ss_sample_id",factor1,lines=F,longnames=F,shortnames=F,ellipse=T){
  taxatab<-rm.0readtaxSam(taxatab)
  taxatab2<-binarise.taxatab(taxatab)
  distance_matrix<-taxatab2bray(taxatab2)
  
  cmds<-cmdscale(distance_matrix,k=4, list. = T, eig = T)
  cor.cmds<-cor(taxatab2,cmds$points)
  VarExplainedPC1<-round(cor(vegan::vegdist(cmds$points[,1],method = "euclidean"),distance_matrix)^2,digits = 2)
  VarExplainedPC2<-round(cor(vegan::vegdist(cmds$points[,2],method = "euclidean"),distance_matrix)^2,digits = 2)
  
  #create loadings for plotting
  loadings<-as.data.frame(cor.cmds)
  cmdspoints<-as.data.frame(cmds$points)
  cmdspoints[,MS_colID]<-rownames(cmdspoints)
  cmdspoints<-merge(cmdspoints,master_sheet,by=MS_colID,all.x = T)
  
  #plot
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  p<-ggplot(cmdspoints,aes(x=V1,y=V2))+
    geom_point(aes(size=1,colour=cmdspoints[,factor1],shape=cmdspoints[,factor1],stroke=1))+
    #geom_jitter(position = )
    scale_shape_manual(values=c(1, 2, 0, 5, 6, 3, 4))+
    scale_color_manual(values = c(cbPalette))+
    xlab(bquote("Variance explained =" ~ .(VarExplainedPC1)))+
    ylab(bquote("Variance explained =" ~ .(VarExplainedPC2))) +
    theme_bw()+
    #labs(shape = factor1)+
    guides(size = FALSE)+
    guides(shape=guide_legend(override.aes = list(size = 4)))+
    theme(legend.title = element_blank())+
    theme(legend.text=element_text(size=12))+
    theme(legend.spacing.x = unit(0.2, 'cm'))+
    theme(axis.title = element_text(size = 12))
  
  message("Principal Component Analysis plot of community simmilarity using Bray-Curtis distances")
  
  if(ellipse){
    p<- p+stat_ellipse(aes(colour=cmdspoints[,factor1] 
                           ,fill=cmdspoints[,factor1]
    )
    ,type = "norm", level=0.90, 
    geom = "polygon",alpha=0.2,
    show.legend = F,segments = 100) +
      scale_fill_manual(values = c(cbPalette))
    
    p$layers<-rev(p$layers)
    
    message("ellipses are drawn with a confidence level of 0.90")
  }
  
  if(lines){
    p<- p +geom_segment(data = loadings, aes(x=0,y=0,xend=V1,yend=V2),arrow=arrow(length=unit(0.1,"cm")))
    if(longnames) if(shortnames) stop("Can only use EITHER long OR short names")
    if(!longnames) if(!shortnames) message("No names added")
    if(longnames) if(!shortnames) p<- p + geom_text(data = loadings, aes(x=V1, y=V2, label=colnames(taxatab2)))
    if(shortnames) if(!longnames){
      taxa<-do.call(rbind,str_split(colnames(taxatab2),";"))
      taxa<-cbind(taxa,do.call(rbind,str_split(taxa[,7]," ")))
      taxa<-as.data.frame(substr(taxa,start = 1,stop = 3))
      taxa<-apply(taxa,MARGIN = 1,FUN = function(x) paste0(x[1],".",x[2],".",x[3],".",x[4],".",x[5],".",x[6],".",x[8],".",x[9]))
      p<- p + geom_text(data = loadings, aes(x=V1, y=V2, label=taxa))
    }
  }
  p
}

read.ecopcrdb<-function(ecopcr.db,ecopcr.tab,obitaxoR){
  message("Converting to tab")
  f<-process$new(command = "obitab", args=c("-o",ecopcr.db), echo_cmd = T,stdout=ecopcr.tab)
  f$wait()
  message("Reading in ecopcrdb")
  ECOPCRDB<-as.data.frame(data.table::fread(file = ecopcr.tab,header = TRUE,sep = "\t"))
  message("Adding taxonomy path")
  taxon.table=get.classic.taxonomy.Bas(ECOPCRDB,obitaxdb = obitaxoR)
  taxon.table$class_name_ok<-as.character(taxon.table$class_name_ok)
  taxon.table$order_name_ok<-as.character(taxon.table$order_name_ok)
  ECOPCRDB$kingdom_name<-taxon.table$kingdom_name_ok
  ECOPCRDB$phylum_name<-taxon.table$phylum_name_ok
  ECOPCRDB$class_name<-taxon.table$class_name_ok
  ECOPCRDB$order_name<-taxon.table$order_name_ok
  ECOPCRDB$family_name<-taxon.table$family_name_ok
  ECOPCRDB$genus_name<-taxon.table$genus_name_ok
  ECOPCRDB$species_name<-taxon.table$species_name_ok
  ECOPCRDB
}


#' Find primers to amplify a group of taxa
#' @title Find primers to amplify a group of taxa.
#' @param ecopcrdb An ecopcrdb containing the sequences to be tested
#' @param max_error Max number of mismatches allowed
#' @param min_length Minimum amplicon length (without primers - I think)
#' @param max_length Maximum amplicon length (without primers - I think)
#' @param strict_match_p Proportion of sequences with identical primer match
#' @param e_match_p Proportion of sequences with primer match within \code{e}
#' @return Default: A list containing the following dataframes
#' \itemize{
#'     \item \code{primer_table}: A table with the following columns
#'         \itemize{
#'             \item \code{PID} Primer ID
#'             \item \code{Pf} Forward primer sequence
#'             \item \code{Pr} Reverse primer seqeunce
#'             \item \code{TmPf} Melting temperature forward primer (without mismatches)
#'             \item \code{MinTmPf} Minimum melting temperature forward primer
#'             \item \code{TmPr} Melting temperature reverse primer (without mismatches)
#'             \item \code{MinTmPr} Minimum melting temperature reverse primer
#'             \item \code{CGPf} Number of Cs or Gs in forward primer
#'             \item \code{CGPr} Number of Cs or Gs in reverse primer
#'             \item \code{Specificity} GG (Good-Good) means that both primer are specific to the target dataset,
#'                 GB or BG (Good-Bad or Bad-Good) means that only one of the two primers is specific to the target dataset
#'             \item \code{Amp_Records} Number of records in the target dataset that are amplified
#'             \item \code{P_Records} Proportion of records in the target dataset that are amplified
#'             \item \code{Amp_taxa} Number of taxa in the target dataset that are amplified
#'             \item \code{Amp_NT_taxa} Number of taxa in the non-target dataset that are amplified
#'             \item \code{Taxa_unique} Number of taxa with unique amplicons
#'             \item \code{P_taxa_unique} Proportion of taxa with unique amplicons
#'             \item \code{Min_length} Minimum amplicon length (excluding primers)
#'             \item \code{Max_length} Maximum amplicon length (excluding primers)
#'             \item \code{Mean_length} Mean amplicon length (excluding primers)}
#'     \item \code{metadata}: Parameters used during ecoPrimers command}
#' @note Specificty and Amp_NT_taxa have meaningless results unless a non-target dataset was specifiec during ecoPrimers command
#´
#' @examples
#' a<-"/media/sf_Documents/WORK/CIBIO/AA_PROJECTS/IRANVERTS/TESTING_PRIMERS_FOR_PREY/ovids_cytB/taxids.cytb_ovids_5-4-19.ecopcrdb"
#' b<-ecoPrimers.Bas(a,max_error = 1,min_length = 100,max_length = 300,strict_match_p = 0.7,e_match_p = 0.8)
#' @export
ecoPrimers.Bas<-function(ecopcrdb,max_error,min_length,max_length,strict_match_p,e_match_p,p3=T){
  cb <- function(line, proc) {cat(line, "\n")}
  a<-processx::run(command = "ecoPrimers",
                   args=c("-d",ecopcrdb,"-e",max_error,"-l",min_length,"-L", max_length,"-q",strict_match_p,"-s",e_match_p,"-3",p3),
                   echo=F,echo_cmd = T,stderr_line_callback = cb)
  
  f<-max(stringr::str_locate_all(string = a$stderr,pattern = ": ")[[1]][,2])
  g<-stringr::str_length(string = a$stderr)
  h<-as.numeric(substr(a$stderr,start = f+1,stop = g-1))
  message(paste("Primer search complete.",h,"primers found"))
  if(!h>0) stop("stop", call. = F)
  
  b<-data.table::fread(input = a$stdout,sep = "\t")
  ecoprimeroutput<-as.data.frame(b[,1:21])
  colnames(ecoprimeroutput)<- c("PID","Pf","Pr","TmPf","MinTmPf","TmPr","MinTmPr","CGPf","CGPr","Specificity","Amp_Records","IGNORE","P_Records",
                                "Amp_taxa","IGNORE2","Amp_NT_taxa","Taxa_unique","P_taxa_unique","Min_length","Max_length","Mean_length")
  ecoprimeroutput<-ecoprimeroutput[,-grep("IGNORE",colnames(ecoprimeroutput))]
  max(stringr::str_locate_all(string = a$stdout,pattern = "#")[[1]][,1])
  d<-substr(x = a$stdout,start = 1,stop=max(stringr::str_locate_all(string = a$stdout,pattern = "#")[[1]][,1]))
  ecoprimeroutput_META<-data.table::fread(input = d,sep = "\n")
  output_list <- list(ecoprimeroutput_META, ecoprimeroutput)
  names(output_list)<-c("metadata","primer_table")
  return(output_list)
}

#' Perform in silico PCR
#' @title Perform in silico PCR.
#' @param ecopcrdb An ecopcrdb containing the sequences to be tested
#' @param max_error Max number of mismatches allowed
#' @param min_length Minimum amplicon length (without primers - I think)
#' @param max_length Maximum amplicon length (without primers - I think)
#' @return A dataframe with in silico PCR results
#'
#' @export
ecoPCR.Bas<-function(Pf,Pr,ecopcrdb,max_error,min_length,max_length,out,buffer=NULL){
  #cb <- function(line, proc) {cat(line, "\n")}
  
  if(length(grep("I",Pf))>0)(Pf<-gsub("I","N",Pf))
  if(length(grep("I",Pr))>0)(Pf<-gsub("I","N",Pr))
  
  if(is.null(buffer)){
    system2(command = "ecoPCR",args=c(Pf, Pr,"-d",ecopcrdb,"-e",max_error,"-l",min_length,"-L", max_length,"-k","-c"),
            stdout = out,wait = T)}
  
  if(!is.null(buffer)){
    system2(command = "ecoPCR",
            args=c(Pf, Pr,"-d",ecopcrdb,"-e",max_error,"-l",min_length,"-L", max_length,"-k","-c","-D",buffer),
            stdout = out,wait = T)}
  
  #remove comment lines to help data.table to read properly
  f<-process$new(command = "sed", args = c("-i","/^#/ d", out), echo_cmd = T)
  f$wait()
  #remove weird hashes also
  f<-process$new(command = "sed", args = c("-i","s/###//g", out), echo_cmd = T)
  f$wait()
  #and hashes within names
  f<-process$new(command = "sed", args = c("-i","s/#/_/g", out), echo_cmd = T)
  f$wait()
  #remove quotes too
  f<-process$new(command = "sed", args = c("-i",'s/"//g', out), echo_cmd = T)
  f$wait()
  #remove quotes too
  f<-process$new(command = "sed", args = c("-i","s/'//g", out), echo_cmd = T)
  f$wait()
  
  b<-data.table::fread(file = out,data.table = F,
                       sep = "|",
                       col.names = c("AC","seq_length","taxid","rank","species","species_name","genus",
                                     "genus_name","family","family_name","superkingdom","superkingdom_name",
                                     "strand","forward_match","forward_mismatch","forward_tm","reverse_match",
                                     "reverse_mismatch","reverse_tm","amplicon_length","sequence","definition"))
  #remove any apparently erroneous entries
  #those that have no amplicon length or forward_mismatch or reverse_mismatch
  b<-b[!is.na(b$amplicon_length),]
  b<-b[!is.na(b$forward_mismatch),]
  b<-b[!is.na(b$reverse_mismatch),]
  
  if(length(b$AC)<1) stop("No hits",call. = F)
  write.table(b,file = out,quote = F,row.names = F,sep = "\t")
}

#modelling taxtabs
binarise.taxatab<-function(taxatab,t=F){
  #transpose to have species as columns
  taxatab2<-as.data.frame(t(taxatab[,-1]))
  colnames(taxatab2)<-taxatab[,1]
  #make binary
  taxatab2[taxatab2>0]<-1
  #if we want same format in output as in input
  if(t==T) {
    taxatab2<-as.data.frame(t(taxatab2))
    taxatab2$taxon<-rownames(taxatab2)
    taxatab2<-taxatab2[,c(length(colnames(taxatab2)),1:(length(colnames(taxatab2))-1))]
    rownames(taxatab2)<-NULL
  }
  return(taxatab2)
}

#make bray distance matrix
taxatab2bray<-function(binarised.taxatab){
  #convert binary matrix to distance matrix
  distance_matrix<-vegan::vegdist(binarised.taxatab,method = "bray",binary = T)
}

#' Download NCBI taxonomy
#' @title Download NCBI taxonomy.
#' @param path Full path to download NCBI taxonomy to. If using current directory use \code{path="here"}
#' @return An unpacked NCBI taxonomy
#' @note Any previous versions of the taxonomy located the same directory are overwritten
#' @examples
#' getNCBItaxonomy("here")
#' @export
getNCBItaxonomy<-function(path){
  if(path!="here"){setwd(path)}
  list.files()
  file.remove("taxdump.tar.gz")
  cb <- function(line, proc) {cat(line, "\n")}
  processx::run(command = "wget",
                args=c('ftp://ftp.ncbi.nlm.nih.gov://pub/taxonomy/taxdump.tar.gz'),echo=F,stderr_line_callback = cb)
  taxdump<-"taxdump.tar.gz"
  processx::run(command = "tar",
                args=c("-xzf",taxdump),echo=F,stderr_line_callback = cb)
  b<-"SUCCESS!"
}

#' Convert NCBI taxonomy to obitools taxonomy
#' @title Convert NCBI taxonomy to obitools taxonomy.
#' @param path Full path to directory containing the (unzipped) NCBI taxonomy files. If using current directory use \code{path="here"}
#' @param out A new name to give the obitaxonomydb
#' @return An obitools taxonomy
#' @examples
#' NCBI2obitaxonomy("here","obitax_26-4-19")
#' @export
NCBI2obitaxonomy<-function(path,out){
  if(path!="here"){setwd(path)}
  a<-getwd()
  cb <- function(line, proc) {cat(line, "\n")}
  processx::run(command = "obitaxonomy", args=c("-t",a,"-d",out),echo=F,stderr_line_callback = cb,echo_cmd = T)
  b<-"SUCCESS!"
}

#' Download megan accession 2 taxonomy map
#' @title Download megan accession 2 taxonomy map.
#' @param path Full path to download acc2tax to. If using current directory use \code{path="here"}
#' @return An unpacked acc2tax map
#' @note This file is c. 1gAny previous versions of the taxonomy located the same directory are overwritten??????
#' @examples
#' get_acc2tax_map("/media/sf_Documents/WORK/CIBIO/STATS_AND_CODE/TAXONOMIES/")
#' @export
get_acc2tax_map<-function(path){
  if(path!="here"){setwd(path)}
  cb <- function(line, proc) {cat(line, "\n")}
  a<-processx::run(command = "wget",
                   args=c("-nc","-r","-nd","--no-parent","-A",'nucl_acc2tax*',"--spider",
                          'http://ab.inf.uni-tuebingen.de/data/software/megan6/download/'),echo=F,spinner = T)
  e<-stringr::str_extract_all(string = a$stderr,pattern =  "nucl_acc2tax(.*?)zip")[[1]][1]
  message(c("Found file from megan6/download/: ",e))
  f<-strsplit(e,split = ".",fixed = T)[[1]][1]
  print(list.files())
  if(length(grep(x = list.files(),pattern = f))>0){
    stop("Extracted file already exists in folder, not downloading",call. = F)}
  
  message("downloading...will take some time (>1Gb)")
  d<-processx::run(command = "wget",
                   args=c("-nc","-r","-nd","--no-parent","-A",'nucl_acc2tax*',
                          'http://ab.inf.uni-tuebingen.de/data/software/megan6/download/'),echo=F,spinner = T)
  message("unpacking...")
  processx::run(command = "unzip", args=c(e),echo=F,stderr_line_callback = cb)
  file.remove(e)
  b<-"SUCCESS!"
}

#' Some colours I like to use for plots (n=28)
#'@export
MyCols <- c("dodgerblue2","#E31A1C", # red
            "green4",
            "#6A3D9A", # purple
            "#FF7F00", # orange
            "gold1",
            "skyblue2","#FB9A99", # lt pink
            "palegreen2",
            "#CAB2D6", # lt purple
            "#FDBF6F", # lt orange
            "gray70", "khaki2",
            "maroon","orchid1","deeppink1","blue1","steelblue4",
            "darkturquoise","green1","yellow4","yellow3",
            "darkorange4","brown","chartreuse2","coral4","darkslategray1","mediumseagreen")

#do xtabs to 1) find levels in factors and 2) see counts
master_xtabs<-function(master_sheet,columns){
  
  coln<-as.numeric(0)
  for(i in 1:length(columns)){
    coln[i]<-grep(paste0(columns[i],"$"), colnames(master_sheet))
  }
  
  xtabs(data = master_sheet[,coln],addNA = T)
}  

#subset a mastersheet based on required values
# e.g
# ms_ss<-subset_mastersheet(master_sheet,
#                           list(experiment_id=c("2018_02")
#                           ,Primer_set=c("12SV5.2","12SV51"),
#                           Sample_Type=c("Field")))
subset_mastersheet<-function(master_sheet,...){
  
  a<-as.list(...)
  
  print(a)
  
  m2<-master_sheet
  
  for(i in 1:length(a)){
    coln<-grep(paste0(names(a)[i],"$"), colnames(m2))
    m2<-m2[m2[,coln] %in% a[[i]],]
  }
  return(m2)
}

check.low.res.results<-function(pathofinterest,bins,btab){
  #first query otus that contributed to pathofinterest
  bins<-bins[bins$path==pathofinterest,"qseqid"]
  
  #next, query otus in btab
  btab<-btab[btab$qseqid %in% bins,]
  
  #find average pident for each contributor taxon
  contributors<-do.call(data.frame,aggregate(x = btab$pident,by=list(btab$path),
                                             FUN=function(x) c(mn = mean(x), n = range(x) )))
  colnames(contributors)<-c("contributors","mean.pident","low.pident","high.pident")
  
  #add pathofinterest for reference
  contributors$pathofinterest<-pathofinterest
  
  #add taxid for later
  btab2<-btab[!duplicated(btab$path),]
  btab3<-btab2[,c("taxids","path")]
  contributors<-merge(contributors,btab3,by.x = "contributors",by.y = "path")
  
  return(contributors)
}

check.low.res.df<-function(filtered.taxatab,filtered_blastfile, binfile,disabledTaxaFile=NULL,
                           spident=NULL,gpident=NULL,fpident=NULL,abspident=NULL){
  
  bins<-data.table::fread(binfile,sep = "\t",data.table = F)
  bins$path<-paste0(bins$K,";",bins$P,";",bins$C,";",bins$O,";",bins$F,";",bins$G,";",bins$S)
  taxatab.tf<-data.table::fread(filtered.taxatab,sep = "\t",data.table = F)
  taxatab.tf<-taxatab.tf[taxatab.tf$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]
  taxatab.tf<-taxatab.tf[taxatab.tf$taxon!="No_hits",]
  taxatab.tf<-taxatab.tf[taxatab.tf$taxon!="NA;NA;NA;NA;NA;NA;NA",]
  btab<-data.table::fread(file = filtered_blastfile,sep = "\t",data.table = F)
  btab$path<-paste0(btab$K,";",btab$P,";",btab$C,";",btab$O,";",btab$F,";",btab$G,";",btab$S)
  
  contributorlist<-list()
  for(i in 1:length(taxatab.tf$taxon)){
    contributorlist[[i]]<-check.low.res.results(pathofinterest = taxatab.tf$taxon[i],bins = bins,btab = btab)
  }
  
  contributordf<-do.call(rbind,contributorlist)
  
  #add number of reads and no. of "pcrs" pathofinterest
  taxatab.tf$readcounts<-rowSums(taxatab.tf[,2:length(colnames(taxatab.tf))])
  taxatab.tf$n.samples<-rowSums(taxatab.tf[,2:(length(colnames(taxatab.tf))-1)]>0)
  contributordf<-merge(contributordf,taxatab.tf[,c("taxon","readcounts","n.samples")],
                       by.x = "pathofinterest",by.y = "taxon")
  
  #add disabled taxa columns
  if(!is.null(disabledTaxaFile)){
    disabledTaxaDf<-read.table(disabledTaxaFile, header=T,sep = "\t")
    if(!"taxids" %in% colnames(disabledTaxaDf)) stop("No column called 'taxids'")
    if(!"disable_species" %in% colnames(disabledTaxaDf)) stop("No column called 'disable_species'")
    if(!"disable_genus" %in% colnames(disabledTaxaDf)) stop("No column called 'disable_genus'")
    if(!"disable_family" %in% colnames(disabledTaxaDf)) stop("No column called 'disable_family'")
    
    disabledSpecies<-disabledTaxaDf[disabledTaxaDf$disable_species==T,"taxids"]
    disabledSpecies<-disabledSpecies[!is.na(disabledSpecies)]
    
    disabledGenus<-disabledTaxaDf[disabledTaxaDf$disable_genus==T,"taxids"]
    disabledGenus<-disabledGenus[!is.na(disabledGenus)]
    
    disabledFamily<-disabledTaxaDf[disabledTaxaDf$disable_family==T,"taxids"]
    disabledFamily<-disabledFamily[!is.na(disabledFamily)]
    
    contributordf$species_disabled<-contributordf$taxids %in% disabledSpecies
    contributordf$genus_disabled<-contributordf$taxids %in% disabledGenus
    contributordf$family_disabled<-contributordf$taxids %in% disabledFamily
    
  } else {contributordf$species_disabled<-"FALSE"
  contributordf$genus_disabled<-"FALSE"
  contributordf$family_disabled<-"FALSE"
  }
  
  #add rank
  temprank<-stringr::str_count(contributordf$pathofinterest,";NA")
  temprank<-gsub(0,"species",temprank)
  temprank<-gsub(1,"genus",temprank)
  temprank<-gsub(2,"family",temprank)
  temprank<-gsub(3,"above_family",temprank)
  
  contributordf$rank<-temprank
  
  if(!is.null(spident) & !is.null(gpident) & !is.null(fpident) & !is.null(abspident)){
    #add may_be_improved
    max_pidents<-aggregate(contributordf$high.pident,by=list(contributordf$pathofinterest),FUN=max)
    best_is_species<-max_pidents[max_pidents$x>spident,]
    if(length(best_is_species$Group.1)>0) best_is_species$best_possible<-"species"
    max_pidents<-max_pidents[!max_pidents$x>spident,]
    best_is_genus<-max_pidents[max_pidents$x>gpident,]
    if(length(best_is_genus$Group.1)>0) best_is_genus$best_possible<-"genus"
    max_pidents<-max_pidents[!max_pidents$x>gpident,]
    best_is_family<-max_pidents[max_pidents$x>fpident,]
    if(length(best_is_family$Group.1)>0) best_is_family$best_possible<-"family"
    max_pidents<-max_pidents[!max_pidents$x>fpident,]
    best_is_above_family<-max_pidents[max_pidents$x>abspident,]
    if(length(best_is_above_family$Group.1)>0) best_is_above_family$best_possible<-"above_family"
    
    best_all<-rbind(best_is_species,best_is_genus,best_is_family,best_is_above_family)
    colnames(best_all)[1]<-"pathofinterest"
    best_all$x=NULL
    
    contributordf<-merge(contributordf,best_all)
    contributordf$may_be_improved<-contributordf$rank!=contributordf$best_possible
    contributordf<-contributordf[contributordf$may_be_improved=="TRUE",]
    contributordf<-contributordf[!contributordf$high.pident<abspident,]
    
  } else {(message("Not adding 'May be improved' column as no pidents provided"))}
  
  
  #if pathofinterest at genus level, dont output contributors from diferent genera
  contributordf$contributors<-as.character(contributordf$contributors)
  
  #1. make column to see whether genera match
  contributordf$contrGenus<-do.call(rbind,stringr::str_split(contributordf$contributors,";"))[,6]
  contributordf$pathGenus<-do.call(rbind,stringr::str_split(contributordf$pathofinterest,";"))[,6]
  contributordf$contr.path.genus.match<-contributordf$contrGenus==contributordf$pathGenus
  
  #2. If rank=genus, and columns dont match, then remove
  contributordf<-contributordf[!(contributordf$rank=="genus" & contributordf$contr.path.genus.match==FALSE),]
  
  #if pathofinterest at family level, dont output contributors from different families
  
  #1. make column to see whether families match
  contributordf$contrFam<-do.call(rbind,stringr::str_split(contributordf$contributors,";"))[,5]
  contributordf$pathFam<-do.call(rbind,stringr::str_split(contributordf$pathofinterest,";"))[,5]
  contributordf$contr.path.fam.match<-contributordf$contrFam==contributordf$pathFam
  
  #2. If rank=family, and columns dont match, then remove
  contributordf<-contributordf[!(contributordf$rank=="family" & contributordf$contr.path.fam.match==FALSE),]
  
  contributordf<-contributordf[order(contributordf$pathofinterest,-contributordf$high.pident),1:14]
  
  write.table(x = contributordf,file=gsub("spliced.txt","spliced.contr.txt",filtered.taxatab),
              sep="\t",quote = F,row.names = F)
}


#read mastersheets and make a processing sheet 
google.make.experiment.sheet<-function(outDir,sheeturls,experiment_id){
  master<-list()
  headers<-c("barcode_id","Primer_set","Primer_F","Primer_R","Min_length","Max_length","ss_sample_id","experiment_id")
  for(i in 1:length(sheeturls)){
    master[[i]]<-google.read.master.url(sheeturls[i])
    if(length(headers)!=sum(headers %in% colnames(master[[i]]))){
      stop (c("one of the following headers missing: ", paste(headers,collapse = " ")))}
    master[[i]]<-master[[i]][,headers]
  }
  
  #make a processing sheet
  experimentsheet<-as.data.frame(data.table::rbindlist(master))
  experimentsheet<-experimentsheet[experimentsheet$experiment_id==experiment_id,]
  
  #write file
  write.table(experimentsheet,paste0(outDir,experiment_id,"_experiment_sheet.txt"),sep = "\t",quote = F,row.names = F)
  message(paste("file saved as",paste0(outDir,experiment_id,"_experiment_sheet.txt")))
  
  return(experimentsheet)
}


filter.blast<-function(blastfile,headers="qseqid evalue staxid pident qcovs",ncbiTaxDir,out,min_qcovs=70,
                       max_evalue=0.001,top=1){
  
  if(length(grep("qcovs",headers))<1) stop("qcovs not in headers")
  if(length(grep("evalue",headers))<1) stop("evalue not in headers")
  if(length(grep("qseqid",headers))<1) stop("qseqid not in headers")
  if(length(grep("pident",headers))<1) stop("pident not in headers")
  if(length(grep("staxid",headers))<1) stop("staxid not in headers")
  
  if(is.null(ncbiTaxDir)) stop("ncbiTaxDir not specified")
  if(is.null(obitaxdb)) stop("obitaxdb not specified")
  if(is.null(out)) stop("out not specified")
  
  message("reading blast results")
  btab<-as.data.frame(data.table::fread(file = blastfile,sep = "\t"))
  colnames(btab)<-strsplit(headers,split = " ")[[1]]
  
  message("applying global min_qcovs threshold")
  btab<-btab[btab$qcovs>min_qcovs,]
  message("applying global max_evalue threshold")
  btab<-btab[btab$evalue<max_evalue,]
  message("applying global top threshold")
  if(length(btab[,1])==0) stop("No hits passing min_qcovs and max_evalue thresholds")
  topdf<-aggregate(x = btab[,colnames(btab) %in% c("qseqid","pident")],by=list(btab$qseqid),FUN = max)
  topdf$min_pident<-topdf$pident-topdf$pident*top/100
  btab<-merge(btab,topdf[,c(2,4)],by = "qseqid", all.y = T) #can definitely do this differently and faster
  btab<-btab[btab$pident>btab$min_pident,]
  
  #add lineage to results
  message("adding taxonomic lineages")
  btab$taxids<-btab$staxid #add.lineage.df requires this colname
  btab<-add.lineage.df(btab,ncbiTaxDir)
  
  #remove crappy hits 
  #1. btab$S contains uncultured
  message("Removing species containing the terms: uncultured, environmental, 
          unidentified,fungal, eukaryote or unclassified")
  if(length(grep("uncultured",btab$S,ignore.case = T))>0) {
    btab<-btab[-grep("uncultured",btab$S,ignore.case = T),]}
  #2. btab$S contains environmental
  if(length(grep("environmental",btab$S,ignore.case = T))>0) {
    btab<-btab[-grep("environmental",btab$S,ignore.case = T),]}
  #3. btab$S contains unclassified
  if(length(grep("unclassified",btab$S,ignore.case = T))>0) {
    btab<-btab[-grep("unclassified",btab$S,ignore.case = T),]}
  #4. btab$S contains "unidentified"
  if(length(grep("unidentified",btab$S,ignore.case = T))>0) {
    btab<-btab[-grep("unidentified",btab$S,ignore.case = T),]}
  #5. btab$S contains "fungal "
  if(length(grep("fungal ",btab$S,ignore.case = T))>0) {
    btab<-btab[-grep("fungal ",btab$S,ignore.case = T),]}
  #6. btab$S contains "eukaryote"
  if(length(grep("eukaryote",btab$S,ignore.case = T))>0) {
    btab<-btab[-grep("eukaryote",btab$S,ignore.case = T),]}
  
  write.table(x = btab,file = out,sep="\t",quote = F,row.names = F)
  
}


#loop glm by taxon
loop.glm<-function(taxatab, master_sheet,factor,grouping,summaries=F){
  if(!"ss_sample_id" %in% colnames(master_sheet)) stop("No column called ss_sample_id")
  message("Currently, response is binary detections in PCR reps.
          For each model, only PCR replicates belonging to a 
          group that had at least one taxon detection are included")
  taxatab<-rm.0readtaxSam(taxatab)
  master_sheet2<-master_sheet[master_sheet$ss_sample_id %in% colnames(taxatab[,-1]),]
  #should only run models for levels in group that had taxon
  
  model.list<-list()
  for(i in 1:length(taxatab$taxon)){
    taxatab2<-taxatab[taxatab$taxon==taxatab$taxon[i],]
    taxatab3<-binarise.taxatab(taxatab2)
    taxatab3$ss_sample_id<-rownames(taxatab3)
    taxatab4<-merge(taxatab3,master_sheet2[,c("ss_sample_id",factor,grouping)],by="ss_sample_id")
    
    #find groups with dxns
    pos.groups<-aggregate(taxatab4[,as.character(taxatab$taxon[i])],by=list(taxatab4[,grouping]),FUN=sum)
    pos.groups<-pos.groups[pos.groups$x>0,]
    
    taxatab4<-taxatab4[taxatab4[,grouping] %in% pos.groups$Group.1,]
    
    species<-stringr::str_split(taxatab$taxon[i],";")[[1]][7]
    
    message(paste("Running glm with",factor,"as predictor and binary detection of",species,  "as response"))
    
    message(paste0("number of ",factor,"s where taxon was detected at least once = ",length(unique(taxatab4[,factor]))))
    if(length(unique(taxatab4[,factor]))<2) {
      model.list[[i]]<-"Cannot run this model as there are less than two groups"
      names(model.list)[i]<-species 
      message("Cannot run this model as there are less than two groups")
    } else{
      message("number of PCRs in model = ",length(taxatab4$ss_sample_id))
      
      model1<-glm(taxatab4[,as.character(taxatab$taxon[i])]~taxatab4[,factor], family = "binomial")
      if(summaries) print(summary(model1))
      a<-capture.output(interp.glm.bin(model1))
      b<-gsub("taxatab4\\[, factor\\]",factor,a)
      message(gsub("taxatab4\\[, as.character\\(taxatab\\$taxon\\[i\\]\\)\\]",species,b))
      
      print(aggregate(taxatab4[,as.character(taxatab$taxon[i])],by=list(taxatab4[,factor]),FUN=sum))
      
      model.list[[i]]<-model1
      names(model.list)[i]<-species
    }
  }
  
  return(model.list)
  
}
