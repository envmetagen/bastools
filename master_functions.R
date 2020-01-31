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
  
  all.taxatabs<-all.taxatabs[order(all.taxatabs$taxon),]
  
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
  message(paste("Using detection filter of",filter_dxn, ": reads removed:",reads1-reads2,"from",reads1, "; detections removed:",dxns1-dxns2,"from",dxns1))
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
remove.contaminant.taxa<-function(master_sheet,taxatab,negatives,group.codes,printcontaminations=T){
  
  totalstart<-sum(taxatab[,-1,drop=F])
  
  for(j in 1:length(negatives)){
    
    negdf<-negatives[[j]]
    
    negative_type<-names(negatives)[j]
    
    if(length(negdf)!=0){
      
      for(i in 2:length(colnames(negdf))){
        neg<-colnames(negdf)[i]
        negtaxa<-negdf[,c("taxon",neg)]
        sumnegtaxa<-sum(negtaxa[,-1,drop=F])
        negtaxa<-as.character(negtaxa[rowSums(negtaxa[,-1,drop=F])>0,"taxon"])
        
        #get group of a sample
        group.id<-master_sheet[grep(neg,master_sheet[,"ss_sample_id"]),group.codes[j]]
        #get other samples in group
        group.samples<-master_sheet[master_sheet[,group.codes[j]]==group.id,"ss_sample_id"]
        group.samples<-group.samples[!is.na(group.samples)]
        
        
        message(paste("Based on",negative_type,neg, ", Removing detections of"))
        print(negtaxa)
        message("if it occurred in any of the following samples")
        print(group.samples)
        message(paste("which all belong to",group.codes[j],":",group.id))
        
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

#keep only xLevel assignments
keep.below.xLevel.assigns<-function(taxatab,xLevel="species"){
  message("Broken without full 7-level path ")
  
  if(xLevel=="species"){
    taxatab<-taxatab[-grep(";NA$",taxatab$taxon),]
    taxatab<-rm.0readtaxSam(taxatab)
  }
  
  if(xLevel=="genus"){
    taxatab<-taxatab[-grep(";NA;NA$",taxatab$taxon),]
    taxatab<-rm.0readtaxSam(taxatab)
  }
  
  if(xLevel=="family"){
    taxatab<-taxatab[-grep(";NA;NA;NA$",taxatab$taxon),]
    taxatab<-rm.0readtaxSam(taxatab)
  }
  
  return(taxatab)
}

aggregate.at.xLevel<-function(taxatab,xLevel){
  
  if(!xLevel %in% c("genus","family","order")) stop("Only allowable at genus, family or order level")
  
  splittaxonomy<-as.data.frame(do.call(rbind,stringr::str_split(taxatab[,1],";")))
  
  if(xLevel=="genus"){
    xPath=paste0(splittaxonomy[,1],";",splittaxonomy[,2],";",splittaxonomy[,3],";",splittaxonomy[,4],";",splittaxonomy[,5],
                 ";",splittaxonomy[,6])
  }
  
  if(xLevel=="family"){
    xPath=paste0(splittaxonomy[,1],";",splittaxonomy[,2],";",splittaxonomy[,3],";",splittaxonomy[,4],";",splittaxonomy[,5])
  }
  
  if(xLevel=="order"){
    xPath=paste0(splittaxonomy[,1],";",splittaxonomy[,2],";",splittaxonomy[,3],";",splittaxonomy[,4])
  }
  
  taxatab<-aggregate(taxatab[,-1],by = list(xPath),FUN=sum)
  
  colnames(taxatab)[1]<-"taxon"
  
  taxatab
  
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

add.lineage.df<-function(df,ncbiTaxDir,as.taxids=F){
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
  if(as.taxids==F){
    system2(command = "taxonkit",args =  c("reformat",taxids_fileB,"-i",3,"--data-dir",ncbiTaxDir)
          ,stdout = taxids_fileC,stderr = "",wait = T)
  } else {
    system2(command = "taxonkit",args =  c("reformat","-t", taxids_fileB,"-i",3,"--data-dir",ncbiTaxDir)
            ,stdout = taxids_fileC,stderr = "",wait = T)
  }
  lineage<-as.data.frame(data.table::fread(file = taxids_fileC,sep = "\t"))
  colnames(lineage)<-gsub("V1","taxids",colnames(lineage))
  colnames(lineage)<-gsub("V2","new_taxids",colnames(lineage))
  if(as.taxids==F){
    colnames(lineage)<-gsub("V5","path",colnames(lineage))
  } else {
    colnames(lineage)<-gsub("V6","path",colnames(lineage))
  }

  #merge with df
  #message("replacing taxids with updated taxids. Saving old taxids in old_taxids.")
  df<-merge(df,lineage[,c("taxids","new_taxids","path")],by = "taxids")
  df$old_taxids<-df$taxids
  df$taxids<-df$new_taxids
  df$new_taxids=NULL
  df<-cbind(df,do.call(rbind, stringr::str_split(df$path,";")))
  colnames(df)[(length(df)-6):length(df)]<-c("K","P","C","O","F","G","S")
  if(as.taxids==F){
    df$K<-as.character(df$K)
    df$P<-as.character(df$P)
    df$C<-as.character(df$C)
    df$O<-as.character(df$O)
    df$F<-as.character(df$F)
    df$G<-as.character(df$G)
    df$S<-as.character(df$S)
  
  #change empty cells to "unknown"
  df[,(length(df)-6):length(df)][df[,(length(df)-6):length(df)]==""]<- "unknown"
  }
  
  df$path=NULL
  unlink(taxids_fileA)
  unlink(taxids_fileB)
  unlink(taxids_fileC)
  return(df)
}

get.children.taxonkit<-function(df,column){
  if(is.null(df$taxids)) {stop("No column called taxids")}
  df$taxids<-as.integer(as.character(df$taxids)) 
  #write taxids to file
  taxids_fileA<-paste0("taxids",as.numeric(Sys.time()),".txt")
  write.table(unique(df$taxids),file = taxids_fileA,row.names = F,col.names = F,quote = F)
  
  #get children 
  taxids_fileB<-paste0("taxids",as.numeric(Sys.time()),".txt")
  system2(command = "taxonkit",args =  c("list","--ids",paste(as.character(df$taxids),collapse = ","),"--indent", "''","--data-dir",ncbiTaxDir)
          ,stdout = taxids_fileB,stderr = "",wait = T)
  
  children<-data.table::fread(file = taxids_fileB,sep = "\t",data.table = F)
  children<-children[!is.na(children$V1),]
  
  unlink(taxids_fileA)
  unlink(taxids_fileB)
  return(children)
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
    a[[i]]<-tryCatch(h[[i]]$get_status(),error=function(e) print("not running, may be finished"))
    if(nchar(a[[i]])!=28) message(gsub("sleeping","running",h[[i]]$get_status()))
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
  
  # #import counts step4 
  #cant figure out how to do this in R!
  current<-getwd()
  setwd(MBCtsvDir)
  tempname<-paste0(as.numeric(Sys.time()),".txt")
  writeLines(c("#!/bin/bash\nfind . -name *none.flash2_merged.vsearch_qfilt.cutadapt.vsearch_uniq.fasta.gz -print0 | xargs -0 zgrep '>' | sed 's/\\.\\///;s/.*\\///;s/.none.*size=/ /'")
             ,paste0(MBCtsvDir,tempname))
  step4<-as.data.frame(system2(command = "sh",args = paste0(MBCtsvDir,tempname),stdout = T,wait = T))
  unlink(paste0(MBCtsvDir,tempname))
  setwd(current)
  colnames(step4)<-"V1"
  step4$counts<-as.numeric(do.call(rbind,stringr::str_split(step4$V1," "))[,2])
  step4$samples<-do.call(rbind,stringr::str_split(step4$V1," "))[,1]
  step4<-aggregate(step4$counts,by=list(step4$samples),FUN=sum)
  colnames(step4)<-c("ss_sample_id","nseqs")
  after_rm.singletons<-sum(step4[step4$ss_sample_id %in% ms_ss$ss_sample_id,"nseqs"])
  
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
  if("NA;NA;NA;NA;NA;NA;NA" %in% taxatab_C$taxon){
    sumNAs<-sum(taxatab_C[taxatab_C$taxon=="NA;NA;NA;NA;NA;NA;NA",-1])
  } else {sumNAs<-0}
  after_blast_filt<-after_blast-sumNAs
  
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
                  "After singleton removal"=after_rm.singletons,
                  "After size selection"=after_size_select,
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


make.blastdb.bas<-function(infasta,makeblastdb_exec="makeblastdb",addtaxidsfasta=F, ncbiTaxDir, dbversion=5){
  library(processx)
  
  message("Reminder: Assumes header includes taxid=taxid;")
  message("Reminder: Assumes there are only spaces between attributes, not within them, e.g. species names should not have spaces")
  
  
  message("Reminder: Only works for blast 2.9.0: Current version:")
  system2(command = makeblastdb_exec,args = c("-version"))
  if(!infasta %in% list.files()) stop("infasta not found in current directory")

  tempfasta<-phylotools::read.fasta(infasta)
  
  tempfasta$ids<-do.call(rbind,strsplit(as.character(tempfasta$seq.name)," "))[,1]
  taxids<-do.call(rbind,strsplit(as.character(tempfasta$seq.name),"taxid="))[,2]
  tempfasta$taxids<-do.call(rbind,strsplit(as.character(taxids),";"))[,1]
  
  #give unique ids
  if(length(unique(tempfasta$ids))!=length(tempfasta$ids)) {
    message("IDs not unique - giving unique ids - ensure to check if mapping back!")
    tempfasta$ids<-make.unique(tempfasta$ids, sep = ".")
  }

  #Add "database name to header
  message("Adding db name to headers")
  tempfasta$db<-gsub(".fasta","",infasta)
  
  #ensure ids are <50 characters
  message("Ensuring ids are <50 characters - ensure to check if mapping back!")
  tempfasta$ids<-stringr::str_trunc(as.character(tempfasta$ids),width = 49)
  #remove any quotes
  tempfasta$ids<-gsub('"',"",tempfasta$ids)
  
  #make final headers & write fasta
  tempfasta$seq.name<-paste0(tempfasta$ids," taxids=",tempfasta$taxids, "; db=",tempfasta$db)
  phylotools::dat2fasta(tempfasta,gsub(".fasta",".blastdbformatted.fasta",infasta))
  
  #make mapping file
  mappingfile<-paste0("mapping",as.numeric(Sys.time()),".txt")
  write.table(tempfasta[,c("ids","taxids")],mappingfile,quote = F,sep = " ",row.names = F,col.names = F)
  
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
  message("Does new db have taxids - column V3?")
  print(data.table::fread(gsub(".fasta",".blastdbformatted.test.blast.txt",infasta)))
  system2(command = "blastdbcheck",args = c("-must_have_taxids","-db",gsub(".fasta","",infasta)))
  
  unlink(gsub(".fasta",".blastdbformatted.test.blast.txt",infasta))
  unlink(gsub(".fasta",".blastdbformatted.test.fasta",infasta))
  unlink(mappingfile)
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
obiconvert.Bas<-function(infile,in_type,taxo,out_type,out,add2existing=F){
  
  if(in_type!="fasta" & out_type=="--ecopcrdb-output") stop("Only fasta files can be converted into an ecopcrdb")
  
  cb <- function(line, proc) {cat(line, "\n")}
  if(out_type=="--fasta-output"){
    b<-processx::run(command = "obiconvert", args=c(infile,"-d",taxo,out_type),echo=F,stderr_line_callback = cb,echo_cmd = T)
    writeLines(b$stdout,con = out)
  }
  
  if(out_type=="--ecopcrdb-output"){
    message("Make sure input fasta file has taxids in sequence headers in the format \" taxid=XXXX;\"")
    if(add2existing==F) {
      if(length(grep(paste0(out,".adx"),list.files())>0)) {
        message(paste("Removing existing datase", out))
        unlink(x = list.files(pattern = out)) 
      }
    }
    
    tempfasta<-paste0(as.numeric(Sys.time()),".fasta")
        #remove "|"s to stop ecopcr errors
    f<-process$new(command = "sed", args = c("s/|/_/g",infile), echo_cmd = T,
                   stdout=tempfasta)
    f$wait()
    system2(command = "obiconvert", args=c(tempfasta,"-d",taxo,paste0(out_type,"=",out),"--skip-on-error"),wait = T)
    d<-"SUCCESS!"
    unlink(tempfasta)
  }
  
  #count input fasta
  l<-paste0("input_read_count=",bascount.fasta(infile))
  #count output ecopcr
  h<-processx::run(command = "obicount", args=out,echo=F,echo_cmd = F,stderr_line_callback = cb)
  j<-stringr::str_locate_all(string = h$stdout,pattern = " " )[[1]][,1]
  k<-substr(h$stdout,1,j-1)
  m<-paste0("output_read_count=",k)
  
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
  taxa<-do.call(rbind,stringr::str_split(taxatab$taxon,";"))
  taxa<-cbind(taxa,do.call(rbind,stringr::str_split(taxa[,7]," ")))
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
  

  if(as.dxns) long<-long[!duplicated(long),]
  
  a<-ggplot2::ggplot(data=long , aes(y=value, x=variable, fill=taxon))+
    theme(legend.title = element_text(size=10), legend.text=element_text(size=10),
          axis.text.x=element_text(size=8,angle=45, hjust=1),legend.position="right",legend.direction="vertical")+
    ggtitle(paste("x axis=",column,"; as.percent=",as.percent,"; as.dxns=",as.dxns, ";facetcol=",facetcol))
  
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
  
  if(MS_colID!="ss_sample_id") {
    taxatab2<-sumreps(taxatab = taxatab,master_sheet = master_sheet,by = MS_colID) 
  } else {
    taxatab2<-taxatab
    }
  
  taxatab2<-rm.0readtaxSam(taxatab2)
  
  taxatab2<-binarise.taxatab(taxatab2)
  
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
  #cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  p<-ggplot(cmdspoints,aes(x=V1,y=V2))+
    geom_point(aes(size=1,colour=cmdspoints[,factor1],stroke=1))+
    #geom_point(aes(size=1,colour=cmdspoints[,factor1],shape=cmdspoints[,factor1],stroke=1))+
    #geom_jitter(position = )
    #scale_shape_manual(values=c(1, 2, 0, 5, 6, 3, 4))+
    #scale_color_manual(values = c(cbPalette))+
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
    message("Note: Ellipses will not be calculated if there are groups with too few data points")
    p<- p+stat_ellipse(aes(colour=cmdspoints[,factor1] 
                           ,fill=cmdspoints[,factor1]
    )
    ,type = "norm", level=0.90, 
    geom = "polygon",alpha=0.2,
    show.legend = F,segments = 100) #+
      #scale_fill_manual(values = c(cbPalette))
    
    p$layers<-rev(p$layers)
    
    message("ellipses are drawn with a confidence level of 0.90")
  }
  
  if(lines){
    p<- p +geom_segment(data = loadings, aes(x=0,y=0,xend=V1,yend=V2),arrow=arrow(length=unit(0.1,"cm")))
    if(longnames) if(shortnames) stop("Can only use EITHER long OR short names")
    if(!longnames) if(!shortnames) message("No names added")
    if(longnames) if(!shortnames) p<- p + geom_text(data = loadings, aes(x=V1, y=V2, label=colnames(taxatab2)))
    if(shortnames) if(!longnames){
      taxa<-do.call(rbind,stringr::str_split(colnames(taxatab2),";"))
      taxa<-cbind(taxa,do.call(rbind,stringr::str_split(taxa[,7]," ")))
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
ecoPCR.Bas<-function(Pf,Pr,ecopcrdb,max_error,min_length,max_length=NULL,out,buffer=NULL){
  #cb <- function(line, proc) {cat(line, "\n")}
  
  if(length(grep("I",Pf))>0)(Pf<-gsub("I","N",Pf))
  if(length(grep("I",Pr))>0)(Pf<-gsub("I","N",Pr))
  
  if(!is.null(max_length)){
  if(is.null(buffer)){
    system2(command = "ecoPCR",args=c(Pf, Pr,"-d",ecopcrdb,"-e",max_error,"-l",min_length,"-L", max_length,"-k"),
            stdout = out,wait = T)}
  
  if(!is.null(buffer)){
    system2(command = "ecoPCR",
            args=c(Pf, Pr,"-d",ecopcrdb,"-e",max_error,"-l",min_length,"-L", max_length,"-k","-D",buffer),
            stdout = out,wait = T)}
  }
  
  if(is.null(max_length)){
    if(is.null(buffer)){
      system2(command = "ecoPCR",args=c(Pf, Pr,"-d",ecopcrdb,"-e",max_error,"-l",min_length,"-k"),
              stdout = out,wait = T)}
    
    if(!is.null(buffer)){
      system2(command = "ecoPCR",
              args=c(Pf, Pr,"-d",ecopcrdb,"-e",max_error,"-l",min_length,"-k","-D",buffer),
              stdout = out,wait = T)}
  }
  
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
    disabledTaxaDf<-data.table::fread(disabledTaxaFile, data.table = T,sep = "\t")
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
    
    #remove duplicated columns
    if(length(grep("\\.\\.\\.",colnames(master[[i]]),value = T))>0) {
      message("Removing duplicated columns")
      a<-grep("\\.\\.\\.",colnames(master[[i]]),value = T)
      b<-do.call(rbind,stringr::str_split(a,"\\.\\.\\."))
      d<-b[duplicated(b[,1]),]
      e<-paste0(d[,1],"...",d[,2])
      master[[i]]<-master[[i]][,!colnames(master[[i]]) %in% e]
      colnames(master[[i]])<-gsub("\\.\\.\\..*","",colnames(master[[i]]))
    }
    
    if(length(headers)!=sum(headers %in% colnames(master[[i]]))){
      stop (c("one of the following headers missing: ", paste(headers,collapse = " ")))
      }
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

taxatab.pieplot<-function(taxatab,hidelegend=F){
  require(plotly)
  a<-summary.dxns.by.taxon(taxatab)
  
  b<-plotly::plot_ly(a, labels = ~taxon, values = ~total.reads, type = 'pie',textposition = 'outside',textinfo = 'label+percent')%>%
       layout(showlegend = FALSE)

  # b<-ggplot(a, aes(x="", y=total.reads, fill=taxon))+
  #   geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)
  # 
  # if(nrow(a)<29) b <- b+scale_fill_manual(values = MyCols)
  # 
  # b+ggrepel::geom_text_repel(aes(x = 1.4, y = 1.4, label = taxon), 
  #                   nudge_x = .3, 
  #                   segment.size = .7, 
  #                   show.legend = FALSE) 
  # 
  # if(hidelegend) {
  #   message("Outputting as a list where first element is plot and second is legend")
  #   plotlist<-list()
  #   plotlist[[1]]<-b+theme(legend.position="none")
  #   plotlist[[2]]<-ggpubr::as_ggplot(cowplot::get_legend(b))
  #   return(plotlist)
  # }else return(b)
  b
}


reads.by.rank<-function(taxatab){
  a<-summary.dxns.by.taxon(taxatab)
  a<-a[a$taxon!="NA;NA;NA;NA;NA;NA;NA",]
  a<-a[a$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]
  
  #add rank
  temprank<-stringr::str_count(a$taxon,";NA")
  temprank<-gsub(0,"species",temprank)
  temprank<-gsub(1,"genus",temprank)
  temprank<-gsub(2,"family",temprank)
  temprank<-gsub(3,"above_family",temprank)
  
  a$rank<-temprank
  
  d<-aggregate(a$total.reads,by=list(a$rank),FUN=sum)
  
  colnames(d)<-c("rank","reads")
  
  d$percent<-round(d$reads/sum(d$reads)*100,digits = 1)
  
  return(d)
}

dxns.by.rank<-function(taxatab){
  a<-summary.dxns.by.taxon(taxatab)
  a<-a[a$taxon!="NA;NA;NA;NA;NA;NA;NA",]
  a<-a[a$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]
  
  #add rank
  temprank<-stringr::str_count(a$taxon,";NA")
  temprank<-gsub(0,"species",temprank)
  temprank<-gsub(1,"genus",temprank)
  temprank<-gsub(2,"family",temprank)
  temprank<-gsub(3,"above_family",temprank)
  
  a$rank<-temprank
  
  d<-aggregate(a$n.samples,by=list(a$rank),FUN=sum)
  
  colnames(d)<-c("rank","dxns")
  
  d$percent<-round(d$dxns/sum(d$dxns)*100,digits = 1)
  
  return(d)
}

sumreps<-function(taxatab,master_sheet,by="Sample_Name"){

  grouping<-master_sheet[match(colnames(taxatab[,-1]),table = master_sheet[,"ss_sample_id"]),by]
  
  taxatab2<-taxatab
  colnames(taxatab2)<-c("taxon",grouping)
  
  summed<-list()
  for(i in 1:length(unique(grouping))){
    taxatab.temp<-as.data.frame(taxatab2[,grep(unique(grouping)[i],colnames(taxatab2)),drop=F])
    if(length(colnames(taxatab.temp))>1) taxatab.temp$sum<-rowSums(taxatab.temp) else taxatab.temp$sum<-taxatab.temp[,1]
    
    taxatab.temp2<-taxatab.temp[,"sum",drop=F]
    
    colnames(taxatab.temp2)<-gsub("sum",unique(grouping)[i],colnames(taxatab.temp2))
    
    summed[[i]]<-taxatab.temp2
  }
  
  taxatab.out<-cbind(taxon=taxatab$taxon,do.call(cbind, summed))
}

names2taxids<-function(vector,ncbiTaxDir){
  message("Reminder - this function may require user input, DO NOT RUN LINES OF SCRIPT AFTER THIS FUNCTION")
  names_fileA<-paste0("names",as.numeric(Sys.time()),".txt")
  vector<-as.character(vector)
  #dont search names with sp. in the first round
  if(length(grep(" sp\\.",vector))>0){
    vector2<-vector[-grep(" sp\\.",vector)]} else {vector2<-vector}
  write.table(unique(vector2),file = names_fileA,row.names = F,col.names = F,quote = F) 
  taxids_fileA<-gsub("names","taxids",names_fileA)
  system2(command = "taxonkit", args = c("name2taxid",names_fileA,"-r","--data-dir",ncbiTaxDir),
          stdout = taxids_fileA,stderr = "",wait = T)
  taxidsA<-read.table(taxids_fileA,sep = "\t")
  
  #for unfound taxids, repeat search using only first word
  namesB<-c(as.character(taxidsA$V1)[is.na(taxidsA$V2)],vector[grep(" sp\\.",vector)])
  if(length(namesB)>0){
    namesB<-gsub(" .*","",namesB)
    names_fileB<-paste0("names",as.numeric(Sys.time()),".txt")
    taxids_fileB<-gsub("names","taxids",names_fileA)
    write.table(unique(namesB),file = names_fileB,row.names = F,col.names = F,quote = F)
    system2(command = "taxonkit", args = c("name2taxid",names_fileB,"-r","--data-dir",ncbiTaxDir),
            stdout = taxids_fileB,stderr = "",wait = T)
    taxidsB<-read.table(taxids_fileB,sep = "\t")}
  
  #Note that if a name matches more than one taxid, taxonkit creates a new row and includes both taxids, 
  #allow user input for choices
  message("Should add choices for unknowns")
  
  #first for good results in taxidsA
  taxidsA2<-taxidsA[!is.na(taxidsA$V2),]
  #find names with multiple matches
  if(length(taxidsA2$V1[duplicated(taxidsA2$V1)])>0){
    
    taxidsA4<-taxidsA2[duplicated(taxidsA2$V1),]
    taxidsA5<-taxidsA2[taxidsA2$V1 %in% taxidsA4$V1,]
    colnames(taxidsA5)<-c("name","taxids","rank")
    #add lineage to make choice easier
    taxidsA5<-add.lineage.df(taxidsA5,ncbiTaxDir)
    taxidsA5$old_taxids=NULL
    #present choice
    message("The following ", length(unique(taxidsA5$name))," taxa had more than one taxid match. Please
            type the full taxid of your choice")
    choicesA_df<-as.data.frame(unique(taxidsA5$name))
    for(i in 1:length(unique(taxidsA5$name))){
      print(taxidsA5[taxidsA5$name %in%  unique(taxidsA5$name)[i],])
      choicesA_df[i,2] <- readline("Type full chosen taxid: ")  
    }
    choicesA_df$V3<-"NA"
    colnames(choicesA_df)<-c("V1","V2","V3")
    #combine good results
    taxidsA6<-taxidsA2[!taxidsA2$V1 %in% unique(taxidsA5$name),]
    taxidsA7<-rbind(taxidsA6,choicesA_df)
  } else {taxidsA7<-taxidsA2}
  
  
  if(length(namesB)>0){
    #for good results in taxidsB
    taxidsB2<-taxidsB[!is.na(taxidsB$V2),]
    
    #find names with multiple matches
    if(length(taxidsB2$V1[duplicated(taxidsB2$V1)])>0){
      taxidsB4<-taxidsB2[duplicated(taxidsB2$V1),]
      taxidsB5<-taxidsB2[taxidsB2$V1 %in% taxidsB4$V1,]
      colnames(taxidsB5)<-c("name","taxids","rank")
      #add lineage to make choice easier
      taxidsB5<-add.lineage.df(taxidsB5,ncbiTaxDir)
      taxidsB5$old_taxids=NULL
      
      #if genus vs subgenus, choose genus
      choosegenusdf<-as.data.frame(unique(taxidsB5$name))
      for(i in 1:length(unique(taxidsB5$name))){
        tempdf<-taxidsB5[taxidsB5$name %in%  unique(taxidsB5$name)[i],]
        if(TRUE %in% duplicated(tempdf[4:10])){
          if(length(tempdf$rank)==2){
            if(length(tempdf$rank[tempdf$rank=="genus"])==1){
              if(length(tempdf$rank[tempdf$rank=="subgenus"])==1){
                choosegenusdf[i,2] <- tempdf$taxids[tempdf$rank=="genus"]
              }
            }
          }
        }
      }
      colnames(choosegenusdf)<-c("name","taxid")
      
      taxidsB5_remaining<-taxidsB5[!taxidsB5$name %in% 
                                     choosegenusdf[!is.na(choosegenusdf$taxid),"name"],]
      
      #for the remaining options present a choice
      message("The following ", length(unique(taxidsB5_remaining$name))," taxa had more than one taxid match. Please
              type the full taxid of your choice")
      choicesB_df<-as.data.frame(unique(taxidsB5_remaining$name))
      for(i in 1:length(unique(taxidsB5_remaining$name))){
        print(taxidsB5_remaining[taxidsB5_remaining$name %in%  unique(taxidsB5_remaining$name)[i],])
        choicesB_df[i,2] <- readline("Type full chosen taxid: ")  
      }
      choicesB_df$V3<-"NA"
      colnames(choicesB_df)<-c("V1","V2","V3")
      #combine good results
      taxidsB6<-taxidsB2[!taxidsB2$V1 %in% unique(taxidsB5_remaining$name),]
      taxidsB7<-rbind(taxidsB6,choicesB_df)
    } else {taxidsB7<-taxidsB2}
    
  }
  
  #combine taxidA, taxidB:choosegenus and taxidB:remaining results
  outdf<-as.data.frame(vector)
  outdf$vec2<-gsub(" sp\\.","",outdf$vector)
  outdf$vec3<-gsub(" .*","",outdf$vector)
  outdf<-merge(outdf,taxidsA7,by.x = "vec2",by.y = "V1",all.x = T,all.y = F)
  if(length(namesB)>0){
    outdf<-merge(outdf,taxidsB7,by.x = "vec3",by.y = "V1",all.x = T,all.y = F)
    outdf$V2.z<-ifelse(!is.na(outdf$V2.x) & !is.na(outdf$V2.y) | is.na(outdf$V2.y),outdf$V2.z<-outdf$V2.x,NA)}
  #add choosegenus
  if(length(namesB)>0){
    outdf<-merge(outdf,choosegenusdf,by.x = "vec3",by.y = "name",all.x = T,all.y = F)
    outdf$V2.w<-ifelse(is.na(outdf$V2.x) & is.na(outdf$V2.z),outdf$V2.w<-outdf$taxid,NA)}
  
  #combine taxids into one column
  if(length(namesB)>0){
    outdf$taxids<-do.call(pmax, c(outdf[,c("V2.z","V2.y","V2.w")], list(na.rm=TRUE)))} else {outdf$taxids<-outdf$V2}
  
  #remove files
  unlink(names_fileA)
  unlink(taxids_fileA)
  if(length(namesB)>0){
    unlink(names_fileB)
    unlink(taxids_fileB)}
  
  #remove extraneous columns
  colnames(outdf)<-gsub("vector","name",colnames(outdf))
  outdf<-outdf[match(vector, outdf$name),]
  outdf<-outdf[,"taxids"]
}

taxids2names<-function(df,ncbiTaxDir){
  taxids_fileA<-paste0("taxids",as.numeric(Sys.time()),".txt")
  write.table(unique(df$taxid),file = taxids_fileA,row.names = F,col.names = F,quote = F)
  taxids_fileB<-paste0("taxids",as.numeric(Sys.time()),".txt")
  g<-process$new(command = "taxonkit", args = c("lineage","-r",taxids_fileA,"--data-dir",ncbiTaxDir),
                 echo_cmd = T,stdout = taxids_fileB)
  g$wait()
  lineage<-read.table(taxids_fileB,sep = "\t")
  
  lineage$V3<-gsub("infraclass","class",lineage$V3)
  lineage$V3<-gsub("subfamily","family",lineage$V3)
  lineage$V3<-gsub("tribe","family",lineage$V3)
  lineage$V3<-gsub("suborder","order",lineage$V3)
  lineage$V3<-gsub("infraorder","order",lineage$V3)
  lineage$V3<-gsub("superfamily","order",lineage$V3)
  lineage$V3<-gsub("subclass","class",lineage$V3)
  lineage$V3<-gsub("subgenus","genus",lineage$V3)
  lineage$V3<-gsub("species subgroup","species",lineage$V3)
  lineage$V3<-gsub("cohort","order",lineage$V3)
  
  df2<-merge(df,lineage[,c(1,3)],by.x = "taxid",by.y = "V1", all.y = T)
  colnames(df)<-gsub("V3","rank",colnames(df))
  df<-cbind(df,do.call(rbind, stringr::str_split(df$path,";")))
  colnames(df)[(length(df)-6):length(df)]<-c("K","P","C","O","F","G","S")
  df[,(length(df)-6):length(df)] <- sapply(df[,(length(df)-6):length(df)],as.character)
  df[,(length(df)-6):length(df)][df[,(length(df)-6):length(df)]==""]<- "unknown"
  
  unlink(taxids_fileA)
  unlink(taxids_fileB)
  
  return(df)
}

build.refs<-function(input.ecopcr.results,output){
  #select one hit per genus
  a<-input.ecopcr.results[!duplicated(input.ecopcr.results$genus),]
  #concatenate the primer forward binding site, sequence, reverse binding site
  a$reverse_matchRC<-insect::rc(z = a$reverse_match)
  a$newseq<-paste0(a$forward_match, a$sequence,a$reverse_matchRC)
  a$definition<-paste0(a$AC," genus=",a$genus_name,"; taxid=",a$genus,";")
  export<-a[,c("definition","newseq")]
  invisible(seqRFLP::dataframe2fas(export,file = output))
}

map2targets<-function(queries.to.map,refs,out){
  #use blast to align seqs to ref
  message("mapping sequences to reference")
  system2(command = "makeblastdb", args=c("-in", refs, "-dbtype", "nucl", "-parse_seqids","-out","refdb"),wait=T)
  
  system2(command = "blastn", args=c("-query", queries.to.map, "-task", "megablast","-db","refdb",
                                     "-outfmt",'"7 qseqid qlen qstart qend slen sstart send length pident qcovs sstrand"',
                                     "-num_threads", "16","-max_target_seqs", "3"),stdout=out,wait = T)
}


add.counts.to.biasfile<-function(ncbiTaxDir,download.fasta,after.minL.fasta,after.checks.fasta,first.ecopcr.hit.table,
                                 mapped.fasta,out_bias_file,long.ecopcr.file){
  
  #catted DLS
  message("counting families in cattedDLS")
  cattedDLS<-count.fams.in.fasta(download.fasta,ncbiTaxDir)
  
  #after minL
  message("counting families after minL")
  afterminL<-count.fams.in.fasta(after.minL.fasta,ncbiTaxDir)
  
  #after rm fams and checks
  message("counting families after checks")
  afterchecks<-count.fams.in.fasta(after.checks.fasta,ncbiTaxDir)
  
  #after first ecopcr
  message("counting families after 1st ecopcr")
  
  firstecopcr<-data.table::fread(first.ecopcr.hit.table,sep = "\t")
  firstecopcr<-firstecopcr[!firstecopcr$amplicon_length>max_length,]
  firstecopcr$taxids<-firstecopcr$taxid
  firstecopcr<-add.lineage.df(firstecopcr,ncbiTaxDir)
  firstecopcr$count<-1
  firstecopcr$path<-paste(firstecopcr$K,firstecopcr$P,firstecopcr$C,firstecopcr$O,firstecopcr$F,sep = ";")
  a<-aggregate(firstecopcr$count,by=list(firstecopcr$path),FUN=sum)
  colnames(a)<-c("Family","nseqs")
  firstecopcr<-a
  
  #after long ecopcr
  message("counting families in long ecopcr")
  
  long.ecopcr<-data.table::fread(long.ecopcr.file,data.table = F)
  long.ecopcr<-long.ecopcr[long.ecopcr$amplicon_length>max_length & long.ecopcr$amplicon_length<long_length,]
  long.ecopcr$taxids<-long.ecopcr$taxid
  long.ecopcr<-add.lineage.df(long.ecopcr,ncbiTaxDir)
  long.ecopcr$count<-1
  long.ecopcr$path<-paste(long.ecopcr$K,long.ecopcr$P,long.ecopcr$C,long.ecopcr$O,long.ecopcr$F,sep = ";")
  a<-aggregate(long.ecopcr$count,by=list(long.ecopcr$path),FUN=sum)
  colnames(a)<-c("Family","nseqs")
  long.ecopcr1<-a
  
  #mean length in long ecopcr
  a<-aggregate(long.ecopcr$amplicon_length,by=list(long.ecopcr$path),FUN=mean)
  colnames(a)<-c("Family","nseqs")
  long.ecopcr2<-a
  
  #after mapping back
  message("counting families after mapping back")
  aftermapping<-count.fams.in.fasta(mapped.fasta,ncbiTaxDir)
  
  #merge all
  message("merging with bias table")
  merged<-merge(cattedDLS,afterminL,by = "Family",all = T)
  colnames(merged)<-gsub("nseqs.x","downloaded",colnames(merged))
  colnames(merged)<-gsub("nseqs.y","after_min_length",colnames(merged))
  merged<-merge(merged,afterchecks,by = "Family",all = T)
  colnames(merged)<-gsub("nseqs","after_checks",colnames(merged))
  merged<-merge(merged,firstecopcr,by = "Family",all = T)
  colnames(merged)<-gsub("nseqs","first_ecopcr",colnames(merged))
  
  merged<-merge(merged,long.ecopcr1,by = "Family",all = T)
  colnames(merged)<-gsub("nseqs","first_ecopcr_but_above_max_length",colnames(merged))
  
  merged<-merge(merged,long.ecopcr2,by = "Family",all = T)
  colnames(merged)<-gsub("nseqs","mean_length_for_above_max_length",colnames(merged))
  merged$mean_length_for_above_max_length<-round(merged$mean_length_for_above_max_length,digits = 0)
  
  merged<-merge(merged,aftermapping,by = "Family",all = T)
  colnames(merged)<-gsub("nseqs","after_mapping_back",colnames(merged))
  
  sum(merged$downloaded,na.rm = T)
  sum(merged$after_min_length,na.rm = T)
  sum(merged$after_checks,na.rm = T)
  sum(merged$first_ecopcr,na.rm = T)
  sum(merged$after_mapping,na.rm = T)
  
  biastemp<-data.table::fread(out_bias_file,sep = "\t")

  mergedbias<-merge(merged,biastemp,by.x ="Family",by.y = "in.odb",all = T)
  mergedbias$nseqs.odb=NULL
  
  #fixing taxonomy columns
  splittaxonomy<-do.call(rbind,stringr::str_split(mergedbias$Family,";"))
  mergedbias$K<-splittaxonomy[,1]
  mergedbias$P<-splittaxonomy[,2]
  mergedbias$C<-splittaxonomy[,3]
  mergedbias$O<-splittaxonomy[,4]
  mergedbias$F<-splittaxonomy[,5]
  
  mergedbias$path<-mergedbias$Family
  mergedbias$Family=NULL
  
  mergedbias<-mergedbias %>% select(path,K,P,C,O,F,everything())
  
  return(mergedbias)
}

count.fams.in.fasta<-function(fasta,ncbiTaxDir){
  message("Reminders: headers must contain \"taxid=taxid;\"")
  tempfasta<-phylotools::read.fasta(fasta)
  tempfasta$taxids<-stringr::str_match(tempfasta$seq.name, "taxid=(.*?);")[,2]
  tempfasta<-add.lineage.df(tempfasta,ncbiTaxDir)
  tempfasta$count<-1
  tempfasta$path<-paste(tempfasta$K,tempfasta$P,tempfasta$C,tempfasta$O,tempfasta$F,sep = ";")
  a<-aggregate(tempfasta$count,by=list(tempfasta$path),FUN=sum)
  colnames(a)<-c("Family","nseqs")
  return(a)
}

#read (and write) a processing sheet 
google.make.exp.sheet.illumina<-function(outDir,sheeturls,experiment_id){
  master<-list()
  headers<-c("Primer_set","Primer_F","Primer_R","Min_length","Max_length","ss_sample_id","experiment_id")
  for(i in 1:length(sheeturls)){
    master[[i]]<-google.read.master.url(sheeturls[i])
    
    #remove duplicated columns
    if(length(grep("\\.\\.\\.",colnames(master[[i]]),value = T))>0) {
      message("Removing duplicated columns")
      a<-grep("\\.\\.\\.",colnames(master[[i]]),value = T)
      b<-do.call(rbind,stringr::str_split(a,"\\.\\.\\."))
      d<-b[duplicated(b[,1]),]
      e<-paste0(d[,1],"...",d[,2])
      master[[i]]<-master[[i]][,!colnames(master[[i]]) %in% e]
      colnames(master[[i]])<-gsub("\\.\\.\\..*","",colnames(master[[i]]))
    }
    
    if(length(headers)!=sum(headers %in% colnames(master[[i]]))){
      stop (c("one of the following headers missing: ", paste(headers,collapse = " ")))}
    
    master[[i]]<-master[[i]][,headers]
  }
  
  #make a processing sheet
  experimentsheet<-as.data.frame(data.table::rbindlist(master))
  experimentsheet<-experimentsheet[experimentsheet$experiment_id %in% experiment_id,]
  
  #write file
  write.table(experimentsheet,paste0(outDir,paste(experiment_id,collapse = "_"),"_experiment_sheet.txt"),sep = "\t",quote = F,row.names = F)
  message(paste("file saved as",paste0(outDir,paste(experiment_id,collapse = "_"),"_experiment_sheet.txt")))
  
  return(experimentsheet)
}

#workaround for getting filenames of fastas to blast
google.get.startingfastas<-function(outDir,sheeturls,experiment_id,usingobiuniq){
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
  if(usingobiuniq){
    for(i in 1:length(primer_combo.bcs)){
      print(paste0(experiment_id,"_",names(primer_combo.bcs[[i]][1]),".uniq.filtlen.wlen.obi.fasta"))
    }}
  
  if(usingobiuniq==F){
    for(i in 1:length(primer_combo.bcs)){
      print(paste0(experiment_id,"_",names(primer_combo.bcs[[i]][1]),".filtlen.wlen.obi.fasta"))
    }}
  
}

ecopcr2refs2<-function(ecopcrfile,outfile,bufferecopcr,min_length=NULL,max_length=NULL){
  message("Reminder: assumes -D option was used for ecopcr")
  #read results
  ecopcroutput<-data.table::fread(ecopcrfile,sep = "\t",data.table = F)
  #remove hits outside desired lengths
  if(!is.null(min_length)) {ecopcroutput<-ecopcroutput[!ecopcroutput$amplicon_length<min_length,]}
  if(!is.null(max_length)) {ecopcroutput<-ecopcroutput[!ecopcroutput$amplicon_length>max_length,]}
  #remove duplicates (i.e. pick one entry per AC, based on lowest mismatches)
  ecopcroutput$total_mismatches<-as.numeric(ecopcroutput$forward_mismatch)+as.numeric(ecopcroutput$reverse_mismatch)
  ecopcroutput <- ecopcroutput[order(ecopcroutput$AC,ecopcroutput$total_mismatches),]
  ecopcroutput<-ecopcroutput[!duplicated(ecopcroutput$AC),]
  
  #the edges removal in next step takes a long time, choosing one seq per genus 
  ecopcroutput <- ecopcroutput[order(ecopcroutput$family_name,ecopcroutput$amplicon_length,decreasing = T),]
  ecopcroutput <-ecopcroutput[!duplicated(ecopcroutput[,c("family")]),] 
  
  #remove seqs with edges less than buffer
  a<-strsplit(ecopcroutput$sequence,split = "")
  ecopcroutput$leftbuffer<-"0"
  ecopcroutput$rightbuffer<-"0"
  
  d<-list()
  #pb = txtProgressBar(min = 0, max = length(ecopcroutput$leftbuffer), initial = 0,style = 3) 
  for(i in 1:length(ecopcroutput$leftbuffer)){
    ecopcroutput$leftbuffer[i]<-match(FALSE,a[[i]] %in% letters)
    d[[i]]<-a[[i]][ecopcroutput$leftbuffer[i]:length(a[[i]])]
    ecopcroutput$rightbuffer[i]<-length(d[[i]])-as.numeric(match(TRUE,d[[i]] %in% letters))
    #setTxtProgressBar(pb,i)
  }
  #remove seqs with left buffer less than buffer
  ecopcroutput<-ecopcroutput[!ecopcroutput$leftbuffer<bufferecopcr,]
  #remove seqs with right buffer less than buffer
  ecopcroutput<-ecopcroutput[!ecopcroutput$rightbuffer<bufferecopcr,]
  
  #output as fasta
  colnames(ecopcroutput)<-gsub("sequence", "seq.text",colnames(ecopcroutput))
  ecopcroutput$seq.text<-toupper(ecopcroutput$seq.text)
  ecopcroutput$seq.name<-paste0(ecopcroutput$AC," taxid=",ecopcroutput$taxid,"; ",ecopcroutput$definition)
  phylotools::dat2fasta(ecopcroutput[,c("seq.name","seq.text")],outfile)
}

add.3pmms<-function(ecopcroutput,Pf,Pr){
  
  message("Adding f_mismatches_3prime: the no. of mismatches in the 3 prime half of the forward primer")
  message("Adding r_mismatches_3prime: the no. of mismatches in the 3 prime half of the reverse primer")
  
  #add 3' mismatches to ecopcroutput
  f_mismatch_table<-mismatch.table(ecopcroutput,Pf,"f")
  f_mismatches_3prime<-as.data.frame(rowSums(f_mismatch_table[,as.integer(nchar(Pf)/2):nchar(Pf)]))
  colnames(f_mismatches_3prime)<-"f_mismatches_3prime"
  r_mismatch_table<-mismatch.table(ecopcroutput, Pr,"r")
  r_mismatches_3prime<-as.data.frame(rowSums(r_mismatch_table[,as.integer(nchar(Pr)/2):nchar(Pr)]))
  colnames(r_mismatches_3prime)<-"r_mismatches_3prime"
  ecopcroutput<-cbind(ecopcroutput,f_mismatches_3prime,r_mismatches_3prime)
  
  message("Adding f_mismatches_3prime6: the no. of mismatches in the 3 prime 6bp of the forward primer")
  message("Adding r_mismatches_3prime6: the no. of mismatches in the 3 primer 6bp of the reverse primer")
  
  #add 3' mismatches to ecopcroutput - 6bp
  f_mismatches_3prime6<-as.data.frame(rowSums(f_mismatch_table[,as.integer(nchar(Pf)-5):nchar(Pf)]))
  colnames(f_mismatches_3prime6)<-"f_mismatches_3prime6"
  r_mismatches_3prime6<-as.data.frame(rowSums(r_mismatch_table[,as.integer(nchar(Pr)-5):nchar(Pr)]))
  colnames(r_mismatches_3prime6)<-"r_mismatches_3prime6"
  ecopcroutput<-cbind(ecopcroutput,f_mismatches_3prime6,r_mismatches_3prime6)
}

add.tm.ecopcroutput<-function(ecopcroutput){
  message("Calculating Tm for forward & reverse primer binding sites: 1) full binding site 2) 3 prime half 3) 3 primer last 6bp")
  ecopcroutput$fTms<-Tm.calc(ecopcroutput$forward_match)
  ecopcroutput$rTms<-Tm.calc(ecopcroutput$reverse_match)
  ecopcroutput$fTms3primehalf<-Tm.calc(substr(x = ecopcroutput$forward_match,
                                              start = as.integer(nchar(as.character(ecopcroutput$forward_match))/2+1),
                                              stop = nchar(as.character(ecopcroutput$forward_match))))
  ecopcroutput$rTms3primehalf<-Tm.calc(substr(x = ecopcroutput$reverse_match,
                                              start = as.integer(nchar(as.character(ecopcroutput$reverse_match))/2+1),
                                              stop = nchar(as.character(ecopcroutput$reverse_match))))
  ecopcroutput$fTms3prime6<-Tm.calc(substr(x = ecopcroutput$forward_match,
                                           start = as.integer(nchar(as.character(ecopcroutput$forward_match))-5),
                                           stop = nchar(as.character(ecopcroutput$forward_match))))
  ecopcroutput$rTms3prime6<-Tm.calc(substr(x = ecopcroutput$reverse_match,
                                           start = as.integer(nchar(as.character(ecopcroutput$reverse_match))-5),
                                           stop = nchar(as.character(ecopcroutput$reverse_match))))
  ecopcroutput
  
}

add.gc.ecopcroutput<-function(ecopcroutput){
  message("Calculating % gc for each primer binding site")
  ecopcroutput$fgc<-gc.calc(ecopcroutput$forward_match)
  ecopcroutput$rgc<-gc.calc(ecopcroutput$reverse_match)
  ecopcroutput
}

#calculating Tm
#calculate Tm for each primer binding site in vector
Tm.calc<-function(primerVec){
  Tms<-vector()
  for (i in 1:length(primerVec)){
    Tms[i]<-2*(stringr::str_count(primerVec[i],"A")+stringr::str_count(primerVec[i],"T")) +
      4*(stringr::str_count(primerVec[i],"G")+stringr::str_count(primerVec[i],"C"))
  }
  return(Tms)
}

#calculating gc%
#calculate gc% for each primer binding site in vector
gc.calc<-function(primerVec){
  gc<-vector()
  for (i in 1:length(primerVec)){
    gc[i]<-(stringr::str_count(primerVec[i],"G")+stringr::str_count(primerVec[i],"C"))/
      (nchar(primerVec[i]))*100
  }
  return(gc)
}

#add taxonomic resolution to ecopcroutput
add.res.Bas<-function(ecopcroutput,obitaxdb){
  message("Adding resolution using Robitools method. Note that it fails when rank=no rank")
  
  if (class(obitaxdb)[1]!="obitools.taxonomy") obitaxdb2=ROBITaxonomy::read.taxonomy(obitaxdb)
  if (class(obitaxdb)[1]=="obitools.taxonomy") obitaxdb2=obitaxdb
  res=ROBIBarcodes::resolution(taxonomy = obitaxdb2,ecopcr = ecopcroutput)
  a<-as.data.frame(res)
  a$res<-gsub(x = a$res,pattern = "infraclass",replacement = "class")
  a$res<-gsub(x = a$res,pattern = "subfamily",replacement = "family")
  a$res<-gsub(x = a$res,pattern = "tribe",replacement = "family")
  a$res<-gsub(x = a$res,pattern = "suborder",replacement = "order")
  a$res<-gsub(x = a$res,pattern = "infraorder",replacement = "order")
  a$res<-gsub(x = a$res,pattern = "superfamily",replacement = "order")
  a$res<-gsub(x = a$res,pattern = "subclass",replacement = "class")
  a$res<-gsub(x = a$res,pattern = "subgenus",replacement = "genus")
  a$res<-gsub(x = a$res,pattern = "species subgroup",replacement = "species")
  a$res<-gsub(x = a$res,pattern = "cohort",replacement = "order")
  
  #K,P,C,O,F,G,S
  
  b<-cbind(ecopcroutput,a)
}

#function to calculate % species with family or lower res
res.fam.or.better.basres<-function(ecopcroutput,column="resolution.bas"){
  (length(ecopcroutput[,column][ecopcroutput[,column]=="family"])+
     length(ecopcroutput[,column][ecopcroutput[,column]=="genus"])+
     length(ecopcroutput[,column][ecopcroutput[,column]=="species"])) /
    (length(ecopcroutput[,column]))
}

res.fam.or.better.obi<-function(ecopcroutput,column="res"){
  (length(ecopcroutput[,column][ecopcroutput[,column]=="family"])+
     length(ecopcroutput[,column][ecopcroutput[,column]=="genus"])+
     length(ecopcroutput[,column][ecopcroutput[,column]=="species"])) /
    (length(ecopcroutput[,column]))
}



# #gc clamp present?
# gc.clamp<-function(primerVec){
#   gc<-vector()
#   for (i in 1:length(primerVec)){
#     gc[i]<-(stringr::str_count(primerVec[i],"G")+stringr::str_count(primerVec[i],"C"))/
#       (nchar(primerVec[i]))*100
#   }
#   return(gc)
# }


#Mean No. 3 prime PRIMER MISMATCHES (last half and last 6 bases) BY FAMILY
calc.3pmms.fam<-function(ecopcroutput,Pf,Pr){
  #mean primer mismatch tables for each base by family
  fam_f_mismatch_table<-family.mean.mismatch(ecopcroutput,Pf,"f")
  fam_r_mismatch_table<-family.mean.mismatch(ecopcroutput,Pr,"r")
  
  #mean 3' mismatches (last half of bases)
  fam_f_mismatches_3prime<-as.data.frame(rowMeans(fam_f_mismatch_table[,as.integer(nchar(Pf)/2+1):nchar(Pf)]))
  fam_f_mismatches_3prime$family<-rownames(fam_f_mismatches_3prime)
  colnames(fam_f_mismatches_3prime)<-c("mean_3prime_mms_f","family")
  fam_r_mismatches_3prime<-as.data.frame(rowMeans(fam_r_mismatch_table[,as.integer(nchar(Pr)/2+1):nchar(Pr)]))
  fam_r_mismatches_3prime$family<-rownames(fam_r_mismatches_3prime)
  colnames(fam_r_mismatches_3prime)<-c("mean_3prime_mms_r","family")
  
  #mean 3' mismatches (last 6 bases)
  fam_f_mismatches_3prime6<-as.data.frame(rowMeans(fam_f_mismatch_table[,as.integer(nchar(Pf)-5):nchar(Pf)]))
  fam_f_mismatches_3prime6$family<-rownames(fam_f_mismatches_3prime6)
  colnames(fam_f_mismatches_3prime6)<-c("mean_3prime6bp_mms_f","family")
  fam_r_mismatches_3prime6<-as.data.frame(rowMeans(fam_r_mismatch_table[,as.integer(nchar(Pr)-5):nchar(Pr)]))
  fam_r_mismatches_3prime6$family<-rownames(fam_r_mismatches_3prime6)
  colnames(fam_r_mismatches_3prime6)<-c("mean_3prime6bp_mms_r","family")
  
  list(fam_f_mismatches_3prime,fam_r_mismatches_3prime,fam_f_mismatches_3prime6,fam_r_mismatches_3prime6)
}

calc.fam.res<-function(ecopcroutput,res="bas2"){
  #calculate percentage species that have tax res to family or better
  famsplit<-split(ecopcroutput,f = ecopcroutput$path)
  if(res=="bas2") {
    b<-lapply(X = famsplit,FUN = res.fam.or.better.basres)
  } else {
    b<-lapply(X = famsplit,FUN = res.fam.or.better.obi)
  }
  
  pc.res<-as.data.frame(t(as.data.frame(b)))
  pc.res$family<-names(famsplit)
  colnames(pc.res)<-gsub("V1","pc_res_fam_or_better",colnames(pc.res))
  pc.res
}

calc.stat.ecopcroutput<-function(ecopcroutput,variable,appliedstat="mean"){
  if(appliedstat=="mean"){
    a<-as.data.frame(aggregate(x =  ecopcroutput[,variable], by = list(ecopcroutput$path),FUN = mean))
    colnames(a)<-c("path",paste0("mean_",variable))
    }
  if(appliedstat=="var"){
    a<-as.data.frame(aggregate(x =  ecopcroutput[,variable], by = list(ecopcroutput$path),FUN = var))
    colnames(a)<-c("path",paste0("var_",variable))
    }
  a
}


modify.ecopcroutput<-function(ecopcrfile,out,min_length=NULL,max_length=NULL){
  #GENERAL CLEANING 
  
  #read results
  ecopcroutput<-data.table::fread(ecopcrfile,sep = "\t")
  message("Removing hits outside desired lengths, if provided")
  if(!is.null(min_length)) ecopcroutput<-ecopcroutput[!ecopcroutput$amplicon_length<min_length,]
  if(!is.null(max_length)) ecopcroutput<-ecopcroutput[!ecopcroutput$amplicon_length>max_length,]
  message("Removing duplicates (i.e. pick one entry per AC, based on lowest mismatches)")
  ecopcroutput$total_mismatches<-as.numeric(ecopcroutput$forward_mismatch)+as.numeric(ecopcroutput$reverse_mismatch)
  ecopcroutput <- ecopcroutput[order(ecopcroutput$AC,ecopcroutput$total_mismatches),]
  ecopcroutput<-ecopcroutput[!duplicated(ecopcroutput$AC),]
  message("Removing weird primer mismatches (only a few usually)")
  ecopcroutput<-ecopcroutput[!nchar(ecopcroutput$forward_match,allowNA = T)<nchar(Pf),]
  ecopcroutput<-ecopcroutput[!nchar(ecopcroutput$reverse_match,allowNA = T)<nchar(Pr),]
  
  #write  file
  write.table(x=ecopcroutput,file = out,quote = F,sep = "\t",row.names = F)
  
}

#DO CALCULATIONS AND STATS
make.primer.bias.table<-function(originaldbtab,
                                  out_bias_file,mod_ecopcrout_file,Pf,Pr, obitaxoR,min_length,
                                  max_length){
  
  ecopcroutput<-data.table::fread(mod_ecopcrout_file,sep = "\t",data.table = F)
  
  ##########################################
  #READ ECOPCRDB and add lineages
  originaldb<-as.data.frame(data.table::fread(originaldbtab,header = TRUE,sep = "\t"))
  colnames(originaldb)<-gsub("taxid","taxids",colnames(originaldb))
  originaldb<-add.lineage.df(originaldb,ncbiTaxDir)
  ##########################################
  #LIST FAMILIES IN ODB, use full path as some families have same names
  originaldb$path<-paste(originaldb$K,originaldb$P,originaldb$C,originaldb$O,originaldb$F,sep = ";")
  all_primer_bias<-data.frame(row.names = 1:length(unique(originaldb$path)))
  all_primer_bias$in.odb<-unique(originaldb$path)[order(as.character(unique(originaldb$path)))]
  ##########################################
  #COUNT NO. SEQS IN ORIGINALDB 
  originaldb$n<-1
  nseqs.odb<-aggregate(originaldb$n,by=list(originaldb$path),FUN = sum)
  colnames(nseqs.odb)<-c("in.odb","nseqs.odb")
  all_primer_bias<-merge(all_primer_bias,nseqs.odb,by = "in.odb",all.x = T)
  ##########################################
  #COUNT No. Taxa IN ORIGINALDB 
  originaldb$path<-paste(originaldb$K,originaldb$P,originaldb$C,originaldb$O,originaldb$F,originaldb$G,originaldb$S,sep = ";")
  originaldb.taxa<-originaldb[!duplicated(originaldb$path),c("path","n")]
  aa<-do.call(rbind,stringr::str_split(originaldb.taxa$path,";"))
  originaldb.taxa$path<-paste(aa[,1],aa[,2],aa[,3],aa[,4],aa[,5],sep = ";")
  nseqs.odb<-aggregate(originaldb.taxa$n,by=list(originaldb.taxa$path),FUN = sum)
  colnames(nseqs.odb)<-c("in.odb","ntaxa.odb")
  all_primer_bias<-merge(all_primer_bias,nseqs.odb,by = "in.odb",all.x = T)
  ##########################################
  #ADD LOGICAL IF FAMILY AMPED
  ecopcroutput$path<-paste(ecopcroutput$K,ecopcroutput$P,ecopcroutput$C,ecopcroutput$O,ecopcroutput$F,sep = ";")
  all_primer_bias$amplified<-all_primer_bias$in.odb %in% unique(ecopcroutput$path) 
  ##########################################
  #COUNT NO. SEQS THAT AMPED 
  ecopcroutput$n<-1
  nseqs.amped<-aggregate(ecopcroutput$n,by=list(ecopcroutput$path),FUN = sum)
  colnames(nseqs.amped)<-c("in.odb","nseqs.amped")
  all_primer_bias<-merge(all_primer_bias,nseqs.amped,by = "in.odb",all.x = T)
  ##########################################
  #COUNT No. Taxa THAT AMPED
  ecopcroutput$path<-paste(ecopcroutput$K,ecopcroutput$P,ecopcroutput$C,ecopcroutput$O,ecopcroutput$F,ecopcroutput$G,ecopcroutput$S,sep = ";")
  ecopcroutput.taxa<-ecopcroutput[!duplicated(ecopcroutput$path),c("path","n")]
  aa<-do.call(rbind,stringr::str_split(ecopcroutput.taxa$path,";"))
  ecopcroutput.taxa$path<-paste(aa[,1],aa[,2],aa[,3],aa[,4],aa[,5],sep = ";")
  nseqs.odb<-aggregate(ecopcroutput.taxa$n,by=list(ecopcroutput.taxa$path),FUN = sum)
  colnames(nseqs.odb)<-c("in.odb","ntaxa.amped")
  all_primer_bias<-merge(all_primer_bias,nseqs.odb,by = "in.odb",all.x = T)
  
  ##########################################
  #percentage seqs amped
  all_primer_bias$pc_seqs_amped<-round(all_primer_bias$nseqs.amped/all_primer_bias$nseqs.odb*100,digits = 2)
  #percentage taxa amped
  all_primer_bias$pc_taxa_amped<-round(all_primer_bias$ntaxa.amped/all_primer_bias$ntaxa.odb*100,digits = 2)
  ##########################################
  
  #for rest of stats keep only unique barcodes for each family
  message("Dereplicating by family. Bias tables are based on dereplicated seqeunces only")
  ecopcroutput$path<-paste(ecopcroutput$K,ecopcroutput$P,ecopcroutput$C,ecopcroutput$O,ecopcroutput$F,sep = ";")
  ecopcroutput<-ecopcroutput[!duplicated(ecopcroutput[,c("path","fullseq")]),]
  
  #COUNT NO. UNIQUE BARCODES THAT AMPED 
  nseqs.amped<-aggregate(ecopcroutput$n,by=list(ecopcroutput$path),FUN = sum)
  colnames(nseqs.amped)<-c("in.odb","n.uniq.brcds.amped")
  all_primer_bias<-merge(all_primer_bias,nseqs.amped,by = "in.odb",all.x = T)
  
  
  ###########################################################################################
  #CALCULATE MEANS
  
  #mean no. primer mismatches fw
  mean_mms_fw<-calc.stat.ecopcroutput(ecopcroutput,variable="forward_mismatch","mean")
  #mean no. primer mismatches rv
  mean_mms_rv<-calc.stat.ecopcroutput(ecopcroutput,variable="reverse_mismatch","mean")
  #mean no. primer mismatches total
  mean_mms_total<-calc.stat.ecopcroutput(ecopcroutput,variable="total_mismatches","mean")
  #Mean no. 3 prime mismatches (last half and last 6 bases) 
  mean.3pmms<-calc.3pmms.fam(ecopcroutput,Pf,Pr) ######LIST
  mean_fmms3Phalf<-mean.3pmms[[1]]
  mean_rmms3Phalf<-mean.3pmms[[2]]
  mean_fmms3P6<-mean.3pmms[[3]]
  mean_rmms3P6<-mean.3pmms[[4]]
  #mean ftm 
  mean_ftm<-calc.stat.ecopcroutput(ecopcroutput,variable="fTms","mean")
  #mean rtm 
  mean_rtm<-calc.stat.ecopcroutput(ecopcroutput,variable="rTms","mean")
  #mean fTm 3' half
  mean_ftm3Phalf<-calc.stat.ecopcroutput(ecopcroutput,variable="fTms3primehalf","mean")
  #mean rTm 3' half
  mean_rtm3Phalf<-calc.stat.ecopcroutput(ecopcroutput,variable="rTms3primehalf","mean")
  #mean fTm for 3' half - 6bp  
  mean_ftm3P6<-calc.stat.ecopcroutput(ecopcroutput,variable="fTms3prime6","mean")  
  #mean rTm for 3' half - 6bp  
  mean_rtm3P6<-calc.stat.ecopcroutput(ecopcroutput,variable="rTms3prime6","mean")  
  #% of unqiue seqs that get to fam or better
  mean_fam.res.bas<-calc.fam.res(ecopcroutput,res = "bas2")
  colnames(mean_fam.res.bas)<-gsub("pc_res_fam_or_better", "res.fam.plus.bas",colnames(mean_fam.res.bas))
  #% of unqiue seqs that get to fam or better
  mean_fam.res.obi<-calc.fam.res(ecopcroutput,res = "obi")
  colnames(mean_fam.res.obi)<-gsub("pc_res_fam_or_better", "res.fam.plus.obi",colnames(mean_fam.res.obi))
  #mean amplicon length
  mean_amplicon.len<-aggregate(x =  nchar(ecopcroutput$sequence), by = list(ecopcroutput$path),FUN = mean)
  colnames(mean_amplicon.len)<-gsub("x","mean_amplicon.len",colnames(mean_amplicon.len))
  #mean gc_fw
  mean_fgc<-calc.stat.ecopcroutput(ecopcroutput,variable="fgc","mean")
  #mean gc_rv
  mean_rgc<-calc.stat.ecopcroutput(ecopcroutput,variable="rgc","mean")
  #mean tm.ta.fw
  mean_tm.ta.fw<-calc.stat.ecopcroutput(ecopcroutput,variable="tm.ta_fw","mean")
  #mean tm.ta.rv
  mean_tm.ta.rv<-calc.stat.ecopcroutput(ecopcroutput,variable="tm.ta_rv","mean")
  #mean diff tm
  mean_diff_tm<-calc.stat.ecopcroutput(ecopcroutput,variable="diff_tm","mean")
  #mean tm_fw_3p6_perc
  mean_tm_fw_3p6_perc<-calc.stat.ecopcroutput(ecopcroutput,variable="tm_fw_3p6_perc","mean")
  #mean tm_rv_3p6_perc
  mean_tm_rv_3p6_perc<-calc.stat.ecopcroutput(ecopcroutput,variable="tm_rv_3p6_perc","mean")
  ##########################################################################################
  #VARIANCES
  #var no. primer mismatches fw
  var_mms_fw<-calc.stat.ecopcroutput(ecopcroutput,variable="forward_mismatch","var")
  #var no. primer mismatches rv
  var_mms_rv<-calc.stat.ecopcroutput(ecopcroutput,variable="reverse_mismatch","var")
  #var no. primer mismatches total
  var_mms_total<-calc.stat.ecopcroutput(ecopcroutput,variable="total_mismatches","var")
  #var amplicon length
  var_amplicon.len<-aggregate(x =  nchar(ecopcroutput$sequence), by = list(ecopcroutput$path),FUN = var)
  colnames(var_amplicon.len)<-gsub("x","var_amplicon.len",colnames(var_amplicon.len))
  #var gc_fw
  var_fgc<-calc.stat.ecopcroutput(ecopcroutput,variable="fgc","var")
  #var gc_rv
  var_rgc<-calc.stat.ecopcroutput(ecopcroutput,variable="rgc","var")
  #var tm.ta_fw
  var_tm.ta.fw<-calc.stat.ecopcroutput(ecopcroutput,variable="tm.ta_fw","var")
  #var tm.ta_rv
  var_tm.ta.rv<-calc.stat.ecopcroutput(ecopcroutput,variable="tm.ta_rv","var")
  #diff tm
  var_diff_tm<-calc.stat.ecopcroutput(ecopcroutput,variable="diff_tm","var")
  #var tm_fw_3p6_perc
  var_tm_fw_3p6_perc<-calc.stat.ecopcroutput(ecopcroutput,variable="tm_fw_3p6_perc","var")
  #var tm_rv_3p6_perc
  var_tm_rv_3p6_perc<-calc.stat.ecopcroutput(ecopcroutput,variable="tm_rv_3p6_perc","var")
  #mms_fw_3p6
  var_mms_fw_3p6<-calc.stat.ecopcroutput(ecopcroutput,variable="f_mismatches_3prime6","var")
  #mms_rv_3p6
  var_mms_rv_3p6<-calc.stat.ecopcroutput(ecopcroutput,variable="r_mismatches_3prime6","var")
  
  ##########################################################################################
  #compile
  h<-list()
  for(i in 1:length(ls(pattern = "mean\\_.*"))){
    h[[i]]<-get(ls(pattern = "mean\\_.*")[i])
  }
  all_means<-do.call(cbind,h)
  h<-list()
  for(i in 1:length(ls(pattern = "var\\_.*"))){
    h[[i]]<-get(ls(pattern = "var\\_.*")[i])
  }
  all_vars<-do.call(cbind,h)
  
  all_mean_and_var<-cbind(all_means,all_vars)
  all_mean_and_var<-all_mean_and_var[,-grep("family",colnames(all_mean_and_var))]
  all_mean_and_var$Group.1.1=NULL
  
  #compile all
  all_primer_bias<-merge(all_primer_bias,all_mean_and_var,all.x = T,by.x = "in.odb", by.y = "Group.1")
  
  #remove extraneuos
  all_primer_bias<-all_primer_bias[,-grep("path",colnames(all_primer_bias))]
  
  #write primer bias file
  write.table(x=all_primer_bias,file = out_bias_file,quote = F,sep = "\t",row.names = F)
  
}

mapTrim2.simple<-function(query,blast.results.file,qc=0.7,out,pident=20){
  #read in query file and count
  n<-phylotools::read.fasta(query)
  n$qseqid = sub(" .*", "", x = n$seq.name)
  n$seq.text<-as.character(n$seq.text)
  #count
  count_queries=length(n$qseqid)
  
  #turn result into table, ignoring comment lines (so ignoring defintions). The table length excludes queries with
  #non-hits
  j<-read.table(file = blast.results.file)
  colnames(j)<-c("qseqid", "qlen", "qstart", "qend","slen", "sstart", "send", "length", "pident", "qcovs","sstrand")
                 
  #calculate scov
  j$scov<-j$length/j$slen
  
  #remove duplicate hsp hits, based on highest scov
  j2 <- j[order(j$qseqid,j$scov,decreasing = T),]
  j2<-j2[!duplicated(j2$qseqid),]
  count_hits<-length(j2$qseqid)
  
  #remove hits less than specified scov
  j2<-j2[j2$scov>qc,]
  count_qc<-length(j2$qseqid)
  
  #remove hits less than specified pident
  j2<-j2[j2$pident>pident,]
  count_pident<-length(j2$pident)
  
  #merge query and blast results
  k<-merge(x = j2,y = n,by = "qseqid",all.y = F) 
  
  message("outputting as fasta")
  message("testing keeping all sequence that passed buffers, rather than extracting sequence")
  k_export<-k[,c("seq.name","seq.text")]
  count_final_db<-length(k_export$seq.name)
  phylotools::dat2fasta(k_export,outfile = out)
  
  ###########################
  message(paste("Done.", "From", count_queries,"sequences,", count_hits, "mapped to a reference,",count_qc,
            "of which had >",qc*100,"% coverage,",count_pident,"of which had greater than",pident,"% percentage identity"))
}



#####################################################################################


mapTrim2<-function(query,buffer,blast.results.file,qc=0.7,out){
  message("reading query file")
  #read in query file and count
  n<-phylotools::read.fasta(query)
  n$qseqid = sub(" .*", "", x = n$seq.name)
  n$seq.text<-as.character(n$seq.text)
  #count
  count_queries=length(n$qseqid)
  
  message("reading mapping results")
  #turn result into table, ignoring comment lines (so ignoring defintions). The table length excludes queries with
  #non-hits
  j<-read.table(file = blast.results.file)
  colnames(j)<-c("qseqid", "qlen", "qstart", "qend",
                 "slen", "sstart", "send", "length", "pident", "qcovs","sstrand")
  #calculate scov
  j$scov<-j$length/j$slen
  
  #remove duplicate hits
  j2 <- j[order(j$qseqid,j$scov),]
  j2<-j2[!duplicated(j2$qseqid),]
  count_hits<-length(j2$qseqid)
  
  #remove hits less than specified scov
  j2<-j2[j2$scov>qc,]
  count_qc<-length(j2$qseqid)
  
  #merge query and blast results
  k<-merge(x = j2,y = n,by = "qseqid",all.y = F)
  
  #find query position that matches first subject position (i.e. first base of primer-binding site,
  #or first base of buffer if ref already included buffer)
  ##the strand problem...
  kplus<-k[k$sstrand=="plus",]
  kminus<-k[k$sstrand=="minus",]
  kplus$q_start_base<-kplus$qstart-kplus$sstart
  kminus$q_start_base<-kminus$qstart-kminus$send
  kboth<-rbind(kplus,kminus)
  
  #subtract the buffer from this (to only accept queries that extend left of primer)
  kboth$q_left_buff<-kboth$q_start_base-buffer
  #find final query base required
  kboth$q_right_buff<-kboth$q_left_buff+kboth$slen+buffer*2
  
  #remove queries with negative starting base
  k2<-kboth[kboth$q_left_buff>-1,]
  #count
  count_left_buffer<-length(k2$qseqid)
  
  #extract the sequence within thresholds
  k2$ex.seq<-substr(k2$seq.text,start = k2$q_left_buff,stop = k2$q_right_buff)
  #remove queries shorter than expected
  k2$newlen<-nchar(k2$ex.seq)
  #caLculate desired length, with slight leeway (to account for 1-3bp difference in calculation)
  k2$desiredlen<-k2$q_right_buff-k2$q_left_buff-3
  k3<-k2[k2$newlen>k2$desiredlen,]
  #count
  count_right_buffer<-length(k3$qseqid)
  
  message("outputting as fasta")
  message("testing keeping all sequence that passed buffers, rather than extracting sequence")
  k3_export<-k3[,c("seq.name","seq.text")]
  #k3_export<-k3[,c("seq.name","ex.seq")]
  #colnames(k3_export)<-c("seq.name","seq.text")
  count_final_db<-length(k3_export$seq.name)
  phylotools::dat2fasta(k3_export,outfile = out)
  
  ###########################
  message(c("Done. ", "From ", count_queries," sequences, ", count_hits, " mapped to a reference, ",count_qc,
            " of which had > ",qc*100,"% coverage, ",
            count_left_buffer," of which passed left buffer and, of those, ",count_right_buffer," passed right buffer."))
}

#' Convert BLAST into rma6. See megan help file (\code{blast2rma.BAS(h=T)}) for full details.
#'     This is just an R wrapper of that function. This does not support all options for the \code{blast2rma} command.
#' @title Convert BLAST into rma6
#' @param infile BLAST results file
#' @param format format of BLAST results file. Common values: "BlastTab" (default),"BlastXML", "BlastText". See h=T for more
#' @param blastMode mode used for BLAST
#' @param reads Optional. fasta file used to perform BLAST
#' @param outfile filename for output file. Recommended to use .rma6 suffix
#' @param top Optional. When binning, consider hits falling within \code{top} percent of BLAST score of the top hit
#' @param mdf Optional. Files containing metadata to be included in RMA6 files ######Need to find the accpeted formats for this
#' @param ms Min score. Default value: 50.0
#' @param me Max expected. Default value: 0.01
#' @param mrc Min percent of read length to be covered by alignments. Default value: 70.0
#' @param ram Set the read assignment mode. Default value: readCount. Legal values: readCount, readLength, alignedBases, readMagnitude
#' @param a2t Optional, bit highly recommended. Accession-to-Taxonomy mapping file
#' @param h Set h=T to show program usage
#' @return An rma6 file
#' @note Apparently, if \code{infile} is xml format, filename must end in ".xml"
#'
#' @examples
#' blast2rma.BAS("primer_16S.uniq.l75L120.c20.xml",outfile = "primer_16S.uniq.l75L120.c20.TEST.blast2rma.rma6",
#'    a2t = "nucl_acc2tax-Nov2018.abin")
#' #inspect file and disable taxa as necessary
#' primer_16S.uniq.l75L120.c20.TEST.taxon.table<-rma2info.BAS("primer_16S.uniq.l75L120.c20.TEST.blast2rma.rma6")
#' final.table<-merge.tab.taxon(obitab.txt = "primer_16S.uniq.l75L120.c20.tab",primer_16S.uniq.l75L120.c20.TEST.taxon.table)
#' @export
blast2rma.BAS<-function(infile,format="BlastTab",blastMode="BlastN",outfile,reads=F,
                        top=10.0, mdf=F,ms=50.0,me=0.01,mrc=70.0, ram="readCount",a2t=F,h=F){
  if(reads==F) message("Running without a specified fasta file, which is fine, but resultant megan file will not contain as much information")
  
  if(h==T){
    processx::run(command = "blast2rma", args="-h",echo = T)}
  
  if(h==F){
    cb <- function(line, proc) {cat(line, "\n")}
    argsBas<-c(" --in ", infile," --format",format," --blastMode ",blastMode,
               " --reads ", reads," --out ",outfile," --topPercent ",top," -mdf ",mdf," -ms ",ms,
               " -me ",me," -mrc ",mrc," -ram ",ram," -a2t ",a2t," -supp ",0)
    argsBasF<-argsBas[-(grep(FALSE,argsBas)-1)]
    argsBasF<-argsBasF[-grep(FALSE,argsBasF)]
    processx::run(command = "blast2rma", args=argsBasF,stderr_line_callback = cb,echo_cmd = T,echo = F)}
  b<-"SUCCESS!" #this is mainly to stop processx::run printing stdout and stderr to screen, which it does because we
  #cannot explicitly redirect stdout because the command has an "out" option
  return(b)
}

#' Convert rma6 file into a table with read names, associated taxon path and a letter denoting the lowest taxonomic
#'     level reached. See megan help file (\code{rma2info.BAS(h=T)}) for full details.
#'     This is just an R wrapper of the \code{rma2info} command, it does not support all options for the \code{rma2info} command.
#' @title Convert rma6 into taxon path table.
#' @param infile rma6 file
#' @param names Report names rather than taxid numbers. Default value: true.
#' @param h Show program usage
#' @return A dataframe consisting of read names in the "id" column, a letter denoting the lowest taxonomic in the "resolution"
#'     column and seven further columns for the taxon path split by major taxonomic ranks (SK,P,C,O,F,G,S).
#'     The function is fixed to output the taxonomy path to lowest level reached.
#' @examples
#' a<-system.file("extdata", "Mblast.c20.uniq.l85L105.PRIMER_16S.rma6", package = "bastools")
#' b<-rma2info.BAS(a)
#' @export
rma2info.BAS<-function(infile,names="true",h=F,out.txt=NULL){
  if(h==T){
    processx::run(command = "rma2info", args="-h",echo = T)}
  
  if(h==F){
    cb <- function(line, proc) {cat(line, "\n")}
    if(is.null(out.txt)){
      argsBas<-c(" --in ", infile,
                 " -r2c ", "Taxonomy", " --names ",names," -u ","true"," -mro ","true",
                 " --ranks ","true"," --paths ","true"," -v ", "true")
      a<-processx::run(command = "rma2info", args=argsBas, echo_cmd = T,stderr_line_callback = cb)
      
      b<-read.delim(text = a$stdout,header=F)
      colnames(b)<-c("id","resolution","path")
      c<-as.data.frame(stringr::str_split(b$path,pattern = ";",simplify = T))
      d<-c[,1:(length(colnames(c))-1)]
      colnames(d)<-c("SK","P","C","O","F","G","S")
      d[] <- lapply(d, function(x) (gsub("\\[.*\\] ", "", x)))
      e<-cbind(b[,c("id","resolution")],d)
      return(e)
    }
    
    if(!is.null(out.txt)){
      argsBas<-c(" --in ", infile,
                 " -r2c ", "Taxonomy", " --names ",names," -u ","true"," -mro ","true",
                 " --ranks ","true"," --paths ","true"," -v ", "true")
      h<-process$new(command = "rma2info", args=argsBas, echo_cmd = T,stdout = out.txt)
      h$wait()
      b<-"SUCCESS!" #this is mainly to stop processx::run printing stdout and stderr to screen, which it does because we
      #cannot explicitly redirect stdout because the command has an "out" option
      return(b)
      #########this will not be in correct format for merge.tab.taxon
    }
  }
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
merge.tab.taxon<-function(obitab.txt,megan.taxa){
  
  if (class(megan.taxa)!="data.frame") {
    read.csv(file = obitab.txt, sep = "\t")->obitab_input
    read.csv(file = megan.taxa, sep = "\t", header = FALSE)->taxon_input
    taxon_input$taxon<-taxon_input$V2
    taxon_input$V2<-NULL
    merged.table<-merge(obitab_input,taxon_input,by.x = "id",by.y = "V1",all.x = TRUE)
  }
  
  if (class(megan.taxa)=="data.frame") {
    read.csv(file = obitab.txt, sep = "\t")->obitab_input
    merged.table<-merge(obitab_input,megan.taxa,by.x = "id",by.y = "id",all.x = TRUE)
  }
  return(merged.table)
}

family.mean.mismatch<-function(ecopcroutput_clean,primer.seq=NULL,primer.direction){
  
  #split query forward primer
  if(primer.direction=="f"){
    sst <- as.data.frame(t(as.data.frame(strsplit(as.character(ecopcroutput_clean$forward_match), ""))))}
  if(primer.direction=="r"){
    sst <- as.data.frame(t(as.data.frame(strsplit(as.character(ecopcroutput_clean$reverse_match), ""))))}
  sst$path<-ecopcroutput_clean$path
  rownames(sst)<-NULL
  
  #split subject primer into list
  pflist<-list()
  for(i in 1:nchar(primer.seq)){
    pflist[[i]]<-substr(primer.seq,start = i,stop = i)
  }
  
  #IUPAC RULES IN LIST
  IUPAC=list(R=c("A","G"),
             Y=c("C","T"),
             S=c("G","C"),
             W=c("A","T"),
             K=c("G","T"),
             M=c("A","C"),
             B=c("C","G","T"),
             D=c("A","G","T"),
             H=c("A","C","T"),
             V=c("A","C","G"),
             N=c("C","G","T","A"))
  
  #replace subject primer baseswith IUPAC ambiguities
  suppressWarnings(for(i in 1:length(pflist)){
    if(pflist[[i]]=="Y") pflist[[i]]<-IUPAC[["Y"]]
    if(pflist[[i]]=="S") pflist[[i]]<-IUPAC[["S"]]
    if(pflist[[i]]=="W") pflist[[i]]<-IUPAC[["W"]]
    if(pflist[[i]]=="K") pflist[[i]]<-IUPAC[["K"]]
    if(pflist[[i]]=="M") pflist[[i]]<-IUPAC[["M"]]
    if(pflist[[i]]=="B") pflist[[i]]<-IUPAC[["B"]]
    if(pflist[[i]]=="D") pflist[[i]]<-IUPAC[["D"]]
    if(pflist[[i]]=="H") pflist[[i]]<-IUPAC[["H"]]
    if(pflist[[i]]=="V") pflist[[i]]<-IUPAC[["V"]]
    if(pflist[[i]]=="N") pflist[[i]]<-IUPAC[["N"]]
  })
  
  #split ecopcroutput_clean by family
  splitdf<-split(sst, f = sst$path)
  
  #check if query primer is in subject primer
  #function to apply %in% to one data.frame
  matchdf<-function(df,pflist){
    out<-data.frame(matrix(ncol=length(pflist),nrow = length(df[,1])))
    for(i in 1:length(pflist)){
      out[,i]<-df[,i] %in% pflist[[i]]
    }
    return(out)
  }
  
  #apply to all dataframe in list
  out2=list()
  for(i in 1:length(splitdf)){
    out2[[i]]<-matchdf(splitdf[[i]],pflist)
  }
  
  #convert true/false to 1/0
  out3=list()
  for(i in 1:length(out2)){
    out3[[i]]<-out2[[i]]*1
  }
  
  #swap 0s for 1s
  for(i in 1:length(out3)){
    out3[[i]][out3[[i]]==1]<-2
    out3[[i]][out3[[i]]==0]<-1
    out3[[i]][out3[[i]]==2]<-0
  }
  
  #calculate mean mismatches for each family for each base
  out4<-lapply(X=out3,FUN = colMeans)
  names(out4)<-names(splitdf)
  out5<-as.data.frame(t(bind_rows(out4)))
  colnames(out5)<-as.character(strsplit(x = primer.seq,"")[[1]])
  
  return(out5)
}




mismatch.table<-function(ecopcroutput_clean,primer.seq=NULL,primer.direction){
  
  #split query forward primer
  if(primer.direction=="f"){
    sst <- as.data.frame(t(as.data.frame(strsplit(as.character(ecopcroutput_clean$forward_match), ""))))}
  if(primer.direction=="r"){
    sst <- as.data.frame(t(as.data.frame(strsplit(as.character(ecopcroutput_clean$reverse_match), ""))))}
  sst$family_name<-ecopcroutput_clean$family_name
  rownames(sst)<-NULL
  
  #split subject primer into list
  pflist<-list()
  for(i in 1:nchar(primer.seq)){
    pflist[[i]]<-substr(primer.seq,start = i,stop = i)
  }
  
  #IUPAC RULES IN LIST
  IUPAC=list(R=c("A","G"),
             Y=c("C","T"),
             S=c("G","C"),
             W=c("A","T"),
             K=c("G","T"),
             M=c("A","C"),
             B=c("C","G","T"),
             D=c("A","G","T"),
             H=c("A","C","T"),
             V=c("A","C","G"),
             N=c("C","G","T","A"))
  
  #replace subject primer baseswith IUPAC ambiguities
  suppressWarnings(for(i in 1:length(pflist)){
    if(pflist[[i]]=="R") pflist[[i]]<-IUPAC[["R"]]
    if(pflist[[i]]=="Y") pflist[[i]]<-IUPAC[["Y"]]
    if(pflist[[i]]=="S") pflist[[i]]<-IUPAC[["S"]]
    if(pflist[[i]]=="W") pflist[[i]]<-IUPAC[["W"]]
    if(pflist[[i]]=="K") pflist[[i]]<-IUPAC[["K"]]
    if(pflist[[i]]=="M") pflist[[i]]<-IUPAC[["M"]]
    if(pflist[[i]]=="B") pflist[[i]]<-IUPAC[["B"]]
    if(pflist[[i]]=="D") pflist[[i]]<-IUPAC[["D"]]
    if(pflist[[i]]=="H") pflist[[i]]<-IUPAC[["H"]]
    if(pflist[[i]]=="V") pflist[[i]]<-IUPAC[["V"]]
    if(pflist[[i]]=="N") pflist[[i]]<-IUPAC[["N"]]
  })
  
  #check if query primer is in subject primer
  #function to apply %in% to one data.frame
  matchdf<-function(df,pflist){
    out<-data.frame(matrix(ncol=length(pflist),nrow = length(df[,1])))
    for(i in 1:length(pflist)){
      out[,i]<-df[,i] %in% pflist[[i]]
    }
    return(out)
  }
  
  #apply
  out2<-matchdf(sst,pflist)
  
  #convert true/false to 1/0
  out3<-out2*1
  #swap 0s for 1s
  
  out3[out3==1]<-2
  out3[out3==0]<-1
  out3[out3==2]<-0
  
  return(out3)
}

#' get classic taxonomy table from taxids
#' @param x A dataframe with a column called "taxid".
#' @param taxo An obitools-formatted taxonomy database as an R object, generated by \code{\link[ROBITaxonomy]{read.taxonomy}}.
#' @return The original dataframe \code{x} with seven new taxonomy columns, kingdom to species.
#' @examples
#' get.classic.taxonomy.BAS(My_df_with_taxid_col, my_obitools-formatted_taxonomy_as_object)
#' @export
get.classic.taxonomy.Bas = function(x, obitaxdb) {
  classic.taxo = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  taxids = x$taxid
  out = as.data.frame(do.call("cbind", lapply(classic.taxo, function(y) {
    ROBITaxonomy::scientificname(obitaxdb, ROBITaxonomy::taxonatrank(obitaxdb,taxids,y))
  })))
  colnames(out) = paste(classic.taxo, "_name_ok", sep="")
  rownames(out) = row.names(x)
  out$scientific_name_ok = ROBITaxonomy::scientificname(obitaxdb, taxids)
  out$taxonomic_rank_ok = ROBITaxonomy::taxonomicrank(obitaxdb, taxids)
  return(out)
}

#' Add taxids to fasta
#' @param infile A fasta file with scientific names in headers.
#' \itemize{
#'     \item By default, the sequence identifier is used.
#'     \itemize{
#'         \item Underscore characters (_) are substituted by spaces.
#'         \item e.g. \code{>Leiopelma_hochstetteri length=658}}
#'     \item Alternatively, the k option specifies an attribute containing the taxon scientific name.
#'     \itemize{
#'         \item e.g.  \code{>seqX species=Leiopelma_hochstetteri;}
#'         \item Note that the \code{space} before the attribute and the \code{;} after is mandatory}}
#' @param taxo An obitools-formatted taxonomy database - a set of files (.adx,.ndx,.rdx,.tdx).
#' @return A new fasta with taxids in header.
#' @note Sequences for which taxids cannot be found are discarded
#' @examples
#' obiaddtaxids.Bas(infile = a, taxo="/media/sf_Documents/WORK/CIBIO/STATS_AND_CODE/TAXONOMIES/obitax_26-4-19",
#' out = "wTaxids.head.test.op2.fasta",k = "species")
#' @export
obiaddtaxids.Bas<-function(infile,taxo,k=NULL,out){
  cb <- function(line, proc) {cat(line, "\n")}
  if(is.null(k)){
    f<-process$new(command = "obiaddtaxids", args=c(infile,"-d",taxo),
                   echo_cmd = T,stdout=out)
    f$wait()
    d<-"SUCCESS!"
  }
  if(!is.null(k)){
    f<-process$new(command = "obiaddtaxids", args=c(infile,"-d",taxo, "-k", k),
                   echo_cmd = T,stdout=out)
    f$wait()
    d<-"SUCCESS!"
  }
}

#' Report taxonomic resolution attained for each taxon amplified in silico
#' @param ecopcroutput A data frame such as that produced by \code{ecoPCR.Bas}.
#' @param taxo An obitools-formatted taxonomy database. This can be an R object, generated by \code{\link[ROBITaxonomy]{read.taxonomy}}
#'     or \code{\link[bastools]{NCBI2obitaxonomy}}, or a set of files (.adx,.ndx,.rdx,.tdx).
#' @param taxLevel Taxonomic level to use to generate list.
#'     Accepted values are "superkingdom", "phylum", "class", "order", "family", "genus", "species" (default), "subspecies".
#' @return A data frame consisting of the species amplified, the taxid and the taxon resolution attained
#'     by the amplified fragment.
#' @examples
#' a2<-system.file("extdata", "45F-63R_4Mis_ALLVERTS_REFSEQ.ecopcroutput", package = "bastools")
#' b2<-substr(system.file("extdata", "obitax_26-4-19.ndx", package = "bastools"),1,str_locate(system.file("extdata", "obitax_26-4-19.ndx", package = "bastools"),"\\.")-1)
#' d2<-ecopcr.hit.table(a2,obitaxdb = b2,taxLevel = "class")
#' e2<-d2$ecopcroutput_uncleaned
#' test.res2<-ecopcroutput.res.Bas(e2,b2,taxLevel = "order")
#' plot.ecopcroutput.res.Bas(test.res2)
#' @export
#'
ecopcroutput.res.Bas<-function(ecopcroutput,obitaxdb,taxLevel="species"){
  if (class(obitaxdb)[1]!="obitools.taxonomy") obitaxdb2=ROBITaxonomy::read.taxonomy(obitaxdb)
  if (class(obitaxdb)[1]=="obitools.taxonomy") obitaxdb2=obitaxdb
  res=ROBIBarcodes::resolution(taxonomy = obitaxdb2,ecopcr = ecopcroutput)
  a<-as.data.frame(res)
  a$res<-gsub(x = a$res,pattern = "infraclass",replacement = "class")
  a$res<-gsub(x = a$res,pattern = "subfamily",replacement = "family")
  a$res<-gsub(x = a$res,pattern = "tribe",replacement = "family")
  a$res<-gsub(x = a$res,pattern = "suborder",replacement = "order")
  a$res<-gsub(x = a$res,pattern = "infraorder",replacement = "order")
  a$res<-gsub(x = a$res,pattern = "superfamily",replacement = "order")
  a$res<-gsub(x = a$res,pattern = "subclass",replacement = "class")
  a$res<-gsub(x = a$res,pattern = "subgenus",replacement = "genus")
  a$res<-gsub(x = a$res,pattern = "species subgroup",replacement = "species")
  a$res<-gsub(x = a$res,pattern = "cohort",replacement = "order")
  
  #K,P,C,O,F,G,S
  
  b<-cbind(ecopcroutput,a)
  
  taxongroup=paste0(taxLevel,"_name")
  e<-as.data.frame(table(b[,taxongroup],b$res))
  colnames(e)<-c("taxon","resolution","Freq")
  f<-e[e$Freq!=0,]
  return(f)
}

#' @export
obiaddtaxids.nodump.Bas<-function(infile,taxo,k=NULL,out){
  cb <- function(line, proc) {cat(line, "\n")}
  if(is.null(k)){
    f<-process$new(command = "obiaddtaxids", args=c(infile,"-d",taxo),
                   echo_cmd = T,stdout=out)
    f$wait()
    d<-"SUCCESS!"
  }
  if(!is.null(k)){
    f<-process$new(command = "obiaddtaxids", args=c(infile,"-d",taxo, "-k", k),
                   echo_cmd = T,stdout=out)
    f$wait()
    d<-"SUCCESS!"
  }
}

ecopcr2fasta<-function(mod_ecopcrout,out){
  message("Reminder: Converts all ecopcr results to fasta. This can probably be greatly reduced!")
  
  ecopcr<-data.table::fread(mod_ecopcrout,sep = "\t",data.table = F)
  
  ecopcr$seq.name<-paste0(ecopcr$AC, " taxid=", ecopcr$taxid,"; organism=", gsub(" ","_",ecopcr$species_name))
  ecopcr$seq.text<-ecopcr$sequence
  
  phylotools::dat2fasta(ecopcr[,c("seq.name","seq.text")],outfile = out)
}

add.res.bas2<-function(mod_ecopcrout_file,makeblastdb_exec="makeblastdb",ncbiTaxDir,blast_exec="blastn",obitaxdb, top=1){
  
  t1<-Sys.time()
  
  message("Adding resolution attained by fragment.
          Note: Runs a self-blast of inserts in ecopcroutput. From blast results;
          1) keeps only top hit (which will always be 100% identical) and anything within user-defined % (default=1) identity of top hit
          2) keeps only hits with >99% query cover
          3) bins remaining reads to lca & checks the rank of the lca
          
          i.e. using default top, 
          if a sequence has multiple hits that are all over 99% identity across 99% of the insert,
           and those hits are from different species within the same genus, the rank will be 
          assigned genus. If those hits are from the same species, the rank will be species.
          
          Note2: Perahps by running this multiple times with different 'tops'  we can plot resolution by % identity to tell
          us ... something.

          Note3: Need to keep this to uniques in future, veyr slow now (122k seqs in 30min)
          
          ")
  
  out="resolution_test_files"
  
  ecopcr2fasta(mod_ecopcrout_file,paste0(out,".fasta"))
  
  make.blastdb.bas(infasta = paste0(out,".fasta"),makeblastdb_exec = makeblastdb_exec,ncbiTaxDir = ncbiTaxDir,dbversion = 5)
  
  blast.min.bas(infastas = paste0(out,".fasta"),refdb = out,blast_exec = blast_exec)
  
  filter.blast(blastfile = paste0(out,".blast.txt"),ncbiTaxDir = ncbiTaxDir,out = paste0(out,".blast.filt.txt"),top=top,min_qcovs = 99)
  
  bin.blast.lite(filtered_blastfile = paste0(out,".blast.filt.txt"),ncbiTaxDir = ncbiTaxDir
                 ,obitaxdb = obitaxdb, out = paste0(out,".bins.txt"))           
  
  #add rank
  bins<-data.table::fread(paste0(out,".bins.txt"),data.table = F)
  bins$path<-paste(bins$K,bins$P,bins$C,bins$O,bins$F,bins$G,bins$S,sep = ";")
  temprank<-stringr::str_count(bins$path,";NA")
  temprank<-gsub(0,"species",temprank)
  temprank<-gsub(1,"genus",temprank)
  temprank<-gsub(2,"family",temprank)
  temprank<-gsub(3,"above_family",temprank)
  temprank<-gsub(4,"above_family",temprank)
  temprank<-gsub(5,"above_family",temprank)
  temprank<-gsub(6,"above_family",temprank)
  bins$resolution.bas<-temprank
  
  #merge with ecopcr output
  ecopcrout<-data.table::fread(mod_ecopcrout_file,data.table = F)  
  ecopcroutput<-merge(ecopcrout,bins[,c("qseqid","resolution.bas")],by.x = "AC",by.y = "qseqid")
  
  #write
  write.table(ecopcroutput,file = mod_ecopcrout_file,sep = "\t",append = F,quote = F,row.names = F)
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  
  message(paste0("Resolution added to ",mod_ecopcrout_file," in ",t3," mins."))
}



run.ecotaxspecificity<-function(infasta,ecopcrdb){
  test<-system2("ecotaxspecificity",c("-d",ecopcrdb, "-e", 1, infasta),wait = T,stdout = T)
}

google.make.MBC.sheet<-function(tokenDir,email,url,out,subsetlist){
  
  a<-getwd()
  setwd(tokenDir)
  googlesheets4::sheets_auth(email)
  setwd(a)
  
  master<-google.read.master.url(url)
  
  #change the common capitals issues
  colnames(master)<-gsub("Index_combination","index_combination",colnames(master))
  colnames(master)<-gsub("Min_length","min_length",colnames(master))
  colnames(master)<-gsub("Max_length","max_length",colnames(master))
  
  #remove duplicated columns
  if(length(grep("\\.\\.\\.",colnames(master),value = T))>0) {
    message("Removing duplicated columns")
    a<-grep("\\.\\.\\.",colnames(master),value = T)
    b<-do.call(rbind,stringr::str_split(a,"\\.\\.\\."))
    d<-b[duplicated(b[,1]),]
    e<-paste0(d[,1],"...",d[,2])
    master<-master[,!colnames(master) %in% e]
    colnames(master)<-gsub("\\.\\.\\..*","",colnames(master))
  }
  
  headers<-c("Sample_ID",	"Sample_Name",	"Sample_Plate",	"Sample_Well",	"I7_Index_ID",	"index",
             "I5_Index_ID",	"index2",	"Sample_Project",	"Description",	"index_combination",	"MID_F",	
             "MID_R",	"Primer_F",	"Primer_R",	"Filenames","min_length","max_length","ss_sample_id")
  
  for(i in 1:length(headers)){
    if(!headers[i] %in% colnames(master)) stop (c(headers[i]," missing: "))
  }
  
  message("Defaulting to putting all index_combinations to 'used'")
  master$index_combination<-"used"
  
  #subset sheet
  for(i in 1:length(subsetlist)){
    if(!names(subsetlist)[i] %in% colnames(master)) stop (c(names(subsetlist)[i]," missing: "))
  }
  
  ms_ss<-subset_mastersheet(master,subsetlist)
  
  #check subset
  message("Check the numbers below make sense")
  print(master_xtabs(ms_ss,columns=names(subsetlist)))
  
  #make ss_sample_id Sample_Name
  ms_ss$Sample_Name<-ms_ss$ss_sample_id
  
  #use only the required headers (exclude ss_sample_id)
  ms_ss<-ms_ss[,headers[-length(headers)]]
  
  #order correctly
  ms_ss<-ms_ss[c(headers[-length(headers)])]
  
  ms_ss$extra_information<-""
  
  #write file
  write.table(ms_ss,out,sep = "\t",quote = F,row.names = F)
  message(paste0("file saved as ",paste0(getwd(),"/",out)))
  
  return(ms_ss)
  
}


add.stats.ecopcroutput<-function(ecopcroutput,ncbiTaxDir,Ta=NULL,add.3pmm=F,Pf,Pr){
  
  #change taxid to taxids, for next function
  colnames(ecopcroutput)<-gsub("taxid","taxids",colnames(ecopcroutput))
  
  #add full amplicon, including primers
  message("Adding full amplicon to fullseq, including primers")
  ecopcroutput$fullseq<-paste0(ecopcroutput$forward_match,ecopcroutput$sequence,ecopcroutput$reverse_match)
  
  #add taxonomy using taxonkit
  message("Adding full taxonomy")
  ecopcroutput<-add.lineage.df(df = ecopcroutput,ncbiTaxDir = ncbiTaxDir,as.taxids = F)
  
  #add 3 prime mismtaches to ecopcroutput 
  #adds 4 columns: 3prime mismatches for half the primer, 3prime mismatches for last 6bp only,
  #for both forward and reverse primers
  if(add.3pmm) ecopcroutput<-add.3pmms(ecopcroutput,Pf,Pr) 
  
  #add tm
  #adds 6 columns: tm for full primer, 3prime half, 3primer6; for both primers. Uses basic equation
  ecopcroutput<-add.tm.ecopcroutput(ecopcroutput)
  
  #add diff tm: the difference in Tm between forward and reverse primers
  message("Calculating the difference in Tm between forward and reverse primer binding sites")
  ecopcroutput$diff_tm<-ecopcroutput$fTms-ecopcroutput$rTms
  
  #calculate gc% for primer binding sites (forward and reverse)
  ecopcroutput<-add.gc.ecopcroutput(ecopcroutput)
  
  #gc clamp present?###########
  message("Add function for presence of GC clamp?")
  # gc.clamp<-function(primerVec){
  #   gc<-vector()
  #   for (i in 1:length(primerVec)){
  #     gc[i]<-(stringr::str_count(primerVec[i],"G")+stringr::str_count(primerVec[i],"C"))/
  #       (nchar(primerVec[i]))*100
  #   }
  #   return(gc)
  # }
  
  #add diff ta tm
  if(!is.null(Ta)){
    message("Calculating difference between primer binding site and annealing temperature used")
    ecopcroutput$tm.ta_fw<-ecopcroutput$fTms-Ta
    ecopcroutput$tm.ta_rv<-ecopcroutput$rTms-Ta  
  }
  
  #add Tm of last 6 bp of fw primer divided by overall Tm (%)
  message("Calculating Tm of last 6 bp of fw primer divided by overall Tm (%)")
  ecopcroutput$tm_fw_3p6_perc<-ecopcroutput$fTms3prime6/ecopcroutput$fTms*100
  ecopcroutput$tm_rv_3p6_perc<-ecopcroutput$rTms3prime6/ecopcroutput$rTms*100
  
  #remove extraneous columns and sort
  ecopcroutput<- subset(ecopcroutput,select = -c(rank,species,species_name,genus,genus_name,family,family_name,
                                                  superkingdom, superkingdom_name,forward_tm,reverse_tm))
  
  ecopcroutput<-ecopcroutput %>% select(AC,taxids,K,P,C,O,F,G,S,amplicon_length,everything())
  ecopcroutput<-ecopcroutput %>% select(-sequence,sequence)
}

remove.google.dups<-function(master_sheet){
  if(length(grep("\\.\\.\\.",colnames(master_sheet),value = T))>0) {
    message("Removing duplicated columns")
    a<-grep("\\.\\.\\.",colnames(master_sheet),value = T)
    b<-do.call(rbind,stringr::str_split(a,"\\.\\.\\."))
    d<-b[duplicated(b[,1]),]
    e<-paste0(d[,1],"...",d[,2])
    master_sheet<-master_sheet[,!colnames(master_sheet) %in% e]
    colnames(master_sheet)<-gsub("\\.\\.\\..*","",colnames(master_sheet))
  }
  master_sheet
}

#changing sample_type too to match functions
fix.google.colnames<-function(master_sheet){
  colnames(master_sheet)<-gsub("Sample_Type","sample_type",colnames(master_sheet))
  master_sheet
}
