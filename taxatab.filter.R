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
  taxatab3<-cbind(taxon=taxatab2$taxon,taxatab2[,-1][,colSums(taxatab2[,-1])!=0])
}

rm.0readOTUSam<-function(taxatab){
  taxatab2<-taxatab[rowSums(taxatab[,-1])!=0,]
  taxatab3<-cbind(OTU=taxatab2$OTU,taxatab2[,-1][,colSums(taxatab2[,-1])!=0])
}

filter.dxns<-function(taxatab,filter_dxn=50){
  taxatab2<-taxatab[,-1]
  taxatab2[taxatab2<filter_dxn] <- 0
  taxatab2<-cbind(taxon=taxatab$taxon,taxatab2)  
  taxatab2<-rm.0readtaxSam(taxatab2)
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

negs.stats<-function(taxatab,ms_ss,real,ex_hominidae=T){
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
      
      if(length(colnames(taxatab.negs.list[grep(unique(negs$sample_type)[i],names(taxatab.negs.list))][[1]]))>0){
      print(taxatab.negs.list[grep(unique(negs$sample_type)[i],names(taxatab.negs.list))][[1]])}
  }

  return(taxatab.negs.list)
  } else message("No negatives found with reads")
  
}




taxatab.filter<-function(taxatab,taxa_pc=0.15,sample_pc=0.15,absolute_in_rep=0,absolute_detection=100,
                       sumreps=T,pcr_neg_pattern="NC",extraction_neg_pattern="EC",real_sample_pattern="F"){
  message("Assumes taxatab is a matrix with samples (in DMP format) as colnames,
          taxa as rownames and reads as values. taxatab.filter2 expects pcr_negs and extraction_negs
          as character vectors, while taxatab.filter expects a single pattern to distinguish these.
          
          THIS IS ABSOLETE")
  
  preprocess.samples<-do.counts(taxatab[,grep(real_sample_pattern,colnames(taxatab))],"preprocess-real")
  preprocess.pcr.negs<-do.counts(taxatab[,grep(pcr_neg_pattern,colnames(taxatab))],"preprocess-pcr.negs")
  preprocess.ext.negs<-do.counts(taxatab[,grep(extraction_neg_pattern,colnames(taxatab))],"preprocess-ext.negs")
  
  #remove detections where less than x% of total reads for taxon
  message("Applying taxa_pc filter")
  taxatab<-taxatab[!rowSums(taxatab)==0,]
  taxatab.PCS<-sweep(taxatab, MARGIN = 1, STATS = rowSums(taxatab), FUN = "/")*100
  taxatab.PCS[taxatab.PCS<taxa_pc]<-0 #change this value
  taxatab[taxatab.PCS==0]<-0
  
  taxa.pc.samples<-do.counts(taxatab[,grep(real_sample_pattern,colnames(taxatab))],"taxa.pc-real")
  taxa.pc.pcr.negs<-do.counts(taxatab[,grep(pcr_neg_pattern,colnames(taxatab))],"taxa.pc-pcr.negs")
  taxa.pc.ext.negs<-do.counts(taxatab[,grep(extraction_neg_pattern,colnames(taxatab))],"taxa.pc-ext.negs")
  
  #remove detections where less than x% of total reads for sample
  message("Applying sample_pc filter")
  taxatab<-taxatab[,!colSums(taxatab)==0]
  taxatab.PCS<-sweep(taxatab, MARGIN = 2, STATS = colSums(taxatab), FUN = "/")*100
  taxatab.PCS[taxatab.PCS<sample_pc]<-0  
  taxatab[taxatab.PCS==0]<-0

  sample.pc.samples<-do.counts(taxatab[,grep(real_sample_pattern,colnames(taxatab))],"sample.pc-real")
  sample.pc.pcr.negs<-do.counts(taxatab[,grep(pcr_neg_pattern,colnames(taxatab))],"sample.pc-pcr.negs")
  sample.pc.ext.negs<-do.counts(taxatab[,grep(extraction_neg_pattern,colnames(taxatab))],"sample.pc-ext.negs")
  
  #Replicate filter
  message("Applying replicate filter")
  taxatab<-rm.detections(df = taxatab,absolute_in_rep,sumreps = sumreps)
  
  rep.filter.samples<-do.counts(taxatab[,grep(real_sample_pattern,colnames(taxatab))],"rep.filter-real")
  rep.filter.pcr.negs<-do.counts(taxatab[,grep(pcr_neg_pattern,colnames(taxatab))],"rep.filter-pcr.negs")
  rep.filter.ext.negs<-do.counts(taxatab[,grep(extraction_neg_pattern,colnames(taxatab))],"rep.filter-ext.negs")
  
  #Absolute filter
  message("Applying absolute filter")
  taxatab[taxatab<absolute_detection]<-0   
  
  absolute.filter.samples<-do.counts(taxatab[,grep(real_sample_pattern,colnames(taxatab))],"absolute.filter-real")
  absolute.filter.pcr.negs<-do.counts(taxatab[,grep(pcr_neg_pattern,colnames(taxatab))],"absolute.filter-pcr.negs")
  absolute.filter.ext.negs<-do.counts(taxatab[,grep(extraction_neg_pattern,colnames(taxatab))],"absolute.filter-ext.negs")
  
  ###Taxa lists
  taxatab<-taxatab[,colSums(taxatab)!=0]
  taxatab<-taxatab[rowSums(taxatab)!=0,]
  #Taxa in real samples
  if(length(colnames((taxatab[,grep(real_sample_pattern,colnames(taxatab))])))>1){
    taxa.list.samples<-as.data.frame(rowSums(taxatab[,grep(real_sample_pattern,colnames(taxatab))]))}
  if(length(colnames((taxatab[,grep(real_sample_pattern,colnames(taxatab))])))==1){
    taxa.list.samples<-as.data.frame(taxatab[,grep(real_sample_pattern,colnames(taxatab))])}
  if(length(colnames((taxatab[,grep(real_sample_pattern,colnames(taxatab))])))==0){
    taxa.list.samples<-as.data.frame("NA")
    rownames(taxa.list.samples)<-"NA"}
  colnames(taxa.list.samples)<-c("count")
  
  #Taxa in pcr negs
  if(length(colnames((taxatab[,grep(pcr_neg_pattern,colnames(taxatab))])))>1){
    taxa.list.pcr.negs<-as.data.frame(rowSums(taxatab[,grep(pcr_neg_pattern,colnames(taxatab))]))}
  if(length(colnames((taxatab[,grep(pcr_neg_pattern,colnames(taxatab))])))==1){
    taxa.list.pcr.negs<-as.data.frame(taxatab[,grep(pcr_neg_pattern,colnames(taxatab))])}
  if(length(colnames((taxatab[,grep(pcr_neg_pattern,colnames(taxatab))])))==0){
    taxa.list.pcr.negs<-as.data.frame("NA")
    rownames(taxa.list.pcr.negs)<-"NA"}
  colnames(taxa.list.pcr.negs)<-c("count")
  
  #Taxa in ext negs
  if(length(colnames((taxatab[,grep(extraction_neg_pattern,colnames(taxatab))])))>1){
    taxa.list.ext.negs<-as.data.frame(rowSums(taxatab[,grep(extraction_neg_pattern,colnames(taxatab))]))}
  if(length(colnames((taxatab[,grep(extraction_neg_pattern,colnames(taxatab))])))==1){
    taxa.list.ext.negs<-as.data.frame(taxatab[,grep(extraction_neg_pattern,colnames(taxatab))])}
  if(length(colnames((taxatab[,grep(extraction_neg_pattern,colnames(taxatab))])))==0){
    taxa.list.ext.negs<-as.data.frame("NA")
    rownames(taxa.list.ext.negs)<-"NA"}
  colnames(taxa.list.ext.negs)<-c("count")
  
  #Combine taxa lists
  taxa.list.final<-as.data.frame(unique(c(rownames(taxa.list.samples),
                                          rownames(taxa.list.pcr.negs),rownames(taxa.list.ext.negs))))
  colnames(taxa.list.final)<-"taxa"
  taxa.list.samples$taxa<-rownames(taxa.list.samples)
  taxa.list.pcr.negs$taxa<-rownames(taxa.list.pcr.negs)
  taxa.list.ext.negs$taxa<-rownames(taxa.list.ext.negs)
  taxa.list.final<-merge(x = taxa.list.final,y = taxa.list.samples,all = T,by = "taxa")
  colnames(taxa.list.final)<-gsub("count","real",colnames(taxa.list.final))
  taxa.list.final<-merge(x = taxa.list.final,y = taxa.list.pcr.negs,all = T,by = "taxa")
  colnames(taxa.list.final)<-gsub("count","pcr.negs",colnames(taxa.list.final))
  taxa.list.final<-merge(x = taxa.list.final,y = taxa.list.ext.negs,all = T,by = "taxa")
  colnames(taxa.list.final)<-gsub("count","ext.negs",colnames(taxa.list.final))
  taxa.list.final<-taxa.list.final[!taxa.list.final$taxa=="NA",]
  taxa.list.final$taxa<-as.character(taxa.list.final$taxa)
  taxa.list.final$real<-as.numeric(taxa.list.final$real)
  taxa.list.final$pcr.negs<-as.numeric(taxa.list.final$pcr.negs)
  taxa.list.final$ext.negs<-as.numeric(taxa.list.final$ext.negs)
  taxa.list.final[is.na(taxa.list.final)]<-0
  taxa.list.final<-taxa.list.final[rowSums(taxa.list.final[,2:4])!=0,]
  
  
  #combine counts
  Totals<-rbind(preprocess.samples, taxa.pc.samples,sample.pc.samples,rep.filter.samples, 
                absolute.filter.samples,
                preprocess.pcr.negs,taxa.pc.pcr.negs,sample.pc.pcr.negs,rep.filter.pcr.negs,
                absolute.filter.pcr.negs,
                preprocess.ext.negs,
                taxa.pc.ext.negs,sample.pc.ext.negs,rep.filter.ext.negs,absolute.filter.ext.negs)
               
  #Output
  output<-list()
  output[[1]]<-Totals
  output[[2]]<-taxa.list.final
  output[[3]]<-paste0("taxa_pc",taxa_pc,";sample_pc",sample_pc,";absolute_in_rep",absolute_in_rep,
                             ";absolute_detection",absolute_detection,";sumreps",sumreps)
  output[[4]]<-taxatab
  return(output)
}



taxatab.filter2<-function(taxatab,taxa_pc=0.15,sample_pc=0.15,absolute_in_rep=6000,absolute_detection=100,
                          sumreps=T,pcr_negs=NULL,extraction_negs=NULL,real_samples){
  message("Assumes taxatab is a dataframe with samples (in DMP format) as colnames,
          taxa as rownames and reads as values. taxatabfilter2 expects pcr_negs and extraction_negs
          as character vectors, while taxatabfilter expects a single pattern to distinguish these.")
  
  preprocess.samples<-do.counts(taxatab[,real_samples],"preprocess-real")
  preprocess.pcr.negs<-do.counts(taxatab[,pcr_negs],"preprocess-pcr.negs")
  preprocess.ext.negs<-do.counts(taxatab[,extraction_negs],"preprocess-ext.negs")
  
  #remove detections where less than x% of total reads for taxon
  message("Applying taxa_pc filter")
  taxatab<-taxatab[!rowSums(taxatab)==0,]
  taxatab.PCS<-sweep(taxatab, MARGIN = 1, STATS = rowSums(taxatab), FUN = "/")*100
  taxatab.PCS[taxatab.PCS<taxa_pc]<-0 
  taxatab[taxatab.PCS==0]<-0
  
  taxa.pc.samples<-do.counts(taxatab[,real_samples],"taxa.pc-real")
  taxa.pc.pcr.negs<-do.counts(taxatab[,pcr_negs],"taxa.pc-pcr.negs")
  taxa.pc.ext.negs<-do.counts(taxatab[,extraction_negs],"taxa.pc-ext.negs")
  
  #remove detections where less than x% of total reads for sample
  message("Applying sample_pc filter")
  taxatab<-taxatab[,!colSums(taxatab)==0]
  taxatab.PCS<-sweep(taxatab, MARGIN = 2, STATS = colSums(taxatab), FUN = "/")*100
  taxatab.PCS[taxatab.PCS<sample_pc]<-0  
  taxatab[taxatab.PCS==0]<-0
  
  real_samples<-real_samples[real_samples %in% colnames(taxatab)]
  pcr_negs<-pcr_negs[pcr_negs %in% colnames(taxatab)]
  extraction_negs<-extraction_negs[extraction_negs %in% colnames(taxatab)]
  
  sample.pc.samples<-do.counts(taxatab[,real_samples],"sample.pc-real")
  sample.pc.pcr.negs<-do.counts(taxatab[,pcr_negs],"sample.pc-pcr.negs")
  sample.pc.ext.negs<-do.counts(taxatab[,extraction_negs],"sample.pc-ext.negs")
  
  #Replicate filter
  message("Applying replicate filter")
  taxatab<-rm.detections(df = taxatab,absolute_in_rep,sumreps = sumreps)
  
  if(sumreps==T){
    real_samples<-unique(gsub("(.*)-(.*)-(.*)-(.*)-(.*)$", "\\1-\\2-\\3", real_samples))
    real_samples<-real_samples[real_samples %in% colnames(taxatab)]
    pcr_negs<-unique(gsub("(.*)-(.*)-(.*)-(.*)-(.*)$", "\\1-\\2-\\3", pcr_negs))
    pcr_negs<-pcr_negs[pcr_negs %in% colnames(taxatab)]
    extraction_negs<-unique(gsub("(.*)-(.*)-(.*)-(.*)-(.*)$", "\\1-\\2-\\3", extraction_negs))
    extraction_negs<-extraction_negs[extraction_negs %in% colnames(taxatab)]
  } 
  
  if(sumreps==F){
    real_samples<-real_samples[real_samples %in% colnames(taxatab)]
    pcr_negs<-pcr_negs[pcr_negs %in% colnames(taxatab)]
    extraction_negs<-extraction_negs[extraction_negs %in% colnames(taxatab)]
  } 
  
  rep.filter.samples<-do.counts(taxatab[,real_samples],"rep.filter-real")
  rep.filter.pcr.negs<-do.counts(taxatab[,pcr_negs],"rep.filter-pcr.negs")
  rep.filter.ext.negs<-do.counts(taxatab[,extraction_negs],"rep.filter-ext.negs")
  
  #Absolute filter
  message("Applying absolute filter")
  taxatab[taxatab<absolute_detection]<-0   
  
  absolute.filter.samples<-do.counts(taxatab[,real_samples],"absolute.filter-real")
  absolute.filter.pcr.negs<-do.counts(taxatab[,pcr_negs],"absolute.filter-pcr.negs")
  absolute.filter.ext.negs<-do.counts(taxatab[,extraction_negs],"absolute.filter-ext.negs")
  
  ###Taxa lists
  taxatab<-taxatab[,colSums(taxatab)!=0]
  taxatab<-taxatab[rowSums(taxatab)!=0,]
  
  real_samples<-real_samples[real_samples %in% colnames(taxatab)]
  pcr_negs<-pcr_negs[pcr_negs %in% colnames(taxatab)]
  extraction_negs<-extraction_negs[extraction_negs %in% colnames(taxatab)]
  
  #Taxa in real samples
  if(length(real_samples)>1){
    taxa.list.samples<-as.data.frame(rowSums(taxatab[,real_samples]))}
  if(length(real_samples)==1){
    taxa.list.samples<-as.data.frame(taxatab[,real_samples])}
  if(length(real_samples)==0){
    taxa.list.samples<-as.data.frame("NA")
    rownames(taxa.list.pcr.negs)<-("NA")
    }
  colnames(taxa.list.samples)<-c("count")
  
  #Taxa in pcr negs
  if(length(pcr_negs)>1){
    taxa.list.pcr.negs<-as.data.frame(rowSums(taxatab[,pcr_negs]))}
  if(length(pcr_negs)==1){
    taxa.list.pcr.negs<-as.data.frame(taxatab[,pcr_negs])
    rownames(taxa.list.pcr.negs)<-rownames(taxatab)}
  if(length(pcr_negs)==0){
    taxa.list.pcr.negs<-as.data.frame("NA")
    rownames(taxa.list.pcr.negs)<-("NA")
  }
  colnames(taxa.list.pcr.negs)<-c("count")
  
  #Taxa in ext negs
  if(length(extraction_negs)>1){
    taxa.list.ext.negs<-as.data.frame(rowSums(taxatab[,extraction_negs]))}
  if(length(extraction_negs)==1){
    taxa.list.ext.negs<-as.data.frame(taxatab[,extraction_negs])
    rownames(taxa.list.ext.negs)<-rownames(taxatab)}
  if(length(extraction_negs)==0){
    taxa.list.ext.negs<-as.data.frame("NA")
    rownames(taxa.list.pcr.negs)<-("NA")
  }
  colnames(taxa.list.ext.negs)<-c("count")
  
  #Combine taxa lists
  taxa.list.final<-as.data.frame(unique(c(rownames(taxa.list.samples),
                                          rownames(taxa.list.pcr.negs),rownames(taxa.list.ext.negs))))
  colnames(taxa.list.final)<-"taxa"
  taxa.list.samples$taxa<-rownames(taxa.list.samples)
  taxa.list.pcr.negs$taxa<-rownames(taxa.list.pcr.negs)
  taxa.list.ext.negs$taxa<-rownames(taxa.list.ext.negs)
  taxa.list.final<-merge(x = taxa.list.final,y = taxa.list.samples,all = T,by = "taxa")
  colnames(taxa.list.final)<-gsub("count","real",colnames(taxa.list.final))
  taxa.list.final<-merge(x = taxa.list.final,y = taxa.list.pcr.negs,all = T,by = "taxa")
  colnames(taxa.list.final)<-gsub("count","pcr.negs",colnames(taxa.list.final))
  taxa.list.final<-merge(x = taxa.list.final,y = taxa.list.ext.negs,all = T,by = "taxa")
  colnames(taxa.list.final)<-gsub("count","ext.negs",colnames(taxa.list.final))
  taxa.list.final<-taxa.list.final[!taxa.list.final$taxa=="NA",]
  taxa.list.final$taxa<-as.character(taxa.list.final$taxa)
  taxa.list.final$real<-as.numeric(taxa.list.final$real)
  taxa.list.final$pcr.negs<-as.numeric(taxa.list.final$pcr.negs)
  taxa.list.final$ext.negs<-as.numeric(taxa.list.final$ext.negs)
  taxa.list.final[is.na(taxa.list.final)]<-0
  taxa.list.final<-taxa.list.final[rowSums(taxa.list.final[,2:4])!=0,]
  
  
  #combine counts
  Totals<-rbind(preprocess.samples, taxa.pc.samples,sample.pc.samples,rep.filter.samples, 
                absolute.filter.samples,
                preprocess.pcr.negs,taxa.pc.pcr.negs,sample.pc.pcr.negs,rep.filter.pcr.negs,
                absolute.filter.pcr.negs,
                preprocess.ext.negs,
                taxa.pc.ext.negs,sample.pc.ext.negs,rep.filter.ext.negs,absolute.filter.ext.negs)
  
  #Output
  output<-list()
  output[[1]]<-Totals
  output[[2]]<-taxa.list.final
  output[[3]]<-paste0("taxa_pc",taxa_pc,";sample_pc",sample_pc,";absolute_in_rep",absolute_in_rep,
                      ";absolute_detection",absolute_detection,";sumreps",sumreps)
  output[[4]]<-taxatab
  return(output)
}

#remove contaminant taxa from samples, based on occurrences in negatives
remove.contaminant.taxa<-function(master_sheet,taxatab,negatives,group.code){
  
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
        contaminations<-cbind(taxon=contaminations$taxon,contaminations[,-1,drop=F][,colSums(contaminations[,-1,drop=F])>0,drop=F])
        taxatab[taxatab$taxon %in% negtaxa,colnames(taxatab) %in% group.samples]<-0
        
        sumafter<-sum(taxatab[,-1])
        
        message(paste(sumb4-sumafter,"reads removed from",sum(contaminations[,-1,drop=F]>0),"detection(s) in",
                      length(colnames(contaminations[,-1,drop=F])),"sample(s), of which",sumnegtaxa,"were in the negative. See table below for details:"))
        print(contaminations)
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

sumreps<-function(taxatab,master_sheet){
  #this gsub will not always work. e.g where sample names were chosen in google sheet, not exactly the same as sample
  #DNA_samples<-unique(gsub("(.*)-(.*)-(.*)-(.*)$", "\\1", colnames(taxatab[,-1])))
  
  DNA_samples<-master_sheet[match(colnames(taxatab[,-1]),table = master_sheet[,"ss_sample_id"]),"Sample_Name"]
  
  taxatab2<-taxatab
  colnames(taxatab2)<-c("taxon",DNA_samples)
  
  summed<-list()
  for(i in 1:length(unique(DNA_samples))){
    taxatab.temp<-as.data.frame(taxatab2[,grep(unique(DNA_samples)[i],colnames(taxatab2)),drop=F])
    if(length(colnames(taxatab.temp))>1) taxatab.temp$sum<-rowSums(taxatab.temp) else taxatab.temp$sum<-taxatab.temp[,1]
    
    taxatab.temp2<-taxatab.temp[,"sum",drop=F]
    
    colnames(taxatab.temp2)<-gsub("sum",unique(DNA_samples)[i],colnames(taxatab.temp2))
    
    summed[[i]]<-taxatab.temp2
  }
    
  taxatab.out<-cbind(taxon=taxatab$taxon,do.call(cbind, summed))
}


adonis.bas<-function(taxatab,master_sheet,factor1,samLevel="ss_sample_id",stratum=NULL){
  
  taxatab2<-binarise.taxatab(taxatab)
  distance_matrix<-taxatab2bray(taxatab2)
  
  if(samLevel=="ss_sample_id") master_sheet2<-master_sheet[master_sheet$ss_sample_id  %in% colnames(taxatab),]
  
  if(samLevel=="Sample_Name") {
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
