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

  