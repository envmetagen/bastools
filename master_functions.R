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
    write.table(taxatables[[i]],gsub(".txt",".tf.txt",files[i]),row.names = F,quote = F,sep = "\t")
  }
}

taxon.filter.solo.df<-function(taxatab,taxonpc=0.1){
  message("Applying taxon_pc filter.")
  taxatab2<-taxatab[,-1]
  a<-sum(taxatab2)
  dxns1<-sum(taxatab2>0)
  taxatab.PCS<-sweep(taxatab2, MARGIN = 1, STATS = rowSums(taxatab2), FUN = "/")*100
  taxatab.PCS[taxatab.PCS<taxonpc]<-0 
  taxatab2[taxatab.PCS==0]<-0
  b<-sum(taxatab2)
  dxns2<-sum(taxatab2>0)
  taxatab3<-cbind(taxon=taxatab$taxon,taxatab2)
  message(paste("Using filter of", taxonpc, "%. reads removed:",a-b,"from",a,"; detections removed:",dxns1-dxns2,"from",dxns1))
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
  dxns1<-sum(taxatab>0)
  taxatab[taxatab.PCS==0]<-0
  b<-sum(taxatab)
  dxns2<-sum(taxatab>0)
  message(paste("Using filter of", samplepc, "%. reads removed:",a-b,"from",a,"; detections removed:",dxns1-dxns2,"from",dxns1))
  
  taxatab.out<-cbind(taxon,taxatab)
  taxatab.out<-rm.0readtaxSam(taxatab.out)
}

rm.0readtaxSam<-function(taxatab){
  taxatab2<-taxatab[rowSums(taxatab[,-1,drop=F])!=0,]
  taxatab3<-cbind(taxon=taxatab2[,1],taxatab2[,-1,drop=F][,colSums(taxatab2[,-1,drop=F])!=0,drop=F])
}

rm.0readOTUSam<-function(taxatab){
  taxatab2<-taxatab[rowSums(taxatab[,-1])!=0,]
  taxatab3<-cbind(OTU=taxatab2$OTU,taxatab2[,-1][,colSums(taxatab2[,-1])!=0])
}

bas.merge.taxatabs<-function(taxatabs){
  
  if(TRUE %in% duplicated(taxatabs)) stop("taxatabs provided are not unique")
  
  require(tidyverse)
  
  
  taxatabs.list<-list()
  counts<-data.frame(file=taxatabs,reads=0,taxa=0,samples=0)
  
  for(i in 1:length(taxatabs)){
    taxatabs.list[[i]]<-data.table::fread(taxatabs[i],sep = "\t",data.table = F)
    colnames(taxatabs.list[[i]])[1]<-"taxon" #this is for handling OTUtabs
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


filter.dxns<-function(taxatab,filter_dxn=50,rm.empty.taxsam=T){
  taxatab2<-taxatab[,-1]
  reads1<-sum(taxatab2)
  dxns1<-sum(taxatab2>0)
  taxatab2[taxatab2<filter_dxn] <- 0
  reads2<-sum(taxatab2)
  dxns2<-sum(taxatab2>0)
  taxatab2<-cbind(taxon=taxatab$taxon,taxatab2)  
  if(rm.empty.taxsam) taxatab2<-rm.0readtaxSam(taxatab2)
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
  rangetab<-suppressWarnings(apply(taxatab2,1,function(x) range(x,na.rm = T)))
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
        print(taxatab.negs.list[grep(unique(negs$sample_type)[i],names(taxatab.negs.list))][[1]],row.names = F,right = F)
        }
      }
    }
    
    return(taxatab.negs.list)
  } else message("No negatives found with reads")
  
}

#remove contaminant taxa from samples, based on occurrences in negatives
remove.contaminant.taxa<-function(master_sheet,taxatab,negatives,group.codes,printcontaminations=T,remove.entire.dataset=T,rm.only.less.than=T){
  
  if(rm.only.less.than==T)  message("Only removing detections which had a read count less than the negative")
  
    taxatabX<-taxatab
  
    totalstart<-sum(taxatabX[,-1,drop=F])
    
    negtaxaList<-list() #for use if removing neg taxa from entire dataset
    
    for(j in 1:length(negatives)){
      
      negdf<-negatives[[j]]
      
      negative_type<-names(negatives)[j]
    
      if(length(negdf)!=0){
    
        for(i in 2:length(colnames(negdf))){
          neg<-colnames(negdf)[i]
          negtaxa<-negdf[,c("taxon",neg)]
          sumnegtaxa<-sum(negtaxa[,-1,drop=F])
          negtaxa<-as.character(negtaxa[rowSums(negtaxa[,-1,drop=F])>0,"taxon"])
          
          negtaxaList[[i]]<-negtaxa
          
          #get group of a sample
          group.id<-master_sheet[grep(neg,master_sheet[,"ss_sample_id"]),group.codes[j]]
          #get other samples in group
          group.samples<-master_sheet[master_sheet[,group.codes[j]]==group.id,"ss_sample_id"]
          group.samples<-group.samples[!is.na(group.samples)]
          
          message(paste("Based on",negative_type,neg, ", Removing detections of"))
          print(negtaxa)
          message("if it occurred in any samples belonging to ",group.codes[j],": ",group.id)
          
          sumb4<-sum(taxatabX[,-1,drop=F])
          dxnsbefore<-sum(taxatabX[,-1,drop=F]>0)
          
          #put taxon counts to 0
          contaminations<-cbind(taxon=taxatabX$taxon[taxatabX$taxon %in% negtaxa],
                                              taxatabX[taxatabX$taxon %in% negtaxa,colnames(taxatabX) %in% group.samples,drop=F])
          contaminations<-cbind(taxon=contaminations[,1],contaminations[,-1,drop=F][,colSums(contaminations[,-1,drop=F])>0,drop=F])
          
            if(rm.only.less.than==T){
        
              negX<-taxatabX[taxatabX$taxon %in% negtaxa,c("taxon",neg)][,-1]
        
                for(h in 1:length(negtaxa)){
                  taxatab2<-taxatabX[taxatabX$taxon==negtaxa[h],colnames(taxatabX) %in% group.samples,drop=F]
                  taxatab2[taxatab2<(negX[h]+1)]<-0
                  taxatab3<-cbind(taxon=taxatabX[taxatabX$taxon==negtaxa[h],]$taxon,taxatab2)
                  taxatabX[taxatabX$taxon==negtaxa[h],colnames(taxatabX) %in% colnames(taxatab3)]<-taxatab3[1,] 
                }
            } else { taxatabX[taxatabX$taxon %in% negtaxa,colnames(taxatabX) %in% group.samples]<-0 }
        
         sumafter<-sum(taxatabX[,-1])
         dxnsafter<-sum(taxatabX[,-1,drop=F]>0)
        
         message(paste(sumb4-sumafter,"reads removed from",dxnsbefore-dxnsafter,"detection(s) in",
                      length(colnames(contaminations[,-1,drop=F])),"sample(s), of which",sumnegtaxa,"were in the negative. See table below for details:"))
        
         if(printcontaminations==T) {
          contaminations2<-as.data.frame(t(contaminations))
          colnames(contaminations2)<-contaminations$taxon
          contaminations2<-contaminations2[-1,]
          print(contaminations2,right = F,row.names=T,quote = F)
         }else(message("Not printing contamination table"))
      }
    }
    }
    
  taxatabX<-rm.0readtaxSam(taxatabX)
  totalend<-sum(taxatabX[,-1,drop=F])
  message(paste("***********A total of", totalstart-totalend, "reads removed"))
  
  if(remove.entire.dataset==T) {
    
    if(rm.only.less.than==T){
      
      message("****Furthermore, removing the following contaminant taxa from entire dataset if fewer reads than in negatives")
      all.neg.taxa<-do.call(c,negtaxaList)
      print(all.neg.taxa)
      
      for(h in 1:length(negtaxa)){
        taxatab2<-taxatabX[taxatabX$taxon==negtaxa[h],-1,drop=F]
        if(nrow(taxatab2)>0){
          taxatab2[taxatab2<(negX[h]+1)]<-0
          taxatab3<-cbind(taxon=taxatabX[taxatabX$taxon==negtaxa[h],]$taxon,taxatab2)
          taxatabX[taxatabX$taxon==negtaxa[h],colnames(taxatabX) %in% colnames(taxatab3)]<-taxatab3[1,] 
        }
      }
      
      taxatabX<-rm.0readtaxSam(taxatabX)
      message(paste("***********A further ", totalend-sum(taxatabX[,-1,drop=F]), "reads removed")) 
      
    } else {
      message("****Furthermore, removing the following contaminant taxa from entire dataset")
      all.neg.taxa<-do.call(c,negtaxaList)
      print(all.neg.taxa)
      taxatabX<-taxatabX[!taxatabX$taxon %in% all.neg.taxa,]
      taxatabX<-rm.0readtaxSam(taxatab)
      message(paste("***********A further ", totalend-sum(taxatabX[,-1,drop=F]), "reads removed")) 
    }
  }
  
  return(taxatabX)
}

#keep only xLevel assignments
keep.below.xLevel.assigns<-function(taxatab,xLevel="species",rm.trailing.NA=F){
  if(!xLevel %in% c("species","genus","family","order")) stop("Only allowable at genus, family, order or species level")
  message("Reminder: this changes 'unknown' and 'collapsed' to 'NA'")
  
  taxatab$taxon<-gsub("unknown","NA",taxatab$taxon)
  taxatab$taxon<-gsub("collapsed","NA",taxatab$taxon)
  
  message("Removing NAs and no_hits")
  taxatab<-taxatab[taxatab$taxon!="NA;NA;NA;NA;NA;NA;NA",]
  taxatab<-taxatab[taxatab$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]

  if(length(grep(";NA;NA;NA;NA;NA;NA$",taxatab$taxon))>0) taxatab<-taxatab[-grep(";NA;NA;NA;NA;NA;NA$",taxatab$taxon),]
  if(length(grep(";NA;NA;NA;NA;NA$",taxatab$taxon))>0) taxatab<-taxatab[-grep(";NA;NA;NA;NA;NA$",taxatab$taxon),]
  if(length(grep(";NA;NA;NA;NA$",taxatab$taxon))>0) taxatab<-taxatab[-grep(";NA;NA;NA;NA$",taxatab$taxon),]
  
  if(xLevel=="family" | xLevel=="genus" | xLevel=="species" )  {
    if(length(grep(";NA;NA;NA$",taxatab$taxon))>0) taxatab<-taxatab[-grep(";NA;NA;NA$",taxatab$taxon),]
  }
  
  if(xLevel=="genus" | xLevel=="species" ) {
    if(length(grep(";NA;NA$",taxatab$taxon))>0)  taxatab<-taxatab[-grep(";NA;NA$",taxatab$taxon),]
  }
  
  if(xLevel=="species") {
    if(length(grep(";NA$",taxatab$taxon))>0) taxatab<-taxatab[-grep(";NA$",taxatab$taxon),]
  }
  
  taxatab<-rm.0readtaxSam(taxatab)
  
  if(rm.trailing.NA) if(xLevel=="family") taxatab$taxon<-gsub(";NA;NA$","",taxatab$taxon)
  if(rm.trailing.NA) if(xLevel=="genus") taxatab$taxon<-gsub(";NA$","",taxatab$taxon)
  
  return(taxatab)
}

aggregate.at.xLevel<-function(taxatab,xLevel,rm.above=F,rm.trailing.NA=F){
  
  if(!xLevel %in% c("genus","family","order","class","phylum")) stop("Only allowable at genus, family, order, class or phylum level")
  
  splittaxonomy<-as.data.frame(do.call(rbind,stringr::str_split(taxatab[,1],";")))
  
  if(xLevel=="genus"){
    xPath=paste0(splittaxonomy[,1],";",splittaxonomy[,2],";",splittaxonomy[,3],";",splittaxonomy[,4],";",splittaxonomy[,5],
                 ";",splittaxonomy[,6])
    leftover<-c(";collapsed")
  }
  
  if(xLevel=="family"){
    xPath=paste0(splittaxonomy[,1],";",splittaxonomy[,2],";",splittaxonomy[,3],";",splittaxonomy[,4],";",splittaxonomy[,5])
    leftover<-c(";collapsed;collapsed")
  }
  
  if(xLevel=="order"){
    xPath=paste0(splittaxonomy[,1],";",splittaxonomy[,2],";",splittaxonomy[,3],";",splittaxonomy[,4])
    leftover<-c(";collapsed;collapsed;collapsed")
  }
  
  if(xLevel=="class"){
    xPath=paste0(splittaxonomy[,1],";",splittaxonomy[,2],";",splittaxonomy[,3])
    leftover<-c(";collapsed;collapsed;collapsed;collapsed")
  }
  
  if(xLevel=="phylum"){
    xPath=paste0(splittaxonomy[,1],";",splittaxonomy[,2])
    leftover<-c(";collapsed;collapsed;collapsed;collapsed;collapsed")
  }
  
  taxatab<-aggregate(taxatab[,-1],by = list(xPath),FUN=sum)
  
  colnames(taxatab)[1]<-"taxon"
  
  taxatab[,1]<-paste(taxatab[,1],leftover,sep = "")
  
  if(rm.above) taxatab<-keep.below.xLevel.assigns(taxatab,xLevel,rm.trailing.NA = rm.trailing.NA)
  
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

add.lineage.df<-function(dframe,ncbiTaxDir,taxCol="taxids",as.taxids=F){
  
  if(!taxCol %in% colnames(dframe)) {stop("No column called taxids")}

  #write taxids to file
  taxids_fileA<-paste0("taxids",as.numeric(Sys.time()),".txt")
  write.table(unique(dframe[,taxCol]),file = taxids_fileA,row.names = F,col.names = F,quote = F)
  
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
  dframe<-merge(dframe,lineage[,c("taxids","new_taxids","path")],by.x = taxCol,by.y = "taxids")
  dframe$old_taxids<-dframe[,taxCol]
  dframe$taxids<-dframe$new_taxids
  dframe$new_taxids=NULL
  dframe<-cbind(dframe,do.call(rbind, stringr::str_split(dframe$path,";")))
  colnames(dframe)[(length(dframe)-6):length(dframe)]<-c("K","P","C","O","F","G","S")
  if(as.taxids==F){
    dframe$K<-as.character(dframe$K)
    dframe$P<-as.character(dframe$P)
    dframe$C<-as.character(dframe$C)
    dframe$O<-as.character(dframe$O)
    dframe$F<-as.character(dframe$F)
    dframe$G<-as.character(dframe$G)
    dframe$S<-as.character(dframe$S)
  
  #change empty cells to "unknown"
  dframe[,(length(dframe)-6):length(dframe)][dframe[,(length(dframe)-6):length(dframe)]==""]<- "unknown"
  }
  
  dframe$path=NULL
  unlink(taxids_fileA)
  unlink(taxids_fileB)
  unlink(taxids_fileC)
  return(dframe)
}

get.children.taxonkit<-function(df,column="taxids",ncbiTaxDir){
  df$taxids<-df[,column]
  #if(is.null(df$taxids)) {stop("No column called taxids")}
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
#' Perform BLAST 
#' @param infastas A vector of filenames of fastas to BLAST.
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
blast.min.bas<-function(infastas,refdb,blast_exec="blastn",wait=T,taxidlimit=NULL,inverse=F,taxidname=NULL,ncbiTaxDir=NULL,overWrite=F
                        ,max_target_seqs=100,task="megablast",more=NULL){
  
  # more should be a a vector with flags (including dashes) and a single item for each flag,e.g. 
  #more=c("-word_size", 6, "-perc_identity", 50, "-qcov_hsp_perc", 90, "-gapopen", 0, "-gapextend", 2, "-reward", 1, "-penalty", -1)
  
  if(!is.null(taxidlimit)) if(is.null(ncbiTaxDir)) stop("to use taxidlimit, ncbiTaxDir must be supplied")
  if(!is.null(taxidlimit)) if(is.null(taxidname)) stop("to use taxidlimit, taxidname must be supplied")
  if(!is.null(taxidlimit)) message("Make sure infastas,taxidlimit & taxidname are in correct order")
  
  t1<-Sys.time()
  
  require(processx)
  
  if(length(infastas)==1 | length(infastas)==2 | length(infastas)==3) threads<-8
  if(length(infastas)==4 | length(infastas)==5 | length(infastas)==6) threads<-4
  if(length(infastas)==7 | length(infastas)==8 | length(infastas)==9) threads<-2
  if(length(infastas)>9) threads<-1
  
  if(overWrite==F) {
  
  continue<-data.frame("file"<-infastas,"response"="y")
  continue$response<-as.character(continue$response)
  
  for(i in 1:length(infastas)){
    if(paste0(gsub(x = infastas[i],pattern = ".fasta",replacement = ".blast.txt")) %in% list.files()){
      continue[i,2]<-readline(paste0("The following file already exists, Overwrite? (y/n):", "
                                     ",gsub(x = infastas[i],pattern = ".fasta",replacement = ".blast.txt")))
    }
  }
  
  
  
  if("n" %in% continue$response) stop("Abandoned blast due to overwrite conflict")
  
  }
  
  if(!is.null(taxidlimit)){
    
    h<-list()
    
    for(i in 1:length(infastas)){
      system2(command = "taxonkit",args = c("list", "--ids", taxidlimit[i], "--indent", '""',"--data-dir",ncbiTaxDir)
                ,wait=T,stdout = paste0(taxidname[i],"_taxidlimit.temp.txt"))
        
        #remove blank row
        taxidlist<-read.table(paste0(taxidname[i],"_taxidlimit.temp.txt"))
        write.table(taxidlist,paste0(taxidname[i],"_taxidlimit.txt"),row.names = F,quote = F,col.names = F)
        
        unlink(paste0(taxidname[i],"_taxidlimit.temp.txt"))
        message(paste("taxidlist saved to",paste0(taxidname[i],"_taxidlimit.txt")))
        
        if(inverse==F){
          h[[i]]<-process$new(command = blast_exec, 
                          args=c("-query", infastas[i], "-task", task,"-db",refdb,"-outfmt",
                                 "6 qseqid evalue staxid pident qcovs","-num_threads", threads, "-taxidlist", 
                                 paste0(taxidname[i],"_taxidlimit.txt"),"-max_target_seqs", max_target_seqs, "-max_hsps","1",more, "-out",
                                 paste0(gsub(x = infastas[i],pattern = "\\.fasta",replacement = ".blast.txt"))),echo_cmd = T,
                          stderr = paste0("blast.error.temp.processx.file",i))
        }
        if(inverse==T){
          h[[i]]<-process$new(command = blast_exec, 
                              args=c("-query", infastas[i], "-task", task,"-db",refdb,"-outfmt",
                                     "6 qseqid evalue staxid pident qcovs","-num_threads", threads, "-negative_taxids", 
                                     paste0(taxidname[i],"_taxidlimit.txt"),"-max_target_seqs", max_target_seqs, "-max_hsps","1",more, "-out",
                                     paste0(gsub(x = infastas[i],pattern = "\\.fasta",replacement = ".blast.txt"))),echo_cmd = T,
                              stderr = paste0("blast.error.temp.processx.file",i))
        }
    }
  }
  
  if(is.null(taxidlimit)){  
    
    h<-list()
    
    for(i in 1:length(infastas)){
      
      h[[i]]<-process$new(command = blast_exec,
                          args=c("-query", infastas[i], "-task", task,"-db",refdb,"-outfmt",
                                 "6 qseqid evalue staxid pident qcovs","-num_threads", threads, "-max_target_seqs", 
                                 max_target_seqs, "-max_hsps","1",more, "-out",
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
      #message(readLines(paste0("blast.error.temp.processx.file",i)))
      #unlink(paste0("blast.error.temp.processx.file",i))
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


make.blastdb.bas<-function(infasta,makeblastdb_exec="makeblastdb",addtaxidsfasta=F, ncbiTaxDir, dbversion=4,do.checks=T){
  require(processx)
  
  message("Reminder: Assumes header includes taxid=taxid;")
  message("Reminder: Assumes there are only spaces between attributes, not within them, e.g. species names should not have spaces")
  
  
  message("Reminder: Only works for blast 2.9.0: Current version:")
  system2(command = makeblastdb_exec,args = c("-version"))
  if(!infasta %in% list.files()) stop("infasta not found in current directory")

  tempfasta<-phylotools::read.fasta(infasta)
  
  tempfasta$ids<-do.call(rbind,strsplit(as.character(tempfasta$seq.name)," "))[,1]
  taxids<-do.call(rbind,strsplit(as.character(tempfasta$seq.name),"taxid="))[,2]
  tempfasta$taxids<-do.call(rbind,strsplit(as.character(taxids),";"))[,1]
  
  #Add "database name to header
  message("Adding db name to headers")
  tempfasta$db<-gsub(".fasta","",infasta)
  
  #remove any quotes
  tempfasta$ids<-gsub('"',"",tempfasta$ids)
  
  #ensure ids are <50 characters
  if(sum(nchar(tempfasta$ids)>50)>0) {
    message("The following ids exceed 49 characters and will be reduced to the first 49:")
    print(tempfasta[nchar(tempfasta$ids)>49,])
    tempfasta$ids<-stringr::str_trunc(as.character(tempfasta$ids),width = 49)
  }

  #give unique ids
  if(length(unique(tempfasta$ids))!=length(tempfasta$ids)) {
    message("The following IDs not unique - giving unique ids - ensure to check if mapping back!")
    print(tempfasta[duplicated(tempfasta$ids),])
    tempfasta$ids<-make.unique(tempfasta$ids, sep = ".")
  }
  
  #make final headers & write fasta
  tempfasta$seq.name<-paste0(tempfasta$ids," taxid=",tempfasta$taxids, "; db=",tempfasta$db)
  phylotools::dat2fasta(tempfasta,gsub(".fasta",".blastdbformatted.fasta",infasta))
  
  #make mapping file
  mappingfile<-paste0("mapping",as.numeric(Sys.time()),".txt")
  write.table(tempfasta[,c("ids","taxids")],mappingfile,quote = F,sep = " ",row.names = F,col.names = F)
  
  message("Creating blastdb")
  system2(command = makeblastdb_exec, 
          args=c("-in", gsub(".fasta",".blastdbformatted.fasta",infasta), 
                 "-dbtype", "nucl", 
                 "-blastdb_version", dbversion,
                 "-parse_seqids","-taxid_map",mappingfile,"-out",
                 gsub(".fasta","",infasta)), stderr = "",wait = T)
  
  #testblast
  if(do.checks){
  message("Running test blast")
  phylotools::dat2fasta(head(tempfasta,n=1),gsub(".fasta",".blastdbformatted.test.fasta",infasta))
  h<-blast.min.bas(gsub(".fasta",".blastdbformatted.test.fasta",infasta),refdb = gsub(".fasta","",infasta))
  check.blasts(gsub(".fasta",".blastdbformatted.test.fasta",infasta),h)
  message("Does new db have taxids - column V3?")
  print(data.table::fread(gsub(".fasta",".blastdbformatted.test.blast.txt",infasta)))
  system2(command = "blastdbcheck",args = c("-must_have_taxids","-db",gsub(".fasta","",infasta)))
  }
  
  unlink(gsub(".fasta",".blastdbformatted.test.blast.txt",infasta))
  unlink(gsub(".fasta",".blastdbformatted.test.fasta",infasta))
  unlink(mappingfile)
  
  message("Note: Final warning messages about number of columns ok")
  message("An error like:
          Testing 1 volume(s).
  /media/sf_Documents/WORK/CIBIO/temp/temp / MetaData:   [ERROR] caught exception.
 Result=FAILURE. 1 errors reported in 1 volume(s).
          Probably means you need to put the taxdb.bti and taxdb.btd files in the path or in the same folder as database.
          To add to path, do cd $BLASTDB, and put those files in there")
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
  
  require(processx)
  
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
  if(infasta==outfasta) stop("infasta and outfasta must have different names")
  require(processx)
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
  
  otutab.bins<-as.data.frame(merge(x = otutab,y = bins[,c("K","P","C","O","F","G","S","OTU_name")], by = "OTU_name"))
  otutab.bins.all<-as.data.frame(merge(x = otutab,y = bins[,c("K","P","C","O","F","G","S","OTU_name")], by = "OTU_name",all.x=T))
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
splice.taxatables<-function(files){
  
  message("Note:sample names must contain project name with dash")
  
  #read files
  taxatables<-list()
  for(i in 1:length(files)){
    taxatables[[i]]<-data.table::fread(files[i],data.table = F,sep = "\t")
  }
  
  #split by project
  projectnames<-list()
  for(i in 1:length(taxatables)){
    projectnames[[i]]<-(unique(do.call(rbind,stringr::str_split(string = colnames(taxatables[[i]][,-1]),pattern = "-"))[,1]))
  }
  
  projectnames<-unique(do.call(c,projectnames))
  
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
bas.krona.plot<-function(taxatable,KronaPath=NULL,out=NULL){
  if(is.data.frame(taxatable)) {
    a<-taxatable
  } else a<-data.table::fread(taxatable,data.table = F)
  
  b<-as.data.frame(do.call(rbind, stringr::str_split(a[,1],";")))
  colnames(b)<-c("K","P","C","O","F","G","S")
  
  if(length(colnames(a))>2) {
    a$all<-rowSums(a[,-1])
  } else {
    a$all<-a[,2]
  }
  
  d<-colnames(a[,-1])
  d<-d[c(length(d),1:(length(d)-1))]
  
  for(i in 1:length(d)){
    sample<-cbind(a[,d[i]],b)
    colnames(sample)[1]<-d[i]
    write.table(sample,row.names = F,file = paste0(d[i],".krona.txt"),quote = F,sep = "\t",col.names = F)
  }
  
  if(!is.null(KronaPath)){
    command<-KronaPath
    } else {command<- "ktImportText"}
  
  if(is.data.frame(taxatable)){
    outfile<-out
  } else {outfile<-paste0(gsub(".txt",".krona.html",taxatable))}
  
  system2(command = command,args = c(paste0(d,".krona.txt"),"-o",outfile) ,stdout = F,stderr = "",wait = T)
  
  file.remove(paste0(d,".krona.txt"))
}

taxatab.stackplot<-function(taxatab,master_sheet=NULL,column=NULL,as.percent=T,as.dxns=F,facetcol=NULL,hidelegend=F,grouping="ss_sample_id",
                            out.tab=F){
  #If column names are not ss_sample_ids using 'grouping' to specify what they are
  if(out.tab==T){
    taxa<-do.call(rbind,stringr::str_split(taxatab$taxon,";"))
    taxa<-cbind(taxa,do.call(rbind,stringr::str_split(taxa[,7]," ")))
    taxa2<-as.data.frame(substr(taxa,start = 1,stop = 3))
    taxatab$taxon<-apply(taxa2,MARGIN = 1,FUN = function(x) paste0(x[1],".",x[2],".",x[3],".",x[4],".",x[5],".",x[6],".",x[8],"_",x[9]))
    }
  
  if(as.dxns==T) taxatab<-binarise.taxatab(taxatab,t=T)
  
  long<-reshape2::melt(taxatab)
  #long<-long[long$value>0,]
  
  if(!is.null(master_sheet)) {
    message("Note: think about which plots make sense. If grouping is biomaterial, then plots using e.g. extraction method
            do not make sense, because there is more than one possibility for each biomaterial")
    
    long[,4]<-ms_ss[match(long$variable,ms_ss[,grouping]),column]
    colnames(long)[4]<-column
    
    if(!is.null(facetcol)){
      for(i in 1:length(facetcol)){
        if(is.null(ms_ss[match(long$variable,ms_ss[,grouping]),facetcol[i]])) stop("facetcol not found, check spelling, not some columns become lower case")
      
        long[,4+i]<-ms_ss[match(long$variable,ms_ss[,grouping]),facetcol[i]]
        colnames(long)[4+i]<-facetcol[i]
      }
    }
    
    long$variable<-long[,match(column,colnames(long))]
    
  } else message("No master_sheet or column specified, making default plot")
  
  #number samples in each group
  #aggregate(long$value,list(long$variable),function(x) length(x)/length(unique(long$taxon)))

  if(as.dxns) if(as.percent==F) long<-long[!duplicated(long),]
  
  a<-ggplot2::ggplot(data=long , aes(y=value, x=as.factor(variable), fill=taxon))+
    theme(legend.title = element_text(size=8), legend.text=element_text(size=8),
          axis.text.x=element_text(size=8,angle=45, hjust=1),legend.position="right",legend.direction="vertical",
          title = element_text(size=6),axis.title = element_text(size=8),panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    ggtitle(paste("x axis=",column,"; as.percent=",as.percent,"; as.dxns=",as.dxns, ";facetcol=",facetcol),subtitle = )+
    xlab("") + scale_y_continuous(labels = scales::comma) 
    
  
  if(as.percent==T) if(as.dxns==F) a<-a+ylab("Proportion of reads")
  if(as.percent==F) if(as.dxns==F) a<-a+ylab("Number of reads")
  if(as.percent==T) if(as.dxns==T) a<-a+ylab(paste0("Number of ",grouping,"s with detections"))
  if(as.percent==T) if(as.dxns==T) if(grouping=="ss_sample_id") a<-a+ylab(paste0("Number of PCRs with detections"))
  if(as.percent==T) if(as.dxns==T) if(grouping=="biomaterial") a<-a+ylab(paste0("Number of field samples with detections"))
  if(as.percent==F) if(as.dxns==T) a<-a+ylab("Number of species detected")
  
  if(length(unique(long$taxon))<29)  a<-a+scale_fill_manual(values = MyCols) 
  
  if(as.dxns==F){
   if(as.percent) a<-a+geom_bar(position="fill", stat="identity") else a<-a+geom_bar(stat = "identity")
  } else a<-a+geom_bar(stat = "identity")
  
  if(!is.null(facetcol))  if(length(facetcol)==1) a<-a+facet_wrap(vars(long[,match(facetcol[1],colnames(long))]),scales = "fixed")
  if(!is.null(facetcol))  if(length(facetcol)>1) a<-a+facet_wrap(facetcol,scales = "fixed")
  
  if(hidelegend) {
    message("Outputting as a list where first element is plot and second is legend")
    plotlist<-list()
    plotlist[[1]]<-a+theme(legend.position="none")
    plotlist[[2]]<-ggpubr::as_ggplot(cowplot::get_legend(a))
    return(plotlist)
  }else return(a)
}

#Plotting bray distance matrix PCA
taxatab.pca.plot.col<-function(taxatab,ms_ss,grouping="ss_sample_id",factors,lines=F,longnames=F,shortnames=F,ellipse=T,hidelegend=F){
  message("Assuming grouping has already been done")
  
  taxatab2<-taxatab
  taxatab2<-rm.0readtaxSam(taxatab2)
  taxatab2<-binarise.taxatab(taxatab2)
  distance_matrix<-taxatab2bray(taxatab2)
  
  cmds<-cmdscale(distance_matrix,list. = T, eig = T)
  cor.cmds<-cor(taxatab2,cmds$points)
  VarExplainedPC1<-round(cor(vegan::vegdist(cmds$points[,1],method = "euclidean"),distance_matrix)^2,digits = 2)
  VarExplainedPC2<-round(cor(vegan::vegdist(cmds$points[,2],method = "euclidean"),distance_matrix)^2,digits = 2)
  
  #create loadings for plotting
  loadings<-as.data.frame(cor.cmds)
  cmdspoints<-as.data.frame(cmds$points)
  cmdspoints[,grouping]<-rownames(cmdspoints)
  cmdspoints<-merge(cmdspoints,ms_ss,by=grouping,all.x = T)

  if(length(factors)>1) {
    cmdspoints$factorx<-do.call(paste,c(cmdspoints[,factors],sep="."))
  } else cmdspoints$factorx<-cmdspoints[,factors]
  
  #plot
  p<-ggplot(cmdspoints,aes(x=V1,y=V2))+
    geom_point(aes(size=1,colour=factorx,stroke=1))+
    xlab(bquote("Variance explained =" ~ .(VarExplainedPC1)))+
    ylab(bquote("Variance explained =" ~ .(VarExplainedPC2))) +
    theme_bw()+
    guides(size = FALSE)+
    guides(shape=guide_legend(override.aes = list(size = 4)))+
    theme(legend.title = element_blank())+
    theme(legend.text=element_text(size=12))+
    theme(legend.spacing.x = unit(0.2, 'cm'))+
    theme(axis.title = element_text(size = 12))
  
  # if(!is.null(facetcol))  if(length(facetcol)==1) p<-p+facet_wrap(vars(cmdspoints[,match(facetcol[1],colnames(cmdspoints))]),scales = "fixed")
  # if(!is.null(facetcol))  if(length(facetcol)>1) p<-p+facet_wrap(facetcol,scales = "fixed")
  
  message("Principal Coordinates Analysis plot of community simmilarity using Bray-Curtis distances")
  
  if(ellipse){
    message("Note: Ellipses will not be calculated if there are groups with too few data points")
    p<- p+stat_ellipse(aes(colour=factorx 
                           ,fill=factorx
    )
    ,type = "norm", level=0.90, 
    geom = "polygon",alpha=0.2,
    show.legend = F,segments = 100) 

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
  
  plotlist<-list()
  
  if(hidelegend) {
    plotlist[[1]]<-p+theme(legend.position = "none")
    plotlist[[2]]<-ggpubr::as_ggplot(cowplot::get_legend(p))
  } else plotlist[[1]]<-p
  
  return(plotlist)
  
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
  
  require(processx)
  
  if(length(grep("I",Pf))>0)(Pf<-gsub("I","N",Pf))
  if(length(grep("I",Pr))>0)(Pr<-gsub("I","N",Pr))
  
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
}

#modelling taxtabs
binarise.taxatab<-function(taxatab,t=F){
  #transpose to have species as columns
  taxatab2<-as.data.frame(t(taxatab[,-1,drop=F]))
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
#' @param path.to.save.to Full or relative path to download NCBI taxonomy to.
#' @return An unpacked NCBI taxonomy
#' @note Any previous versions of the taxonomy located the same directory are overwritten
#' @examples
#' getNCBItaxonomyDump(".")
#' @export
getNCBItaxonomyDump<-function(path.to.save.to){
    
  require(processx)
  
  current.path<-getwd()
  
  if(dir.exists(path.to.save.to)==F) {
    message("directory ", path.to.save.to, " does not exist. Creating directory")
    dir.create(path.to.save.to) 
  }
  
  setwd(path.to.save.to)
  
  if("taxdump.tar.gz" %in% list.files()) file.remove("taxdump.tar.gz")
  
  cb <- function(line, proc) {cat(line, "\n")}
  processx::run(command = "wget",
                args=c('ftp://ftp.ncbi.nlm.nih.gov://pub/taxonomy/taxdump.tar.gz'),echo=F,stderr_line_callback = cb)
  processx::run(command = "tar",
                args=c("-xvzf","taxdump.tar.gz"),echo=F,stderr_line_callback = cb)
  
  file.remove("taxdump.tar.gz")
  
  setwd(current.path)
  
  message("SUCCESS! Point to ", path.to.save.to, " to use this taxonomy dump with other programs")
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
MyCols <- c("dodgerblue2","darkred", # red
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
  
  if(nrow(taxatab.tf)>0){ 
  
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
      disabledTaxaDf<-data.table::fread(disabledTaxaFile, data.table = F,sep = "\t")
      if(!"taxids" %in% colnames(disabledTaxaDf)) stop("No column called 'taxids'")
      if(!"disable_species" %in% colnames(disabledTaxaDf)) stop("No column called 'disable_species'")
      if(!"disable_genus" %in% colnames(disabledTaxaDf)) stop("No column called 'disable_genus'")
      if(!"disable_family" %in% colnames(disabledTaxaDf)) stop("No column called 'disable_family'")
      
      #convert to logical
      disabled.cols<-c("disable_species","disable_genus","disable_family")
      for(i in 1:3){
        disabledTaxaDf[,disabled.cols[i]]<-as.logical(disabledTaxaDf[,disabled.cols[i]])
        disabledTaxaDf[,disabled.cols[i]][is.na(disabledTaxaDf[,disabled.cols[i]])]<-"FALSE"
      }

      disabledSpecies<-disabledTaxaDf[disabledTaxaDf$disable_species==T,"taxids"]
      disabledSpecies<-disabledSpecies[!is.na(disabledSpecies)]
      
      disabledGenus<-disabledTaxaDf[disabledTaxaDf$disable_genus==T,"taxids"]
      disabledGenus<-disabledGenus[!is.na(disabledGenus)]
      
      disabledFamily<-disabledTaxaDf[disabledTaxaDf$disable_family==T,"taxids"]
      disabledFamily<-disabledFamily[!is.na(disabledFamily)]
      
      contributordf$species_disabled<-contributordf$taxids %in% disabledSpecies
      contributordf$genus_disabled<-contributordf$taxids %in% disabledGenus
      contributordf$family_disabled<-contributordf$taxids %in% disabledFamily
      
    } else {
      contributordf$species_disabled<-"FALSE"
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
      
    } else {(message("Not adding 'May be improved' column, as no pidents provided"))}
    
    
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
    
    write.table(x = contributordf,file=gsub(".txt",".contr.txt",filtered.taxatab),
                sep="\t",quote = F,row.names = F)
  } else message("
                 ERROR: No sequences assigned to any taxon, not making contributor table
                 ")
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


stats.by.rank<-function(taxatab,grouphtf=T){
  a<-summary.dxns.by.taxon(taxatab)
  a<-a[a$taxon!="NA;NA;NA;NA;NA;NA;NA",]
  a<-a[a$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]
  
  a$rank<-bas.get.ranks(taxatab = a,grouphtf = grouphtf)
  
  d<-aggregate(a$n.samples,by=list(a$rank),FUN=sum)
  
  e<-aggregate(a$total.reads,by=list(a$rank),FUN=sum)
  
  colnames(d)<-c("rank","dxns")
  colnames(e)<-c("rank","total.reads")
  
  d$percent.dxns<-round(d$dxns/sum(d$dxns)*100,digits = 1)
  e$percent.reads<-round(e$total.reads/sum(e$total.reads)*100,digits = 1)
  
  f<-cbind(d,reads=e$total.reads,percent.reads=e$percent.reads)
  
  return(f)
}

sumreps<-function(taxatab,ms_ss,grouping="Sample_Name",discard=T,current.grouping="ss_sample_id"){
  
  require(tidyverse)
  
  if(discard) { 
    message("If only one rep, will discard the rep completely")
  } else { message("If only one rep, will keep that rep") }
  
  #get group ids
  mapping<-ms_ss[match(colnames(taxatab[,-1]),ms_ss[,current.grouping]),grouping]
  
  clean<-list()
  for(i in 1:length(unique(mapping))){
    df2<-taxatab[, c(FALSE,mapping %in% unique(mapping)[i]),drop=F]

    if(length(colnames(df2))==1) if(discard==F) {
      df3<-df2
    } else {
      df3=NULL
    }
    
    
    if(length(colnames(df2))>1) df3<-as.data.frame(rowSums(df2))
    
    if(!is.null(df3)) colnames(df3)<-unique(mapping)[i]
    clean[[i]]<-df3
  }
  
  clean<-clean %>% discard(is.null)
  
  out<-cbind(taxon=taxatab$taxon,do.call(cbind,clean))
  
  out<-rm.0readtaxSam(taxatab = out)

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
  
  
  ### should use qcov setting - No, I asm interested in 'scov'
  ### should use blastn - optionally I think, recommended for shorter fragments (<100?)
  ### should use bowtie?
  ### max targets 3 because files get huge otherwise
  ### max 
  
  system2(command = "makeblastdb", args=c("-in", refs, "-dbtype", "nucl", "-parse_seqids","-out","refdb"),wait=T)
  
  system2(command = "blastn", args=c("-query", queries.to.map, "-task", "megablast","-db","refdb",
                                     "-outfmt",'"7 qseqid qlen qstart qend slen sstart send length pident qcovs sstrand"',
                                     "-num_threads", "16","-max_target_seqs", "3"),stdout=out,wait = T)
  
}


count.nseqs.in.fams.in.fasta<-function(fasta,ncbiTaxDir){
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

count.taxa.in.fams.in.fasta<-function(fasta,ncbiTaxDir){
  message("Reminders: headers must contain \"taxid=taxid;\"")
  tempfasta<-phylotools::read.fasta(fasta)
  tempfasta$taxids<-stringr::str_match(tempfasta$seq.name, "taxid=(.*?);")[,2]
  tempfasta<-add.lineage.df(tempfasta,ncbiTaxDir)
  tempfasta$count<-1
  tempfasta$path<-paste(tempfasta$K,tempfasta$P,tempfasta$C,tempfasta$O,tempfasta$F,tempfasta$G,tempfasta$S,sep = ";")
  a<-aggregate(tempfasta$count,by=list(tempfasta$path),FUN=sum)
  a$x=1
  b<-as.data.frame(do.call(rbind,stringr::str_split(a$Group.1,";"))[,1:5])
  a$fampath<-paste(b[,1],b[,2],b[,3],b[,4],b[,5],sep = ";")
  b<-aggregate(a$x,by=list(a$fampath),FUN=sum)
  colnames(b)<-c("Family","ntaxa")
  return(b)
}

ecopcr2refs<-function(ecopcrfile,outfile,bufferecopcr,Pf,Pr,selection="genus"){
  message("Reminder: assumes -D option was used for ecopcr")
  message("Assumes clean.ecopcrouput was used first")
  #read results
  a<-data.table::fread(ecopcrfile,sep = "\t",data.table = F)
  
  message(nrow(a)," total hits")
  
  #remove seqs with edges less than buffer - why necessary?, because ecopcr with buffer includes hits where buffer extended past length of seq
  message("Removing sequences with edges less than buffer")
  #location of first base of insert
  a$start<-stringr::str_locate(a$sequence,"[A-Z]")[1] #all left edges are 25, which makes sense
  a$sequence2<-substr(x = a$sequence, start = a$start, stop = nchar(a$sequence))
  #location of last base of insert
  a$end<-stringr::str_locate(a$sequence2,"[a-z]")[,1]-1
  #noticed one or two seqs that did not have bufers for some reason, resulting in NAs
  if(sum(is.na(a$end))>0) a<-a[!is.na(a$end),]
  
  ##a$sequence3<-substr(x = a$sequence2, start = 1, stop = a$end)
  #length of left edge
  a$leftedge<-a$start-1
  #length of right edge
  a$rightedge<-nchar(a$sequence2)-a$end
  #remove seqs with left buffer less than buffer + primer  
  before<-nrow(a)
  a<-a[!a$leftedge<(bufferecopcr+nchar(Pf)),]
  after<-nrow(a)
  message(before-after," sequences removed for failing left buffer")
  #remove seqs with right buffer less than buffer
  a<-a[!a$rightedge<(bufferecopcr+nchar(Pr)),] ####################introduces NA to AC
  after2<-nrow(a)
  message(after-after2," sequences removed for failing right buffer")
  
  #run a check on if there were any buffer that were too long...shouldnt be
  if(length(table(a$rightedge))!=1) message("Error:  Some left buffers are too long...check reason")
  if(length(table(a$leftedge))!=1)  message("Error: Some left buffers are too long...check reason")
  
  #choosing one seq per genus
  if(!is.null(selection)){
    message("Choosing one sequence per ",selection)
    colnames(a)<-gsub("taxid","taxids",colnames(a))
    a$path<-paste(a$superkingdom_name,a$family_name,a$genus_name,sep = ";")
    a <- a[order(a$path,a$amplicon_length,decreasing = T),]
    a <-a[!duplicated(a[,c(paste0(selection,"_name"))]),] 
    
    message(nrow(a)," selected")
  }
  
  #output as fasta
  colnames(a)<-gsub("sequence", "seq.text",colnames(a))
  #refs should comprise buffer,binding site, insert, binding site, buffer        
          
  a$seq.text<-toupper(a$seq.text)
  a$seq.name<-paste0(a$AC," taxid=",a$taxid,"; definition=",a$definition,"; strand=",a$strand)
  phylotools::dat2fasta(a[,c("seq.name","seq.text")],outfile)
}

add.3pmms<-function(ecopcroutput,Pf,Pr){
  
  message("Adding f_mismatches_3prime: the no. of mismatches in the 3 prime half of the forward primer")
  
  #add 3' mismatches to ecopcroutput
  f_mismatch_table<-mismatch.table(ecopcroutput_clean = ecopcroutput,primer.seq = Pf,primer.direction = "f")
  f_mismatches_3prime<-as.data.frame(rowSums(f_mismatch_table[,as.integer(nchar(Pf)/2):nchar(Pf)]))
  colnames(f_mismatches_3prime)<-"f_mismatches_3prime"
  
  message("Adding r_mismatches_3prime: the no. of mismatches in the 3 prime half of the reverse primer")
  r_mismatch_table<-mismatch.table(ecopcroutput_clean = ecopcroutput,primer.seq = Pr,primer.direction = "r")
  r_mismatches_3prime<-as.data.frame(rowSums(r_mismatch_table[,as.integer(nchar(Pr)/2):nchar(Pr)]))
  colnames(r_mismatches_3prime)<-"r_mismatches_3prime"
  ecopcroutput<-cbind(ecopcroutput,f_mismatches_3prime,r_mismatches_3prime)
  
  message("Adding f_mismatches_3prime6: the no. of mismatches in the 3 prime 6bp of the forward primer")
  
  
  #add 3' mismatches to ecopcroutput - 6bp
  f_mismatches_3prime6<-as.data.frame(rowSums(f_mismatch_table[,as.integer(nchar(Pf)-5):nchar(Pf)]))
  colnames(f_mismatches_3prime6)<-"f_mismatches_3prime6"
  message("Adding r_mismatches_3prime6: the no. of mismatches in the 3 primer 6bp of the reverse primer")
  r_mismatches_3prime6<-as.data.frame(rowSums(r_mismatch_table[,as.integer(nchar(Pr)-5):nchar(Pr)]))
  colnames(r_mismatches_3prime6)<-"r_mismatches_3prime6"
  ecopcroutput<-cbind(ecopcroutput,f_mismatches_3prime6,r_mismatches_3prime6)
}

add.tm.ecopcroutput<-function(ecopcroutput){
  message("Calculating Tm for forward & reverse primer binding sites: 1) full binding site 2) 3 prime half 3) 3 prime last 6bp")
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
res.fam.or.better.basres<-function(ecopcroutput,column="res"){
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


clean.ecopcroutput<-function(ecopcrfile,Pf,Pr,min_length=NULL,max_length=NULL,rm.buffer.insert=F,buffer.used=T){
  #GENERAL CLEANING 
  
  b<-data.table::fread(file = ecopcrfile,data.table = F,
                       sep = "|",
                       col.names = c("AC","seq_length","taxid","rank","species","species_name","genus",
                                     "genus_name","family","family_name","superkingdom","superkingdom_name",
                                     "strand","forward_match","forward_mismatch","forward_tm","reverse_match",
                                     "reverse_mismatch","reverse_tm","amplicon_length","sequence","definition"))
  #remove any apparently erroneous entries
  #those that have no amplicon length or forward_mismatch or reverse_mismatch
  
  step1<-nrow(b)
  message("A total of ",step1," hits")
  
  message("removing hits that have no amplicon length or forward_mismatch or reverse_mismatch, apparent errors")
  b<-b[!is.na(b$amplicon_length),]
  b<-b[!is.na(b$forward_mismatch),]
  b<-b[!is.na(b$reverse_mismatch),]
  step2<-nrow(b)
  message(step1-step2, " hits removed")
  
  if(nrow(b)<1) stop("No hits",call. = F)

  message("Removing hits outside desired lengths, if provided")
  if(!is.null(min_length)) b<-b[!b$amplicon_length<min_length,]
  if(!is.null(max_length)) b<-b[!b$amplicon_length>max_length,]
  step3<-nrow(b)
  message(step2-step3, " hits removed")
  
  message("Removing duplicates (i.e. pick one entry per AC, based on lowest mismatches)")
  b$total_mismatches<-as.numeric(b$forward_mismatch)+as.numeric(b$reverse_mismatch)
  b <- b[order(b$AC,b$total_mismatches),]
  b<-b[!duplicated(b$AC),]
  step4<-nrow(b)
  message(step3-step4, " hits removed")
  
  message("Removing weird primer mismatches, apparent errors, only a few usually")
  b<-b[!nchar(b$forward_match,allowNA = T)<nchar(Pf),]
  b<-b[!nchar(b$reverse_match,allowNA = T)<nchar(Pr),]
  step5<-nrow(b)
  message(step4-step5, " hits removed")
  
  message("Adding full seqs")
  if(buffer.used==F){
    message("As buffer.used=F, fullseq will be primer-insert-primer")
    b$sequence2<-gsub("[a-z]","",b$sequence)
    b$fullseq<-paste0(b$forward_match,b$sequence2,b$reverse_match)
    b$sequence2<-NULL
  } else {
    message("As buffer.used=T, fullseq will be buffer-primer-insert-primer-buffer")
    b$fullseq<-b$sequence
  }
  
  message("Removing entries with Ns in full seqs")
  if(length(grep("N",ignore.case = T,x = b$fullseq))>0)  b<-b[-grep("N",ignore.case = T,x = b$fullseq),]
  step5a<-nrow(b)
  message(step5-step5a, " hits removed")
  
  message("Removing entries with ambiguities in full seqs")
  ambigs<-c("R","Y","S","W","K","M","B","D","H","V")
  for(i in 1:length(ambigs)){
    if(length(grep(ambigs[i],ignore.case = T,x = b$fullseq))>0)  b<-b[-grep(ambigs[i],ignore.case = T,x = b$fullseq),]
  }
  step5b<-nrow(b)
  message(step5a-step5b, " hits removed")
  
  if(rm.buffer.insert==T) {
    message("Removing the buffer portions of the inserts, including primers")
    b$sequence<-gsub("[a-z]","",b$sequence)
  }
  
  #write  file
  write.table(x=b,file = paste0(ecopcrfile,".clean"),quote = F,sep = "\t",row.names = F,append=F)
  
  message("Cleaned ecopcr results saved in ", paste0(ecopcrfile,".clean"))
  
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
  message("Dereplicating by family. Bias tables are based on dereplicated sequences only")
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
  #mean_fam.res.obi<-calc.fam.res(ecopcroutput,res = "obi")
  #colnames(mean_fam.res.obi)<-gsub("pc_res_fam_or_better", "res.fam.plus.obi",colnames(mean_fam.res.obi))
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
  j$scovs<-j$length/j$slen
  
  #remove duplicate hsp hits, based on highest scov
  j2 <- j[order(j$qseqid,j$scovs,decreasing = T),]
  j2<-j2[!duplicated(j2$qseqid),]
  count_hits<-length(j2$qseqid)
  
  #remove hits less than specified scov
  j3<-j2[j2$scovs>scov,]
  count_qc<-length(j2$qseqid)
  
  #remove hits less than specified pident
  j2<-j2[j2$pident>pident,]
  count_pident<-length(j2$pident)
  
  #merge query and blast results
  k<-merge(x = j2,y = n,by = "qseqid",all.y = F) 
  
  message("outputting as fasta")
  message("keeping full sequence for queries that passed buffers, rather than extracting sequence")
  k_export<-k[,c("seq.name","seq.text")]
  count_final_db<-length(k_export$seq.name)
  phylotools::dat2fasta(k_export,outfile = out)
  
  ###########################
  message(paste("Done.", "From", count_queries,"sequences,", count_hits, "mapped to a reference,",count_qc,
            "of which had >",scov*100,"% coverage,",count_pident,"of which had greater than",pident,"% percentage identity"))
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




mismatch.table<-function(ecopcroutput_clean,primer.seq,primer.direction){
  
  #split query forward primer
  if(primer.direction=="f"){
    sst = split.primer.into.cols(primer.seq = primer.seq,ecopcroutput = ecopcroutput_clean,primer.direction = "f")
  }
  
  if(primer.direction=="r"){
    sst = split.primer.into.cols(primer.seq = primer.seq,ecopcroutput = ecopcroutput_clean,primer.direction = "r")
  }
   
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
  message("Note: If primer binding site sequences contains an ambiguity, this will count as a non-match to the subject primer")
  
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
  #check if query in subject primer
  out<-list()
  for(i in 1:length(pflist)){
    out[[i]]<-sst[,i] %in% pflist[[i]]
  }
  
  out2<-as.data.frame(do.call(cbind,out))  
  
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
  message("Think this is wrong!")
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


add.stats.ecopcroutput<-function(ecopcroutput,ncbiTaxDir,Ta=NULL,Pf,Pr){
  
  #change taxid to taxids, for next function
  colnames(ecopcroutput)<-gsub("taxid","taxids",colnames(ecopcroutput))
  
  #add full amplicon, including primers
  message("Adding full amplicon to fullseq, including primers")
  ecopcroutput$fullseq<-paste0(ecopcroutput$forward_match,ecopcroutput$sequence,ecopcroutput$reverse_match)
  
  #add taxonomy using taxonkit
  message("Adding full taxonomy")
  ecopcroutput<-add.lineage.df(df = ecopcroutput,ncbiTaxDir = ncbiTaxDir,as.taxids = F)
  
  #add 3 prime mismatches to ecopcroutput 
  #adds 4 columns: 3prime mismatches for half the primer, 3prime mismatches for last 6bp only,
  #for both forward and reverse primers
  ecopcroutput<-add.3pmms(ecopcroutput,Pf,Pr) 
  
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
  message("Calculating Tm of last 6 bp of primer divided by overall Tm (%)")
  ecopcroutput$tm_fw_3p6_perc<-ecopcroutput$fTms3prime6/ecopcroutput$fTms*100
  ecopcroutput$tm_rv_3p6_perc<-ecopcroutput$rTms3prime6/ecopcroutput$rTms*100
  
  #remove extraneous columns and sort
  ecopcroutput<- subset(ecopcroutput,select = -c(rank,species,species_name,genus,genus_name,family,family_name,
                                                  superkingdom, superkingdom_name,forward_tm,reverse_tm))
  
  ecopcroutput<-ecopcroutput %>% select(AC,taxids,K,P,C,O,F,G,S,amplicon_length,everything())
  ecopcroutput<-ecopcroutput %>% select(-sequence,sequence)
}

#remove detection in less than 2 reps
rm.single.rep.dxn<-function(taxatab,ms_ss,grouping="Sample_Name",negatives=NULL){
  
  message("If a detection was identified in more than one PCR replicate, all detections are kept.
  If a detection was identified in only one PCR replicate, that detection is put to zero.
          If only one PCR replicate is available,read counts are put to zero.")
  
  reads<-sum(taxatab[,-1])
  dxns<-length( which( taxatab[,-1] > 0 ) )
  
  if(!is.null(negatives)){
  message("Reminder: Negatives is the list output by negs.stats function
          Assumes that 'PCR_negative' were not done in replicate, so are exempt")
  taxatab.exempt<-taxatab[,colnames(negatives$PCR_negative)]
  taxatab<-taxatab[,!colnames(taxatab) %in% colnames(taxatab.exempt)[-1]]
  }
  
  #get group ids
  mapping<-ms_ss[match(colnames(taxatab[,-1]),ms_ss$ss_sample_id),grouping]
  
  clean<-list()
  for(i in 1:length(unique(mapping))){
    df2<-taxatab[, c(FALSE,mapping %in% unique(mapping)[i]),drop=F]

    if(length(colnames(df2))==1) df2[,1] = 0 #if only one rep, put all to 0
    if(length(colnames(df2))>1) {
      if(length(colnames(df2))>3) message(paste("
                                                Warning: This group has more than three replicates, not usual:",unique(mapping)[i]))
      for(j in 1:nrow(df2)){
        if(sum(df2[j,1:length(colnames(df2))] > 0)<2) df2[j,1:length(colnames(df2))]<-0  #if less than two reps have reads, put all to zero 
      }
    }
    
    clean[[i]]<-df2
  }
  
  if(!is.null(negatives)){
    if(ncol(taxatab.exempt)>0) {
      out<-cbind(taxatab.exempt,do.call(cbind,clean))
    } else {out<-cbind(taxon=taxatab$taxon,do.call(cbind,clean))}
  } else {out<-cbind(taxon=taxatab$taxon,do.call(cbind,clean))}
  
  out<-rm.0readtaxSam(taxatab = out)
  
  reads2<-sum(out[,-1])
  dxns2<-length( which( out[,-1] > 0 ) )
  
  message(paste(reads-reads2,"reads and",dxns-dxns2,"detections removed, from", reads,"reads and",dxns,"dxns"))
  
  out
}

rm.taxa<-function(taxatab,taxaGrep){
  taxatab2<-taxatab
  for(i in 1:length(taxaGrep)){
    removals<-grep(taxaGrep[i],taxatab2$taxon, value = T)
    
    if(length(removals)>0) {
      message("Removing taxa:")
      print(removals)
      taxatab2<-taxatab2[! taxatab2$taxon %in% removals,]
    }
  }
  taxatab2
}

taxatab.sumStats<-function(taxatab,stepname="this_step"){
  
  sumStats.list<-list()
  
  taxatab1<-summary.dxns.by.taxon(taxatab)
  
  taxatab1$rank<-bas.get.ranks(taxatab1)
  
  dfsum<-data.frame(detections=sum(taxatab1$n.samples),reads=sum(taxatab1$total.reads),
                    taxa=length(taxatab1$taxon),samples=length(colnames(taxatab[,-1])))
  
  taxa<-taxatab1[,c("taxon","total.reads","rank")]
  
  sumStats.list[[1]]<-dfsum
  sumStats.list[[2]]<-taxa
  
  names(sumStats.list)<-c(paste0(stepname,"_counts"),paste0(stepname,"_taxa"))
  
  sumStats.list
}
  
#plot negs vs real for each taxon
plot.negs.vs.real<-function(taxatab,ms_ss,real){
  require(tidyverse)
  
  long<-reshape2::melt(taxatab)
  long<-long[long$value>0,,drop=F]
  
  long$sample_type<-ms_ss[match(long$variable,ms_ss$ss_sample_id),"sample_type"]
  long$taxon<-as.character(long$taxon)
  
  splitdf<-split(long, f = long$taxon)
  
  #only include dfs with at least one negs and real
  
  ##neg sample types
  negtypes<-unique(ms_ss$sample_type)[!unique(ms_ss$sample_type) %in% real]
  
  splitdf2<-list()
  
  ##dfs with negs
  for(i in 1:length(splitdf)){
    if(!TRUE %in% (negtypes %in% splitdf[[i]]$sample_type)) {
      splitdf2[[i]]<-NULL
    } else {
      splitdf2[[i]]<-splitdf[[i]]
    }
  }
  
  splitdf2<-splitdf2 %>% discard(is.null)
  
  splitdf3<-list()
  
  ##dfs with real
  for(i in 1:length(splitdf2)){
    if(!TRUE %in% (real %in% splitdf2[[i]]$sample_type)) {
      splitdf3[[i]]<-NULL
    } else {
      splitdf3[[i]]<-splitdf2[[i]]
    }
  }
  
  splitdf3<-splitdf3 %>% discard(is.null)
  
  plotlist<-list()
  
  if(length(splitdf3)!=0) {
  
    for(i in 1:length(splitdf3)){
    
      df<-splitdf3[[i]]  
    
      df$sum<-sum(df$value)
    
      df$pc<-round(df$value/df$sum*100,digits = 1)
    
      plotlist[[i]]<-ggplot(data = df,aes(x=variable, y=value,color=sample_type)) + 
      geom_bar(fill="white", alpha=0.5, stat="identity") +
      theme(axis.text.x=element_text(size=8,angle=45, hjust=1),legend.title = element_blank()) +
      geom_text(aes(label=pc), position=position_dodge(width=0.9), vjust=-0.25) +
      ggtitle(unique(df$taxon)) +
      ylab("reads") +
      xlab("labels indicate % of reads in each sample of total reads for taxon")
    } 
  } else message("No taxa appearing in both negatives and positives")
   
  plotlist
}

qplot.taxatab<-function(taxatab,rm.nohits=T,rm.NA=F,specReadCount=NULL,qs=c(0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.99)){
  
  require(dplyr)
  
  if(!is.null(specReadCount)) if(specReadCount==0) {
    message("Ignoring specReadCount as it was set to zero") 
    specReadCount=NULL
  }
  
  
  if(rm.nohits) {
    message("Not counting no_hit detections")
    taxatab<-taxatab[taxatab$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]
  }
  if(rm.NA) {
    message("Not counting NA detections")
    taxatab<-taxatab[taxatab$taxon!="NA;NA;NA;NA;NA;NA;NA",]
  }
  
  long<-reshape2::melt(taxatab)
  long<-long[long$value!=0,]
  colnames(long)<-c("taxon","sample","reads")
  
  quants<-data.frame(quants=as.factor(qs*100),qdata=0,qcount=0)
  
  for(i in 1:length(qs)) {
    quants$qdata[i]<-round(quantile(long$reads,qs[i]),digits = 0)
    quants$qcount[i]<-nrow(long[long$reads<=quants$qdata[i],])
  }
  
  quants$label<-paste(quants$quants,quants$qdata,quants$qcount,sep = "_")
  
  a<-ggplot(data = long,aes(x=reads))+geom_histogram(bins = 100)+scale_x_log10()+xlab("reads/detection")+ylab("No. of detections")+
    geom_vline(data = quants,aes(xintercept = qdata),colour = "red", linetype = "dashed")
  
  b<-a+geom_line(data = quants,aes(x = qdata,y = (qcount/nrow(long))*max(ggplot_build(a)$data[[1]]$ymax)),colour="green")
    
  b<-b+geom_text(data=quants, mapping=aes(x=qdata, y=max(ggplot_build(a)$data[[1]]$ymax), label=label), size=4, angle=90, vjust=-0.6, hjust=1,colour="red")+
    ggtitle("Histogram of number of reads vs detections",subtitle = "Labels indicate %Dxns_readCount_NDxns. Green line = cumulative %Dxns (scale is 0-100%)")
  
  if(!is.null(specReadCount)) {
    specN<-data.frame(quants=0)
    specN$qcount<-nrow(long[long$reads<=specReadCount,])
    specN$quants<-round(specN$qcount/nrow(long)*100,digits=0)
    specN$qdata<-specReadCount
    specN$label<-paste(specN$quants,specN$qdata,specN$qcount,sep = "_")
    
    b<-b+geom_vline(data = specN,aes(xintercept = qdata),colour = "blue", linetype = "dashed")+
     geom_text(data=specN, mapping=aes(x=qdata, y=max(ggplot_build(a)$data[[1]]$ymax), label=label), size=4, angle=90, vjust=-0.6, hjust=1,colour="blue")
  }
  
  return(b)
  
}




#plot taxa occurring in both negatives and positives for each negative
plot.negs.vs.taxa<-function(taxatab,ms_ss,real){
  
  require(tidyverse)
  
  long<-reshape2::melt(taxatab)
  long<-long[long$value>0,,drop=F]
  
  long$sample_type<-ms_ss[match(long$variable,ms_ss$ss_sample_id),"sample_type"]
  long$taxon<-as.character(long$taxon)
  
  splitdf<-split(long, f = long$taxon)
  
  #only include dfs with at least one negs and real
  
  ##neg sample types
  negtypes<-unique(ms_ss$sample_type)[!unique(ms_ss$sample_type) %in% real]
  
  splitdf2<-list()
  
  ##dfs with negs
  for(i in 1:length(splitdf)){
    if(!TRUE %in% (negtypes %in% splitdf[[i]]$sample_type)) {
      splitdf2[[i]]<-NULL
    } else {
      splitdf2[[i]]<-splitdf[[i]]
    }
  }
  
  splitdf2<-splitdf2 %>% discard(is.null)
  
  splitdf3<-list()
  
  ##dfs with real
  for(i in 1:length(splitdf2)){
    if(!TRUE %in% (real %in% splitdf2[[i]]$sample_type)) {
      splitdf3[[i]]<-NULL
    } else {
      splitdf3[[i]]<-splitdf2[[i]]
    }
  }
  
  splitdf3<-splitdf3 %>% discard(is.null)
  
  if(length(splitdf3)!=0) {
  
    #get negatives that had a taxon in themselves and real samples
    splitdf4<-list()
    for(i in 1:length(splitdf3)){
      splitdf4[[i]]<-splitdf3[[i]][splitdf3[[i]]$sample_type %in% neg.types,c("variable","sample_type")]
    }
    final.negs<-unique(do.call(rbind,splitdf4))
    
    #get problem taxa (for %s later)
    splitdf3taxa<-list()
    for(i in 1:length(splitdf3)){
      splitdf3taxa[[i]]<-splitdf3[[i]][1,1]
    }
    
    final.taxa<-do.call(rbind,splitdf3taxa)
    
    #make data frame of problem samples
    
    finaltable<-taxatab[,c("taxon",as.character(final.negs[,"variable"]))]
    finaltable<-rm.0readtaxSam(taxatab = finaltable)
    
    long<-reshape2::melt(finaltable)
    long<-long[long$value>0,,drop=F]
    
    long$sample_type<-ms_ss[match(long$variable,ms_ss$ss_sample_id),"sample_type"]
    long$taxon<-as.character(long$taxon)
    
    splitdf<-split(long, f = long$variable)
    
    plotlist<-list()
    for(i in 1:length(splitdf)){
      df<-splitdf[[i]]
      df$sum<-sum(df$value)
      df$pc<-round(df$value/df$sum*100,digits = 1)
      
      plotlist[[i]]<-ggplot(data = df,aes(x=taxon, y=value,color=sample_type)) + 
        geom_bar(fill="white", alpha=0.5, stat="identity") +
        theme(axis.text.x=element_text(size=8,angle=45, hjust=1),legend.title = element_blank()) +
        geom_text(aes(label=pc), position=position_dodge(width=0.9), vjust=-0.25) +
        ggtitle(unique(df$variable)) +
        ylab("reads") +
        xlab("labels indicate % of reads in each taxon of total reads for sample")
    }
  } else message("No taxa appearing in both negatives and positives")
  
  plotlist
}

google.overlord<-function(url,for.MBC=F,for.post.illscript2=T,tokenDir=NULL,email=NULL,use.lengths=F,subsetList=NULL){
  
  if(sum(for.MBC,for.post.illscript2)!=1) stop("Only one output type must be TRUE")
  
  # Download master sheet
  # a<-getwd()
  # setwd(tokenDir)
  # googlesheets4::sheets_auth(email)
  # setwd(a)
  # 
  master<-google.read.master.url(url)
  
  #fix names (lowercase mainly)
  master<-fix.google.colnames(master)
  
  #remove downloaded dups
  if(length(grep("\\.\\.\\.",colnames(master),value = T))>0) {
    master<-remove.google.dups(master)
  }
  
  #subset
  if(!is.null(subsetList)){
    message("subsetting mastersheet")
    master<-subset_mastersheet(master_sheet = master,subsetList)
    message("checking subset, including sample type by default")
    #check
    print(master_xtabs(master,c("sample_type",names(subsetList))))
  }
  
  
  headers<-c("Sample_ID",	"Sample_Name",	"Sample_Plate",	"Sample_Well",	"I7_Index_ID",	"index",
             "I5_Index_ID",	"index2",	"Sample_Project",	"Description",	"index_combination","Primer_F",	"Primer_R",	"Filenames",
             "min_length","max_length","ss_sample_id",
             
             "Primer_set","experiment_id") #only required by illumina script
  
  #format for MBC option
  if(for.MBC==T){
    headers<-headers[1:(length(headers)-2)]
    
    if(FALSE %in% (headers %in% colnames(master))) stop ("Missing: ", headers[!headers %in% colnames(master)])
    
    message("Defaulting to putting all index_combinations to 'used'")
    master$index_combination<-"used"
    
    message("Defaulting to no MIDs used")
    master$MID_F<-""
    master$MID_R<-""
    
    #make ss_sample_id Sample_Name
    master$Sample_Name<-master$ss_sample_id
    
    #use only the required headers (exclude ss_sample_id)
    if(use.lengths==T){
     headers.MBC<-c("Sample_ID",	"Sample_Name",	"Sample_Plate",	"Sample_Well",	"I7_Index_ID",	"index",
               "I5_Index_ID",	"index2",	"Sample_Project",	"Description",	"index_combination","MID_F","MID_R","Primer_F",	"Primer_R",	"Filenames",
               "min_length","max_length","ss_sample_id")
    }else{
      headers.MBC<-c("Sample_ID",	"Sample_Name",	"Sample_Plate",	"Sample_Well",	"I7_Index_ID",	"index",
                     "I5_Index_ID",	"index2",	"Sample_Project",	"Description",	"index_combination","MID_F","MID_R","Primer_F",	"Primer_R",	"Filenames",
                     "ss_sample_id")
    }
    
    master<-master[,headers.MBC[-length(headers.MBC)]]
    
    master$extra_information<-""
    
    message("Putting NAs to empty")
    
    master[is.na(master)]=""
    
  }
  
  if(for.post.illscript2==T){
    
    headers<-c("primer_set","Primer_F","Primer_R","min_length","max_length","ss_sample_id","experiment_id","biomaterial","Sample_Name")
    
    if(FALSE %in% (headers %in% colnames(master))) stop ("Missing: ", headers[!headers %in% colnames(master)])
    
    message("Checks only include the mostly used headers, please check to see all desired headers exist")
  }
  
  return(master)
  
}

#changing sample_type too to match functions
fix.google.colnames<-function(master_sheet){
  colnames(master_sheet)<-gsub("Sample_Type","sample_type",colnames(master_sheet))
  colnames(master_sheet)<-gsub("Primer_set","primer_set",colnames(master_sheet))
  colnames(master_sheet)<-gsub("Index_combination","index_combination",colnames(master_sheet))
  colnames(master_sheet)<-gsub("Min_length","min_length",colnames(master_sheet))
  colnames(master_sheet)<-gsub("Max_length","max_length",colnames(master_sheet))
  master_sheet
}

google.read.master.url<-function(sheeturl,out=NULL,ws="Master_Samplesheet"){
  url2<-stringr::str_split(sheeturl,"/d/")[[1]][2]
  url2<-stringr::str_split(url2,"/")[[1]][1]
  ss_info<-googlesheets4::gs4_get(ss = url2)
  if(ws == "ENA_sample_data") {
    ss_data<-googlesheets4::read_sheet(ss = url2,sheet = ws,col_types = "c") 
    colnames(ss_data)<-ss_data[2,]
    ss_data<-ss_data[3:length(ss_data$sample_alias),]
    ss_data<-as.data.frame(ss_data[!is.na(ss_data$sample_alias),])
  } else{
    if(ws == "ENA_library_data" | ws == "ENA_library_data_Single_End") { 
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

remove.google.dups<-function(master_sheet){
  if(length(grep("\\.\\.\\.",colnames(master_sheet),value = T))>0) {
    message("Removing duplicated columns")
    a<-grep("\\.\\.\\.",colnames(master_sheet),value = T)
    b<-as.data.frame(do.call(rbind,stringr::str_split(a,"\\.\\.\\.")))
    d<-b[duplicated(b[,1]),,drop=F]
    e<-paste0(d[,1],"...",d[,2])
    master_sheet<-master_sheet[,!colnames(master_sheet) %in% e]
    colnames(master_sheet)<-gsub("\\.\\.\\..*","",colnames(master_sheet))
  }
  master_sheet
}

## inspect negatives, print plots for only those taxa detected in both negatives and real samples
full.negative.inspection<-function(taxatab,ms_ss,real){
  
  message("Checking negs")
  negatives<-negs.stats(taxatab=taxatab,ms_ss = ms_ss,real=real,ex_hominidae=F)
  
  if(is.null(negatives)) {message("No negatives with reads, skipping negatives report")
    } else {
  
      message("Making quantile plots to help deciding dxn_filter threshold")
      for(i in 1:length(negatives) ) {
        df<-negatives[[i]]
        if(length(colnames(df))>0){
          title<-paste0(names(negatives)[i],"; n_dxns=",sum(df[,-1]>0), "; n_samples=",length(colnames(df))-1)
          plot1<-qplot.taxatab(taxatab = df)
          plot1<-plot1 + ggtitle(label = title)
          print(plot1)
        }
      }
      
      message("Making barplots for only those taxa detected in both negatives and real samples - to help decide taxon filter and rm.contaminants filter")
      
      plotlist1<-plot.negs.vs.real(taxatab = taxatab,ms_ss = ms_ss,real = real)
      
      if(length(plotlist1)>0){
        for(i in 1:length(plotlist1)){
          print(plotlist1[[i]])
        }
      }
      
      message("Making barplots for only those taxa detected in both negatives and real samples - to help decide sample filter")
      plotlist2<-plot.negs.vs.taxa(taxatab = taxatab,ms_ss = ms_ss,real = real)
      
      if(length(plotlist2)>0){
        for(i in 1:length(plotlist2)){
          print(plotlist2[[i]])
        }
      }
    }
}

split.primer.into.cols<-function(primer.seq,ecopcroutput,primer.direction){
  pflist<-list()
  for(i in 1:nchar(primer.seq)){
    pflist[[i]]<-substr(primer.seq,start = i,stop = i)
  }
  
  splitprimerList<-list()
  
  if(primer.direction=="f") a <- ecopcroutput$forward_match
  if(primer.direction=="r") a <- ecopcroutput$reverse_match
  
  for(i in 1:nchar(primer.seq)){
    splitprimerList[[i]]<-substr(a, start=i, stop=i)
    names(splitprimerList)[i]<-substr(primer.seq,start = i,stop = i)
  }
  
  splitprimer<-as.data.frame(do.call(cbind,splitprimerList))
  colnames(splitprimer)<-names(splitprimerList)
  splitprimer
}

#make an experiment sheet for MBC. 
#subsetList is usually the experiment_id e.g. list(experiment_id=c("2020_02"))
# urls<-c("https://docs.google.com/spreadsheets/d/1SBBTKa6ZlbYQIWp5Xz2PrhA2v9fYXXlzpWDdn72g3lg/edit#gid=0",
#         "https://docs.google.com/spreadsheets/d/1VaMRnezBhMCul8Xk-yGAxYZ9tV7inn1QoT5doOLns4s/edit#gid=0")
#out="/mnt/Disk1/BASTIAN_RAW/2020_02/TOAD_SCORPION_mastersheet.txt"
google.make.MBC.expSheet<-function(urls,subsetList=NULL,out=NULL){
  
  masterList<-list()
  for(i in 1:length(urls)){
    masterList[[i]]<-google.overlord(url = urls[i]
                                     ,for.MBC = T,for.post.illscript2 = F,subsetList = subsetList,use.lengths = T)
  }
  
  #combining sheets:
  masterC<-do.call(rbind,masterList)
  
  if(!is.null(out)) {
    write.table(masterC,out,append = F,quote = F,sep = "\t",row.names = F)
    message("experiment sheet saved to ",out)
  }
  
  return(masterC)
}


add.target.positions<-function(cattedDLS.checked,refs,out){

  a<-phylotools::read.fasta(refs)
  ref.ids<-do.call(rbind,stringr::str_split(a$seq.name," "))[,1]
  
  write.table(ref.ids,quote = F,row.names = F,append = F,file = "ref.ids.tmp")
  
  system2(command = "seqkit", args=c("grep", "-f","ref.ids.tmp",cattedDLS.checked),stdout = out,wait = T)
  
  unlink("ref.ids.tmp")
  
  #find each seq in original
  b<-phylotools::read.fasta(out)
  
  b2<-b[order(b$seq.name),]
  a2<-a[order(a$seq.name),]
  colnames(a2)<-c("target.seq.name","target.seq.text")
  
  a3<-cbind(a2,b2)
  a3$target.seq.text<-as.character(a3$target.seq.text)
  a3$seq.text<-as.character(a3$seq.text)
  
  a3[grepl("strand=R",a3$target.seq.name),"target.seq.text"]<-insect::rc(z=a3[grepl("strand=R",a3$target.seq.name),"target.seq.text"])
  
  a3$start<-stringr::str_locate(a3$seq.text,a3$target.seq.text)[,1]
  a3$end<-stringr::str_locate(a3$seq.text,a3$target.seq.text)[,2]
  
  a3$seq.name<-paste0(a3$seq.name," target.start=",a3$start,"; target.end=",a3$end,";")
  
  phylotools::dat2fasta(dat = a3[,c("seq.name","seq.text")],outfile = out)
}


plot.ecopcroutput.bubble<-function(ecopcrdf,level="F",parent="C",parent2=NULL,plotvar="amplicon_length"){
  #parent is the facet plot. e.g. put to somethng above family if level="F", parent2 is above parent
  
  require(dplyr)
  require(ggplot2)
  
  if(!plotvar %in% c("amplicon_length","n_mismatches","res")) stop("Please choose a plotting variable: 'amplicon_length', 'res' or 'n_mismatches'")
  
  heir<-c("K","P","C","O","F","G","S")
  if(!level %in% heir)  stop("level must be one of 'K','P','C','O','F','G','S'")
  
  if(is.null(parent2)) parent2<-level
  if(is.null(parent)) parent<-level
  
  if(plotvar %in% c("amplicon_length","res")){
    ecopcr.freq <- ecopcrdf %>%
      dplyr::group_by(parent2=ecopcrdf[,parent2],parent=ecopcrdf[,parent],level=ecopcrdf[,level], plotvar=ecopcrdf[,plotvar],total_mismatches) %>%
      dplyr::summarise(freq = length(level))
    
    ecopcr.freq<-as.data.frame(ecopcr.freq)
    
    parrange <- ggplot(data = ecopcr.freq, aes(x=level, y=plotvar, size=freq,color=ecopcr.freq$total_mismatches+0.001))
    
    if(plotvar %in% "amplicon_length") {
      min <- min(ecopcr.freq$plotvar)
      max <- max(ecopcr.freq$plotvar)
      parrange<-parrange+scale_y_continuous(limits = c(min-2, max+2), breaks = seq(min-2, max+2, by = 20))

    }
    
    scales="free_x"
  }
  
  if(plotvar=="n_mismatches"){
    ecopcr.freq <- ecopcrdf %>%
      dplyr::group_by(parent=ecopcrdf[,parent],level=ecopcrdf[,level], forward_mismatch, reverse_mismatch,total_mismatches) %>%
      dplyr::summarise(freq = length(level))
    
    ecopcr.freq<-as.data.frame(ecopcr.freq)
    
    parrange <- ggplot(data = ecopcr.freq, aes(x=forward_mismatch, y=reverse_mismatch, size=freq,color=ecopcr.freq$total_mismatches+0.001))
    
    scales="fixed"
  }
    
   parrange2<- parrange  + 
      geom_point(alpha=0.8) + 
      #theme_grey() + 
      theme(legend.position = "none", strip.text = element_text(size=8), axis.title = element_text(size=10), 
            axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1, size = 7), title = element_text(face = "bold", size=15)) +
     scale_colour_gradientn(colours = heat.colors(11,rev = T))
   
   if(parent==level) parent=NULL
   if(parent2==level) parent2=NULL
   
   if(!is.null(parent)){
     if(!parent %in% heir)  stop("parent must be one of 'K','P','C','O','F','G','S'")
     if(!is.null(parent2)){
       if(!parent2 %in% heir)  stop("parent must be one of 'K','P','C','O','F','G','S'")
       parrange3<-parrange2+facet_wrap(ecopcr.freq$parent2~ecopcr.freq$parent,scales = scales,nrow = 4) 
     } else parrange3<-parrange2+facet_wrap(~ecopcr.freq$parent,scales = scales,nrow = 4)
   } else parrange3<-parrange2
   
  return(parrange3)
}

bas.get.ranks<-function(taxatab,grouphtf=T){
  
  taxa<-taxatab[,1,drop=F]  
  
  taxa$kin.res<-TRUE
  taxa$phy.res<-TRUE
  taxa$cla.res<-TRUE
  taxa$ord.res<-TRUE
  taxa$fam.res<-TRUE
  taxa$gen.res<-TRUE
  taxa$sp.res<-TRUE
  
  #change NA and collapsed to unknown
  taxa$taxon<-gsub("NA","unknown",taxa$taxon)
  taxa$taxon<-gsub("collapsed","unknown",taxa$taxon)

  if(grouphtf){
    
  if(length(grep(";unknown;unknown;unknown$",taxa$taxon))>0) taxa$fam.res<-!taxa$taxon %in% taxa$taxon[grep(";unknown;unknown;unknown$",taxa$taxon)]
  if(length(grep(";unknown;unknown$",taxa$taxon))>0) taxa$gen.res<- !taxa$taxon %in% taxa$taxon[grep(";unknown;unknown$",taxa$taxon)]
  if(length(grep(";unknown$",taxa$taxon))>0) taxa$sp.res<- !taxa$taxon %in% taxa$taxon[grep(";unknown$",taxa$taxon)]
  
  taxa$res<-"htf" #higher than family
  taxa$res[taxa$fam.res==T]<-"family"
  taxa$res[taxa$gen.res==T]<-"genus"
  taxa$res[taxa$sp.res==T]<-"species"
  } else{
    if(length(grep("unknown;unknown;unknown;unknown;unknown;unknown;unknown$",taxa$taxon))>0) {
      taxa$kin.res<-!taxa$taxon %in% taxa$taxon[grep("unknown;unknown;unknown;unknown;unknown;unknown;unknown$",taxa$taxon)]}
    if(length(grep(";unknown;unknown;unknown;unknown;unknown;unknown$",taxa$taxon))>0) {
      taxa$phy.res<-!taxa$taxon %in% taxa$taxon[grep(";unknown;unknown;unknown;unknown;unknown;unknown$",taxa$taxon)]}
    if(length(grep(";unknown;unknown;unknown;unknown;unknown$",taxa$taxon))>0) {
      taxa$cla.res<-!taxa$taxon %in% taxa$taxon[grep(";unknown;unknown;unknown;unknown;unknown$",taxa$taxon)]}
    if(length(grep(";unknown;unknown;unknown;unknown$",taxa$taxon))>0) {
      taxa$ord.res<-!taxa$taxon %in% taxa$taxon[grep(";unknown;unknown;unknown;unknown$",taxa$taxon)]}
    
    if(length(grep(";unknown;unknown;unknown$",taxa$taxon))>0) taxa$fam.res<-!taxa$taxon %in% taxa$taxon[grep(";unknown;unknown;unknown$",taxa$taxon)]
    if(length(grep(";unknown;unknown$",taxa$taxon))>0) taxa$gen.res<- !taxa$taxon %in% taxa$taxon[grep(";unknown;unknown$",taxa$taxon)]
    if(length(grep(";unknown$",taxa$taxon))>0) taxa$sp.res<- !taxa$taxon %in% taxa$taxon[grep(";unknown$",taxa$taxon)]
    
    taxa$res<-"above_kingdom" #higher than family
    taxa$res[taxa$kin.res==T]<-"kingdom"
    taxa$res[taxa$phy.res==T]<-"phylum"
    taxa$res[taxa$cla.res==T]<-"class"
    taxa$res[taxa$ord.res==T]<-"order"
    taxa$res[taxa$fam.res==T]<-"family"
    taxa$res[taxa$gen.res==T]<-"genus"
    taxa$res[taxa$sp.res==T]<-"species"
  }
  
  return(taxa$res)
  
}

bas.plot.specaccum<-function(taxatab,xlabel){
  require(vegan)
  
  message("removing no_hits and NAs")
  taxatab<-taxatab[taxatab$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]
  taxatab<-taxatab[taxatab$taxon!="NA;NA;NA;NA;NA;NA;NA",]
  
  #transpose
  ttaxatab<-taxatab
  rownames(ttaxatab)<-ttaxatab$taxon  
  ttaxatab$taxon<-NULL  
  ttaxatab<-as.data.frame(t(as.matrix(ttaxatab)))
  
  #build the species accumulation curve 
  taxatab.specaccum <- vegan::specaccum(ttaxatab,method = "exact")
  
  #plot the curve with some predefined settings
  plot(taxatab.specaccum,ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",xlab = xlabel,ylab = "number of taxa")
}

bas.plot.specrich<-function(taxatab){
  require(vegan)
  #transpose
  ttaxatab<-taxatab
  rownames(ttaxatab)<-ttaxatab$taxon  
  ttaxatab$taxon<-NULL  
  ttaxatab<-as.data.frame(t(as.matrix(ttaxatab)))
  
  pool.taxatab <- vegan::poolaccum(ttaxatab)
  plot(pool.taxatab)
}

#blast loop function

threshold.blast<-function(ecopcr.clean.df,ncbiTaxDir,out,
                          TaxlevelTest,blast_exec,makeblastdb_exec,task="megablast",all.sp.in.db=F){
  
  message("ecopcr.clean.df must have K,P,C,O,F,G,S taxonomy and 'taxids' column, as added by add.lineage.df")
  
  t1<-Sys.time()
  
  #function will overwrite 'out'
  
  #TaxLevelTest is the taxonomic level to test at
  #--if "F" then seq_i, belonging to genus_A will be blasted against everything excluding other seqs belonging to genus_A
  
  if(!TaxlevelTest %in% c("F","G","S")) stop("thresholdlevel must be one of F,G,S")
  if(TaxlevelTest=="F") ex.seqid.group<-"G"
  if(TaxlevelTest=="G") ex.seqid.group<-"S"
  if(TaxlevelTest=="S") ex.seqid.group<-NULL 
  
  a2<-ecopcr.clean.df
  a2$seq.name<-paste0(a2$AC," ","taxid=",a2$taxids,";")
  a2$seq.text<-a2$sequence
  
  message("Creating blastDB")
  phylotools::dat2fasta(a2[,c("seq.name","seq.text")],"temp.fasta")
  
  make.blastdb.bas(infasta = "temp.fasta",makeblastdb_exec = makeblastdb_exec,
                   addtaxidsfasta = F,ncbiTaxDir,do.checks=T,dbversion = 4)
  
  unlink(out)
  
  if(TaxlevelTest=="S" | all.sp.in.db==T){
    #blast
    message("straight BLAST running...")
    system2(command = blast_exec,
            args=c("-query", "temp.fasta", "-task", task,"-db","temp","-outfmt",
                   "'6 qseqid sseqid evalue staxid pident qcovs'","-evalue","100", "-max_target_seqs", 
                   "100", "-max_hsps","1","-word_size", "6","-perc_identity", "10","-qcov_hsp_perc","70",
                   "-gapopen", "0", "-gapextend", "2", "-reward", "1", "-penalty", "-1","-out", out), wait = T)
  } else {
  
    #loop blast
    for(i in 1:nrow(a2)){
      
      message("loop ",i)
      
      #make query fasta
      b<-a2[a2$AC==a2$AC[i],]
      phylotools::dat2fasta(b[,c("seq.name","seq.text")],"temp.seq.fasta")
      
      #get seqids of query group
      ex.seqids<-a2[a2[,ex.seqid.group]==b[,ex.seqid.group],"AC"]  
      write.table(ex.seqids,file = "ex.seqids.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = F)
      
      #blast
      system2(command = blast_exec,
              args=c("-query", "temp.seq.fasta", "-task", task,"-db","temp","-outfmt",
                     "'6 qseqid sseqid evalue staxid pident qcovs'","-evalue","100","-num_threads", "8", "-max_target_seqs", 
                     "100", "-max_hsps","1","-word_size", "6","-perc_identity", "10","-qcov_hsp_perc","70",
                     "-gapopen", "0", "-gapextend", "2", "-reward", "1", "-penalty", "-1","-negative_seqidlist", "ex.seqids.txt", "-out",
                     "temp.seq.blast.txt"), wait = T)
      
      #store blast results
      blastResults<-data.table::fread("temp.seq.blast.txt",sep = "\t",data.table = F)
      
      write.table(blastResults,file = out,append = T,quote = F,row.names = F,sep = "\t",col.names = F)
    }
  }
  
  t2<-Sys.time()
  
  message("Loop blast for TaxlevelTest=", TaxlevelTest," complete in ", round(difftime(t2,t1,units = "mins"),digits = 2)," mins")
  message("Output saved in ",out)
  
  unlink("temp.blastdbformatted.fasta")
  unlink("temp.n*")
  unlink("temp.seq.blast.txt")
  unlink("temp.seq.fasta")
  unlink("ex.seqids.txt")
  
  message("Warnings usually ok, just happens when reading empty blast results")
}

path.at.level<-function(pathvector,level){
  
  require(stringr)
  if(!level %in% c("K","P","C","O","F","G","S")) stop("level must be one of K,P,C,O,F,G,S")
  
  mat<-do.call(rbind,stringr::str_split(pathvector,";"))
  
  if(level=="K") newpath<-mat[,1]
  if(level=="P") newpath<-paste(mat[,1],mat[,2],sep = ";")
  if(level=="C") newpath<-paste(mat[,1],mat[,2],mat[,3],sep = ";")
  if(level=="O") newpath<-paste(mat[,1],mat[,2],mat[,3],mat[,4],sep = ";")
  if(level=="F") newpath<-paste(mat[,1],mat[,2],mat[,3],mat[,4],mat[,5],sep = ";")
  if(level=="G") newpath<-paste(mat[,1],mat[,2],mat[,3],mat[,4],mat[,5],mat[,6],sep = ";")
  if(level=="S") newpath<-paste(mat[,1],mat[,2],mat[,3],mat[,4],mat[,5],mat[,6],mat[,7],sep = ";")
  
  return(newpath)
}

# filtering and binning loop
threshold.bin.blast<-function(blastfile,ecopcr.clean.df,headers = "qseqid sseqid evalue staxid pident qcovs",
                              ncbiTaxDir = ncbiTaxDir,max_evalue = 0.001,min_qcovs = 70,top = c(1,5,10),
                              TaxlevelTest="F", pidents=c(90,80,70)){
  
  t1<-Sys.time()
  
  message(blastfile)
  
  message("ecopcr.clean.df must have K,P,C,O,F,G,S taxonomy and 'taxids' column, as added by add.lineage.df")
  
  if(!TaxlevelTest %in% c("F","G","S")) stop("thresholdlevel must be one of F,G,S")
  
  #start making final table
  ecopcr.clean.df$path<-paste(ecopcr.clean.df$K,ecopcr.clean.df$P,ecopcr.clean.df$C,ecopcr.clean.df$O,ecopcr.clean.df$F,
                              ecopcr.clean.df$G,ecopcr.clean.df$S,sep = ";")
  
  if(TaxlevelTest=="F") ecopcr.clean.df$origpath<-path.at.level(ecopcr.clean.df$path,level = "F")
  if(TaxlevelTest=="G") ecopcr.clean.df$origpath<-path.at.level(ecopcr.clean.df$path,level = "G")
  if(TaxlevelTest=="S") ecopcr.clean.df$origpath<-ecopcr.clean.df$path
  
  #filter blast results at each top
  final.table.list<-list()
  counts.list<-list()
  
  for(j in 1:length(top)){
    
    blast.filt.out<-paste(blastfile,"maxe",max_evalue,"qcov",min_qcovs,"top",top[j],sep = "_")
    
    filter.blast(blastfile = blastfile,ncbiTaxDir = ncbiTaxDir,out = blast.filt.out
                 ,headers = "qseqid sseqid evalue staxid pident qcovs",max_evalue = max_evalue,min_qcovs = min_qcovs,
                 top = top[j],rm.unclassified = F)
    
    btab.unfilt<-data.table::fread(blastfile,sep="\t",data.table = F)
    btab<-data.table::fread(blast.filt.out,sep="\t",data.table = F)
    
    final.table<-ecopcr.clean.df[,c("AC","path","taxids","origpath")]
    final.table$had.hits<-final.table$AC %in% btab.unfilt$V1
    final.table$after.filt<-final.table$AC %in% btab$qseqid
    
    counts<-data.frame(correct=rep(0,length(pidents)),
                       above=rep(0,length(pidents)),
                       incorrect=rep(0,length(pidents)),
                       fail.filt=rep(0,length(pidents)),
                       no.hits=rep(0,length(pidents)),
                       fail.bin=rep(0,length(pidents)),
                       pident=pidents,
                       file=blast.filt.out)
    
    for(i in 1:length(pidents)){
      
      pident<-pidents[i]
      message("binning with pident=",pident) 
      
      if(TaxlevelTest=="F"){
        btabf<-btab[!(btab$F=="unknown" & btab$G=="unknown" & btab$S=="unknown"),]   
        btabf<-btabf[!(btabf$F=="unknown" & btabf$G=="unknown"),] 
      }
      
      if(TaxlevelTest=="G") btabf<-btab[btab$G!="unknown",] 
      
      if(TaxlevelTest=="S") btabf<-btab
      
      btabf<-btabf[btabf$pident>pident,]
      
      if(nrow(btabf)>0){
        
        btabf$path<-paste(btabf$K,btabf$P,btabf$C,btabf$O,btabf$F,btabf$G,btabf$S,sep = ";")
        lcaf = aggregate(btabf$path, by=list(btabf$qseqid),function(x) lca(x,sep=";"))
        colnames(lcaf)<-c("qseqid","binpath")
        lcaf<-add.unknown.lca(lcaf)
        mat<-do.call(rbind,stringr::str_split(lcaf$binpath,";"))
        lcaf<-as.data.frame(cbind(lcaf$qseqid,mat[,1],mat[,2],mat[,3],mat[,4],mat[,5],mat[,6],mat[,7]))
        colnames(lcaf)<-c("qseqid","K","P","C","O","F","G","S")
        
      } else {
        lcaf<-data.frame(matrix(nrow=1,ncol = 8))
        colnames(lcaf)<-c("qseqid","K","P","C","O","F","G","S")
      }
      
      if(TaxlevelTest=="F") lcaf$binpath<-paste(lcaf$K,lcaf$P,lcaf$C,lcaf$O,lcaf$F,sep = ";")
      if(TaxlevelTest=="G") lcaf$binpath<-paste(lcaf$K,lcaf$P,lcaf$C,lcaf$O,lcaf$F,lcaf$G,sep = ";")
      if(TaxlevelTest=="S") lcaf$binpath<-paste(lcaf$K,lcaf$P,lcaf$C,lcaf$O,lcaf$F,lcaf$G,lcaf$S,sep = ";")
      
      colnames(lcaf)<-gsub("binpath",paste0("binpath_",pident),colnames(lcaf))
      
      final.table<-merge(final.table, lcaf[,c("qseqid",paste0("binpath_",pident))], by.x="AC",by.y = "qseqid",all.x = T)
      
      #correct 
      final.table$holder<-"holder"
      colnames(final.table)<-gsub("holder",paste0("correct_",pident),colnames(final.table))
      final.table[,paste0("correct_",pident)]<-final.table$origpath==final.table[,paste0("binpath_",pident)]
      counts$correct[i]<-sum(final.table[,paste0("correct_",pident)],na.rm = T)
      
      #above desired rank
      final.table$holder<-"holder"
      colnames(final.table)<-gsub("holder",paste0("rank_",pident),colnames(final.table))
      if(TaxlevelTest=="F"){
        final.table[,paste0("rank_",pident)]<-bas.get.ranks(data.frame(taxon=paste0(final.table[,paste0("binpath_",pident)],";NA;NA")))
      }
      if(TaxlevelTest=="G"){
        final.table[,paste0("rank_",pident)]<-bas.get.ranks(data.frame(taxon=paste0(final.table[,paste0("binpath_",pident)],";NA")))
      }
      if(TaxlevelTest=="S"){
        final.table[,paste0("rank_",pident)]<-bas.get.ranks(data.frame(taxon=final.table[,paste0("binpath_",pident)]))
      }
      
      #distinguish fail because of no_hits, fail_filt and fail_bin
      final.table[final.table$had.hits==F,paste0("rank_",pident)]<-"no.hits"
      final.table[final.table$had.hits==T & final.table$after.filt==F,paste0("rank_",pident)]<-"fail.filt"
      final.table[final.table$had.hits==T & final.table$after.filt==T & is.na(final.table[,paste0("binpath_",pident)]),
                  paste0("rank_",pident)]<-"fail.bin"
      
      if(TaxlevelTest=="F"){
        counts$above[i]<-nrow(final.table[final.table[,paste0("rank_",pident),]=="htf",])
      }
      if(TaxlevelTest=="G"){
        counts$above[i]<-nrow(final.table[final.table[,paste0("rank_",pident),]=="htf" | final.table[,paste0("rank_",pident),]=="family",])
      }
      if(TaxlevelTest=="S"){
        counts$above[i]<-nrow(final.table[final.table[,paste0("rank_",pident),]=="htf" | final.table[,paste0("rank_",pident),]=="family"
                                          | final.table[,paste0("rank_",pident),]=="genus",])
      }
      
      counts$fail.filt[i]<-nrow(final.table[final.table[,paste0("rank_",pident),]=="fail.filt",])
      counts$fail.bin[i]<-nrow(final.table[final.table[,paste0("rank_",pident),]=="fail.bin",])
      counts$no.hits[i]<-nrow(final.table[final.table[,paste0("rank_",pident),]=="no.hits",])
      
      #incorrect 
      final.table$holder<-"holder"
      colnames(final.table)<-gsub("holder",paste0("incorrect_",pident),colnames(final.table))
      
      if(TaxlevelTest=="F"){
        final.table[,paste0("incorrect_",pident)]<-
          final.table$origpath!=final.table[,paste0("binpath_",pident)] & final.table[,paste0("rank_",pident)]=="family"
      }
      if(TaxlevelTest=="G"){
        final.table[,paste0("incorrect_",pident)]<-
          final.table$origpath!=final.table[,paste0("binpath_",pident)] & final.table[,paste0("rank_",pident)]=="genus"
      }
      if(TaxlevelTest=="S"){
        final.table[,paste0("incorrect_",pident)]<-
          final.table$origpath!=final.table[,paste0("binpath_",pident)] & final.table[,paste0("rank_",pident)]=="species"
      }
      
      counts$incorrect[i]<-sum(final.table[,paste0("incorrect_",pident)],na.rm = T)
      
      #put all results in rank
      final.table[final.table[paste0("correct_",pident)]==T & !is.na(final.table[paste0("correct_",pident)]),
                  paste0("rank_",pident)]<-"correct"
      final.table[final.table[paste0("incorrect_",pident)]==T & !is.na(final.table[paste0("incorrect_",pident)]),
                  paste0("rank_",pident)]<-"incorrect"
      
      if(TaxlevelTest=="F") final.table[final.table[,paste0("rank_",pident)]=="htf",paste0("rank_",pident)]<-"above"
      if(TaxlevelTest=="G") final.table[final.table[,paste0("rank_",pident),]=="htf" | 
                                        final.table[,paste0("rank_",pident),]=="family",paste0("rank_",pident)]<-"above"
      if(TaxlevelTest=="S") final.table[final.table[,paste0("rank_",pident),]=="htf" | 
                                        final.table[,paste0("rank_",pident),]=="family" |
                                        final.table[,paste0("rank_",pident),]=="genus",paste0("rank_",pident)]<-"above"
      
      final.table$file<-blast.filt.out
      
      print(counts[i,])
      
      if(sum(counts[i,1:6])!=nrow(ecopcr.clean.df)) stop("Something wrong with counts, check function")
      
      }
    
    counts.list[[j]]<-counts
    final.table.list[[j]]<-final.table
  }
  
  counts.complete<-do.call(rbind,counts.list)
  final.table.complete<-do.call(rbind,final.table.list)
  output.list<-list(counts.complete,final.table.complete)
  
  t2<-Sys.time()
  
  message("Loop binning complete in ", round(difftime(t2,t1,units = "mins"),digits = 2)," mins")
  
  message("NOTE TO CLEAN FILES")
  
  return(output.list)
}

#' Bin reads into taxa based on blast results
#'
#' @param blastfile blast tabular output, e.g. -outfmt 6. Can have any columns, but must have:
#'     qcovs, evalue, staxid, qseqid, pident
#' @param headers A space-separated string of the header names of \code{blastfile}.
#'     Usually this is the same as the blast command used.
#' @param ncbiTaxDir The full path to the directory containing ncbi taxonomy. 
#' @param obitaxdb The full path to the obitools-formatted ncbi taxonomy. 
#' @param out specify file to output results
#' @param min_qcovs The minimum query coverage for all hits
#' @param max_evalue The maximum evalue for all hits
#' @param top The percentage pident to consider hits. i.e. For each query find the top hit pident (% identity),
#'     then find all pidents within "top" of this value, discard others - e.g. if top set at 1, and top hit for
#'     query A has pident=99.8, then any hits to query A with pident above 98.8 are kept
#' @param spident The minimum pident for binning hits at species level
#' @param gpident The minimum pident for binning hits at genus level
#' @param fpident The minimum pident for binning hits at family level
#' @param abspident An absolute pident for binning hits above family level. Hits below this threshold will be
#'     excluded completely
#'
#' @return A tab-separated file consisting of ESV name and associated taxonomic path
#' @note It works like this:
#' \itemize{
#'     \item 1. Apply min_qcovs (query coverage) to whole table
#'     \item 2. Apply max_evalue (blast evalue) to whole table
#'     \item 3. For each query find the top hit pident (% identity),
#'     then find all pidents within "top" of this value, discard others - e.g. if top set at 1,
#'     and top hit for query A has pident=99.8, then any hits with pident above 98.8 are kept.
#'     \item 4. Add taxon path to all remaining hits
#'     \item 5. For species-level binning: first remove hits that do not have species-level information
#'     (i.e. where database entries were only labelled at genus or higher), then remove hits with pident<spident.
#'     Then, for each query, find the lowest-common-ancestor (lca) of the remaining hits.
#'     \item 6. Repeat above for genus and family level binning
#'     \item 7. Do a final binning for hits that had taxa identified at higher-than-family level
#'     \item 6. Final report=if a query had species-level lca, use that. If not, use genus. If not genus, use family.
#'     If not family use absolute.
#'
#' @examples
#' blastfile<-"blast.results.txt"
#' headers<-"qseqid sacc sseqid sallseqid sallacc bitscore score staxid staxids sscinames scomnames salltitles length pident qcovs"
#' ncbiTaxDir<-"/Documents/TAXONOMIES/"
#' obitaxdb<-"/Documents/obitax_26-4-19"
#' bin.blast.bas(blastfile,headers,ncbiTaxDir,obitaxdb,out="blast_bins.txt",min_qcovs=70,max_evalue=0.001,top=1,spident=99,gpident=97,fpident=93)
#' @export
bin.blast2<-function(filtered_blastfile,ncbiTaxDir,
                     out,spident=98,gpident=95,fpident=92,abspident=80,
                     disabledTaxaFiles=NULL,disabledTaxaOut=NULL,
                     force=F,full.force=F,consider_sp.=F){
  t1<-Sys.time()
  
  if(is.null(out)) stop("out not specified")
  if(is.null(ncbiTaxDir)) stop("ncbiTaxDir not specified")

  require(treemap)
  
  ###################################################
  
  btab<-data.table::fread(filtered_blastfile,sep="\t",data.table = F)
  
  #preparing some things for final step
  total_hits<-length(btab$qseqid) #for info later
  total_queries<-length(unique(btab$qseqid))
  qseqids<-as.data.frame(unique(btab$qseqid))
  qseqids$qseqid<-qseqids$`unique(btab$qseqid)`
  qseqids$`unique(btab$qseqid)`=NULL
  
  #read and check disabled taxa file(s) 
  if(!is.null(disabledTaxaFiles)){
    
    disabledTaxaDf<-merge.and.check.disabled.taxa.files(disabledTaxaFiles,disabledTaxaOut,force = force,full.force = full.force)
    
    #get taxids at 7 levels
    disabledTaxaDf<-add.lineage.df(disabledTaxaDf,ncbiTaxDir,as.taxids = T)
    
    disabledSpecies<-disabledTaxaDf[disabledTaxaDf$disable_species==T,]
    disabledSpecies<-disabledSpecies[!is.na(disabledSpecies$taxids),]
    #get children of all disabled species
    disabledSpecies$taxids<-disabledSpecies$S
    if(nrow(disabledSpecies)>0)  {
      childrenS<-get.children.taxonkit(disabledSpecies) 
    } else {
      childrenS<-NULL
    } 
    
    disabledGenus<-disabledTaxaDf[disabledTaxaDf$disable_genus==T,]
    disabledGenus<-disabledGenus[!is.na(disabledGenus$taxids),]
    #get children of all disabled genera
    disabledGenus$taxids<-disabledGenus$G
    if(nrow(disabledGenus)>0)  {
      childrenG<-get.children.taxonkit(disabledGenus) 
    } else {
      childrenG<-NULL
    }
    
    
    disabledFamily<-disabledTaxaDf[disabledTaxaDf$disable_family==T,]
    disabledFamily<-disabledFamily[!is.na(disabledFamily$taxids),]
    #get children of all disabled families
    disabledFamily$taxids<-disabledFamily$F
    if(nrow(disabledFamily)>0)  {
      childrenF<-get.children.taxonkit(disabledFamily) 
    } else {
      childrenF<-NULL
    }
    
    message("The following taxa are disabled at species level")
    if(!is.null(childrenS)){
      childrenAlldf<-as.data.frame(unique(c(childrenS,childrenG,childrenF)))
      colnames(childrenAlldf)<-"taxids"
      childrenNames<-add.lineage.df(childrenAlldf,ncbiTaxDir)
      childrenNames$pathString<-paste("Family",childrenNames$F,childrenNames$G,childrenNames$S,sep = "/")
      childrenNames$pathString<-lapply(childrenNames$pathString, gsub, pattern = "unknown", replacement = "", fixed = TRUE)
      disabledtree <- data.tree::as.Node(childrenNames)
      print(disabledtree,limit = NULL)
    } else (message("No species disabled"))
    
    message("The following taxa are disabled at genus level")
    if(!is.null(childrenG)){
      childrenAlldf<-as.data.frame(unique(c(childrenG,childrenF)))
      colnames(childrenAlldf)<-"taxids"
      childrenNames<-add.lineage.df(childrenAlldf,ncbiTaxDir)
      childrenNames$pathString<-paste("Family",childrenNames$F,childrenNames$G,sep = "/")
      childrenNames$pathString<-lapply(childrenNames$pathString, gsub, pattern = "unknown", replacement = "", fixed = TRUE)
      disabledtree <- data.tree::as.Node(childrenNames)
      print(disabledtree,limit = NULL)
    } else (message("No genera disabled"))
    
    
    message("The following taxa are disabled at family level")
    if(!is.null(childrenF)){
      childrenAlldf<-as.data.frame(unique(c(childrenF)))
      colnames(childrenAlldf)<-"taxids"
      childrenNames<-add.lineage.df(childrenAlldf,ncbiTaxDir)
      childrenNames$pathString<-paste("Family",childrenNames$F,sep = "/")
      childrenNames$pathString<-lapply(childrenNames$pathString, gsub, pattern = "unknown", replacement = "", fixed = TRUE)
      disabledtree <- data.tree::as.Node(childrenNames)
      print(disabledtree,limit = NULL)
    } else (message("No families disabled"))
  }
  
  #species-level assignments
  message("binning at species level")
  
  btabsp<-btab[btab$S!="unknown",]
  
  if(!is.null(disabledTaxaFiles)){
    btabsp<-btabsp[!btabsp$taxids %in% unique(c(childrenS,childrenG,childrenF)),]
  }
  
  if(consider_sp.==F){
    message("Not considering species with 'sp.', numbers or more than one space")
    if(length(grep(" sp\\.",btabsp$S,ignore.case = T))>0) btabsp<-btabsp[-grep(" sp\\.",btabsp$S,ignore.case = T),]
    if(length(grep(" .* .*",btabsp$S,ignore.case = T))>0) btabsp<-btabsp[-grep(" .* .*",btabsp$S,ignore.case = T),]
    if(length(grep("[0-9]",btabsp$S))>0) btabsp<-btabsp[-grep("[0-9]",btabsp$S),]
  } else(message("Considering species with 'sp.', numbers or more than one space"))
  
  btabsp<-btabsp[btabsp$pident>spident,]
  
  if(nrow(btabsp)>0){
    
    btabsp$path<-paste(btabsp$K,btabsp$P,btabsp$C,btabsp$O,btabsp$F,btabsp$G,btabsp$S,sep = ";")
    lcasp = aggregate(btabsp$path, by=list(btabsp$qseqid),function(x) lca(x,sep=";"))
    colnames(lcasp)<-c("qseqid","binpath")
    lcasp<-add.unknown.lca(lcasp)
    mat<-do.call(rbind,stringr::str_split(lcasp$binpath,";"))
    lcasp<-as.data.frame(cbind(lcasp$qseqid,mat[,1],mat[,2],mat[,3],mat[,4],mat[,5],mat[,6],mat[,7]))
    colnames(lcasp)<-c("qseqid","K","P","C","O","F","G","S")
    
  } else {
    lcasp<-data.frame(matrix(nrow=1,ncol = 8))
    colnames(lcasp)<-c("qseqid","K","P","C","O","F","G","S")
  }
  
  rm(btabsp)
  
  #genus-level assignments
  message("binning at genus level") 
  
  btabg<-btab[btab$G!="unknown",]  
  #reason - can have g=unknown and s=known (e.g. Ranidae isolate), these should be removed
  #can have g=unknown and s=unknown (e.g. Ranidae), these should be removed
  #can have g=known and s=unknown (e.g. Leiopelma), these should be kept
  
  if(!is.null(disabledTaxaFiles)){
    btabg<-btabg[!btabg$taxids %in% unique(c(childrenG,childrenF)),]
  }
  
  btabg<-btabg[btabg$pident>gpident,]  
  
  if(nrow(btabg)>0){
    
    btabg$path<-paste(btabg$K,btabg$P,btabg$C,btabg$O,btabg$F,btabg$G,btabg$S,sep = ";")
    lcag = aggregate(btabg$path, by=list(btabg$qseqid),function(x) lca(x,sep=";"))
    colnames(lcag)<-c("qseqid","binpath")
    lcag<-add.unknown.lca(lcag)
    mat<-do.call(rbind,stringr::str_split(lcag$binpath,";"))
    lcag<-as.data.frame(cbind(lcag$qseqid,mat[,1],mat[,2],mat[,3],mat[,4],mat[,5],mat[,6],mat[,7]))
    colnames(lcag)<-c("qseqid","K","P","C","O","F","G","S")
    
  } else {
    lcag<-data.frame(matrix(nrow=1,ncol = 8))
    colnames(lcag)<-c("qseqid","K","P","C","O","F","G","S")
  }
  
  rm(btabg)
  
  #family-level assignments
  message("binning at family level") 
  #can have f=known, g=unknown, s=unknown, these should be kept
  #can have f=unknown, g=known, s=known, these should be kept
  #can have f=unknown, g=known, s=unknown, these should be kept
  #can have f=known, g=known, s=unknown, these should be kept
  #can have f=known, g=known, s=known, these should be kept
  #can have f=unknown, g=known, s=unknown, these should be kept
  #can have f=known, g=unknown, s=known, these should be kept
  
  #can have f=unknown, g=unknown, s=known, these should be removed - 
  #assumes that this case would be a weird entry (e.g. Ranidae isolate)
  
  #can have f=unknown, g=unknown, s=unknown, these should be removed
  
  btabf<-btab[!(btab$F=="unknown" & btab$G=="unknown" & btab$S=="unknown"),]  ####line changed 
  btabf<-btabf[!(btabf$F=="unknown" & btabf$G=="unknown"),] ####line changed 
  
  #this is ok, but when ncbi taxonomy does not have family-level assignment, we only get to order, stupid, or rather we do get to 
  #family or subfamily as lca, but final report puts it to "unknown"
  
  if(!is.null(disabledTaxaFiles)){
    btabf<-btabf[!btabf$taxids %in% unique(c(childrenF)),]
  }
  
  btabf<-btabf[btabf$pident>fpident,]
  
  if(nrow(btabf)>0){
    
    btabf$path<-paste(btabf$K,btabf$P,btabf$C,btabf$O,btabf$F,btabf$G,btabf$S,sep = ";")
    lcaf = aggregate(btabf$path, by=list(btabf$qseqid),function(x) lca(x,sep=";"))
    colnames(lcaf)<-c("qseqid","binpath")
    lcaf<-add.unknown.lca(lcaf)
    mat<-do.call(rbind,stringr::str_split(lcaf$binpath,";"))
    lcaf<-as.data.frame(cbind(lcaf$qseqid,mat[,1],mat[,2],mat[,3],mat[,4],mat[,5],mat[,6],mat[,7]))
    colnames(lcaf)<-c("qseqid","K","P","C","O","F","G","S")
    
  } else {
    lcaf<-data.frame(matrix(nrow=1,ncol = 8))
    colnames(lcaf)<-c("qseqid","K","P","C","O","F","G","S")
  }
  
  rm(btabf)
  
  #higher-than-family-level assignments
  message("binning at higher-than-family level")
  btabhtf<-btab[btab$K!="unknown",]
  
  # if(!is.null(disabledTaxaFiles)){
  #   btabhtf<-btabhtf[!btabhtf$taxids %in% childrenAll,]
  # }
  
  btabhtf<-btabhtf[btabhtf$pident>abspident,]
  
  if(nrow(btabhtf)>0){
    
    btabhtf$path<-paste(btabhtf$K,btabhtf$P,btabhtf$C,btabhtf$O,btabhtf$F,btabhtf$G,btabhtf$S,sep = ";")
    lcahtf = aggregate(btabhtf$path, by=list(btabhtf$qseqid),function(x) lca(x,sep=";"))
    colnames(lcahtf)<-c("qseqid","binpath")
    lcahtf<-add.unknown.lca(lcahtf)
    mat<-do.call(rbind,stringr::str_split(lcahtf$binpath,";"))
    lcahtf<-as.data.frame(cbind(lcahtf$qseqid,mat[,1],mat[,2],mat[,3],mat[,4],mat[,5],mat[,6],mat[,7]))
    colnames(lcahtf)<-c("qseqid","K","P","C","O","F","G","S")
    
  } else {
    lcahtf<-data.frame(matrix(nrow=1,ncol = 8))
    colnames(lcahtf)<-c("qseqid","K","P","C","O","F","G","S")
  }
  
  rm(btabhtf)
  
  ###################################################
  #combine
  sp_level<-lcasp[lcasp$S!="unknown",]
  g_level<-lcag[lcag$G!="unknown",]
  if(nrow(g_level)>0) g_level$S<-NA
  g_level<-g_level[!g_level$qseqid %in% sp_level$qseqid,]
  f_level<-lcaf[lcaf$F!="unknown",]
  if(nrow(f_level)>0) f_level$G<-NA
  if(nrow(f_level)>0) f_level$S<-NA
  f_level<-f_level[!f_level$qseqid %in% sp_level$qseqid,]
  f_level<-f_level[!f_level$qseqid %in% g_level$qseqid,]
  
  abs_level<-lcahtf
  if(nrow(abs_level)>0) abs_level$G<-NA
  if(nrow(abs_level)>0) abs_level$S<-NA
  if(nrow(abs_level)>0) abs_level$F<-NA
  abs_level<-abs_level[!abs_level$qseqid %in% sp_level$qseqid,]
  abs_level<-abs_level[!abs_level$qseqid %in% g_level$qseqid,]
  abs_level<-abs_level[!abs_level$qseqid %in% f_level$qseqid,]
  
  com_level<-rbind(sp_level,g_level,f_level,abs_level)
  com_level<-merge(x=qseqids, y = com_level, by = "qseqid",all.x = T)
  
  #info
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  
  write.table(x = com_level,file = out,sep="\t",quote = F,row.names = F)
  
  message(c("Complete. ",total_hits, " hits from ", total_queries," queries processed in ",t3," mins."))
  
  message("Note: if all hits for a particular OTU are removed due to filters, 
        the results will be NA for all taxon levels.
        If the lca for a particular OTU is above kingdom, e.g. cellular organisms or root, 
        the results will be unknown for all taxon levels.")
}

#############################decided to split function
filter.blast<-function(blastfile,headers="qseqid evalue staxid pident qcovs",ncbiTaxDir,out,min_qcovs=70,
                       max_evalue=0.001,top=1,rm.unclassified=T){
  
  if(length(grep("qcovs",headers))<1) stop("qcovs not in headers")
  if(length(grep("evalue",headers))<1) stop("evalue not in headers")
  if(length(grep("qseqid",headers))<1) stop("qseqid not in headers")
  if(length(grep("pident",headers))<1) stop("pident not in headers")
  if(length(grep("staxid",headers))<1) stop("staxid not in headers")
  
  if(is.null(ncbiTaxDir)) stop("ncbiTaxDir not specified")
  if(is.null(out)) stop("out not specified")
  
  message("reading blast results")
  btab<-as.data.frame(data.table::fread(file = blastfile,sep = "\t"))
  colnames(btab)<-strsplit(headers,split = " ")[[1]]
  
  message("applying global min_qcovs threshold")
  start<-length(unique(btab$qseqid))
  btab<-btab[btab$qcovs>min_qcovs,]
  after_min_qcovs<-length(unique(btab$qseqid))
  message(start-after_min_qcovs, " queries removed")
  
  message("applying global max_evalue threshold")
  btab<-btab[btab$evalue<max_evalue,]
  after_max_eval<-length(unique(btab$qseqid))
  message(after_min_qcovs-after_max_eval, " queries removed")
  
  message("applying global top threshold")
  if(length(btab[,1])==0) stop("No hits passing min_qcovs and max_evalue thresholds")
  topdf<-aggregate(x = btab[,"pident"],by=list(btab$qseqid),FUN = max)
  colnames(topdf)<-c("qseqid","pident")
  topdf$min_pident<-topdf$pident-top
  btab$min_pident<-topdf[match(btab$qseqid,topdf$qseqid),"min_pident"]
  btab<-btab[btab$pident>btab$min_pident,]
  
  #add lineage to results
  message("adding taxonomic lineages")
  btab$taxids<-btab$staxid #add.lineage.df requires this colname
  btab<-add.lineage.df(btab,ncbiTaxDir)
  
  btab$K<-gsub(" ","_",btab$K)
  btab$P<-gsub(" ","_",btab$P)
  btab$C<-gsub(" ","_",btab$C)
  btab$O<-gsub(" ","_",btab$O)
  btab$F<-gsub(" ","_",btab$F)
  btab$G<-gsub(" ","_",btab$G)
  
  #remove crappy hits 
  if(rm.unclassified==T){
    #1. btab$S contains uncultured
    message("Removing species containing the terms: uncultured, environmental, 
            unidentified,fungal, eukaryote or unclassified")
    if(length(grep("uncultured",btab$S,ignore.case = T))>0) btab<-btab[-grep("uncultured",btab$S,ignore.case = T),]
    if(length(grep("environmental",btab$S,ignore.case = T))>0) btab<-btab[-grep("environmental",btab$S,ignore.case = T),]
    if(length(grep("unclassified",btab$S,ignore.case = T))>0) btab<-btab[-grep("unclassified",btab$S,ignore.case = T),]
    if(length(grep("unidentified",btab$S,ignore.case = T))>0) btab<-btab[-grep("unidentified",btab$S,ignore.case = T),]
    if(length(grep("fungal ",btab$S,ignore.case = T))>0) btab<-btab[-grep("fungal ",btab$S,ignore.case = T),]
    if(length(grep("eukaryote",btab$S,ignore.case = T))>0) btab<-btab[-grep("eukaryote",btab$S,ignore.case = T),]
  }
  
  write.table(x = btab,file = out,sep="\t",quote = F,row.names = F)
  
}


merge.and.check.disabled.taxa.files<-function(disabledTaxaFiles,disabledTaxaOut,force=F,full.force=F){
  
  message("Note: Use force=T to ignore any contributor entries where no levels were disabled when consistency checking.
                 Use force=F to do more thorough consistency checks")
  
  disabledTaxaDFList<-list()
  
  for(i in 1:length(disabledTaxaFiles)){
    
    disabledTaxaDFList[[i]]<-data.table::fread(disabledTaxaFiles[i], data.table = F,sep = "\t")
    
    if(!"taxids" %in% colnames(disabledTaxaDFList[[i]]))  stop(paste("No column called 'taxids' in", disabledTaxaFiles[i]))
    if(!"disable_species" %in% colnames(disabledTaxaDFList[[i]])) stop(paste("No column called 'disable_species' in", disabledTaxaFiles[i]))
    if(!"disable_genus" %in% colnames(disabledTaxaDFList[[i]])) stop(paste("No column called 'disable_genus' in", disabledTaxaFiles[i]))
    if(!"disable_family" %in% colnames(disabledTaxaDFList[[i]])) stop(paste("No column called 'disable_family' in", disabledTaxaFiles[i]))
    if(!"contributors" %in% colnames(disabledTaxaDFList[[i]])) stop(paste("No column called 'contributors' in", disabledTaxaFiles[i]))
    
    disabledTaxaDFList[[i]]<-disabledTaxaDFList[[i]][,c("contributors","taxids","disable_species","disable_genus","disable_family")]
    disabledTaxaDFList[[i]]$file<-disabledTaxaFiles[i]
  }
  
  disabledTaxaDf<-do.call(rbind,disabledTaxaDFList)
  
  #convert to logical
  disabled.cols<-c("disable_species","disable_genus","disable_family")
  for(i in 1:3){
    disabledTaxaDf[,disabled.cols[i]]<-as.logical(disabledTaxaDf[,disabled.cols[i]])
    disabledTaxaDf[,disabled.cols[i]][is.na(disabledTaxaDf[,disabled.cols[i]])]<-"FALSE"
  }
  
  disabledTaxaDf<-add.lineage.df(disabledTaxaDf,ncbiTaxDir,as.taxids = T)
  
  if(full.force) force=T
  
  if(force) {
    message("Using force=T")
    disabledTaxaDf$trues<-rowSums(disabledTaxaDf[,c("disable_species","disable_genus","disable_family")])
    disabledTaxaDf<-disabledTaxaDf[disabledTaxaDf$trues>0,-match("trues",colnames(disabledTaxaDf))]
    
    if(full.force) {
      message("Using full.force=T")
      
      #change families marked as TRUE to TRUE 
      splitdf<-split(disabledTaxaDf, f = as.factor(disabledTaxaDf$F))
      
      for (i in 1:(length(splitdf))) {
        if(TRUE %in% splitdf[[i]]$disable_family) splitdf[[i]]$disable_family<-TRUE
      }
      
      disabledTaxaDf<-do.call(rbind, splitdf)
      
      #change genus marked as TRUE to TRUE 
      splitdf<-split(disabledTaxaDf, f = as.factor(disabledTaxaDf$G))
      
      for (i in 1:(length(splitdf))) {
        if(TRUE %in% splitdf[[i]]$disable_genus) splitdf[[i]]$disable_genus<-TRUE
      }
      
      disabledTaxaDf<-do.call(rbind, splitdf)
      
      #change species marked as TRUE to TRUE 
      splitdf<-split(disabledTaxaDf, f = as.factor(disabledTaxaDf$S))
      
      for (i in 1:(length(splitdf))) {
        if(TRUE %in% splitdf[[i]]$disable_species) splitdf[[i]]$disable_species<-TRUE
      }
      
      disabledTaxaDf<-do.call(rbind, splitdf)
    }
    
  } else { message("Using force=F")}
  
  
  #check that identical paths have not been treated differently
  shouldstop1<-list()
  for(i in 1:length(unique(disabledTaxaDf$contributors))){
    temp<-disabledTaxaDf[disabledTaxaDf$contributors==unique(disabledTaxaDf$contributors)[i],]
    if(length(temp$contributors)>1) if(sum(duplicated(temp[,c("disable_species","disable_genus","disable_family")]))!=
                                       length(temp$contributors)-1) {
      message("inconsistent taxid disabling detected")
      print(temp)
      shouldstop1[[i]]<-temp
    }
  }
  
  #check that identical families have not been treated differently
  shouldstop2<-list()
  for(i in 1:length(unique(disabledTaxaDf$F))){
    temp<-disabledTaxaDf[disabledTaxaDf$F==unique(disabledTaxaDf$F)[i],]
    if(length(temp$contributors)>1) if(sum(duplicated(temp[,"disable_family"]))!=length(temp$contributors)-1) {
      message("inconsistent family disabling detected:")
      print(temp)
      shouldstop2[[i]]<-temp
    }
  }
  
  #check that identical genera have not been treated differently
  shouldstop3<-list()
  for(i in 1:length(unique(disabledTaxaDf$G))){
    temp<-disabledTaxaDf[disabledTaxaDf$G==unique(disabledTaxaDf$G)[i],]
    if(length(temp$contributors)>1) if(sum(duplicated(temp[,"disable_genus"]))!=length(temp$contributors)-1) {
      message("inconsistent genus disabling detected")
      print(temp)
      shouldstop3[[i]]<-temp
    }
  }
  
  #check that identical species have not been treated differently
  shouldstop4<-list()
  for(i in 1:length(unique(disabledTaxaDf$S))){
    temp<-disabledTaxaDf[disabledTaxaDf$S==unique(disabledTaxaDf$S)[i],]
    if(length(temp$contributors)>1) if(sum(duplicated(temp[,"disable_species"]))!=length(temp$contributors)-1) {
      message("inconsistent species disabling detected")
      print(temp)
      shouldstop4[[i]]<-temp
    }
  }
  
  if(length(shouldstop1)>0) for(i in 1:length(shouldstop1)) if(!is.null(shouldstop1[[i]])) stop("fix inconsistencies")
  if(length(shouldstop2)>0) for(i in 1:length(shouldstop2)) if(!is.null(shouldstop2[[i]])) stop("fix inconsistencies")
  if(length(shouldstop3)>0) for(i in 1:length(shouldstop3)) if(!is.null(shouldstop3[[i]])) stop("fix inconsistencies")
  if(length(shouldstop4)>0) for(i in 1:length(shouldstop4)) if(!is.null(shouldstop4[[i]])) stop("fix inconsistencies")
  
  #check that if genus disabled, species within that genus also should be 
  temp<-disabledTaxaDf[(disabledTaxaDf$disable_genus==T & disabledTaxaDf$disable_species==F),"contributors"]
  if(length(temp)>0) {
    print(temp)
    stop("Genus was disabled for this contributor, but species was not")
  }
  
  #check that if family disabled, species within that family also should be 
  temp<-disabledTaxaDf[(disabledTaxaDf$disable_family==T & disabledTaxaDf$disable_species==F),"contributors"]
  if(length(temp)>0) {
    print(temp)
    stop("Family was disabled for this contributor, but species was not")
  }
  
  #check that if family disabled, genus within that family also should be 
  temp<-disabledTaxaDf[(disabledTaxaDf$disable_family==T & disabledTaxaDf$disable_genus==F),"contributors"]
  if(length(temp)>0) {
    print(temp)
    stop("Family was disabled for this contributor, but genus was not")
  }
  
  if(!is.null(disabledTaxaOut)) write.table(disabledTaxaDf,disabledTaxaOut,col.names = T,row.names = F,quote = F,sep = "\t")
  
  return(disabledTaxaDf)
}

bin.blast.lite<-function(filtered_blastfile,ncbiTaxDir,obitaxdb,out){
  
  t1<-Sys.time()
  
  if(is.null(out)) stop("out not specified")
  if(is.null(ncbiTaxDir)) stop("ncbiTaxDir not specified")
  if(is.null(obitaxdb)) stop("obitaxdb not specified")
  
  ###################################################
  #read in obitools taxonomy
  obitaxdb2=ROBITaxonomy::read.taxonomy(obitaxdb)
  
  btab<-data.table::fread(filtered_blastfile,sep="\t",data.table = F)
  
  #preparing some things for final step
  total_hits<-length(btab$qseqid) #for info later
  total_queries<-length(unique(btab$qseqid))
  qseqids<-as.data.frame(unique(btab$qseqid))
  qseqids$qseqid<-qseqids$`unique(btab$qseqid)`
  qseqids$`unique(btab$qseqid)`=NULL
  
  message("binning")
  btabsp<-btab[btab$S!="unknown",]
  
  message("Considering species with 'sp.', numbers or more than one space")
  
  if(length(btabsp$taxids)>0){
    lcasp = aggregate(btabsp$taxids, by=list(btabsp$qseqid),function(x) ROBITaxonomy::lowest.common.ancestor(obitaxdb2,x))
    
    #get lca names
    colnames(lcasp)<-gsub("x","taxids",colnames(lcasp))
    if(sum(is.na(lcasp$taxids))>0){
      message("************
              ERROR: Some taxids were not recognized by ROBITaxonomy::lowest.common.ancestor, probably need to update obitaxdb using NCBI2obitaxonomy
              *************")
      problem.contributors<-btabsp[!duplicated(btabsp[lcasp[is.na(lcasp$taxids),1] %in% btabsp$qseqid,c("taxids","K","P","C","O","F","G","S")]),
                                   c("taxids","K","P","C","O","F","G","S")]
      for(i in 1:length(problem.contributors$taxids)){
        problem.contributors$validate[i]<-ROBITaxonomy::validate(obitaxdb2,problem.contributors$taxids[i])
      }
      print(problem.contributors[is.na(problem.contributors$validate),])
    }
    lcasp<-add.lineage.df(df = lcasp,ncbiTaxDir)
    colnames(lcasp)<-gsub("Group.1","qseqid",colnames(lcasp))
  } else {
    lcasp<-data.frame(matrix(nrow=1,ncol = 10))
    colnames(lcasp)<-c("taxids","qseqid","old_taxids","K","P","C","O","F","G","S")
  }
  
  
  ###################################################
  
  com_level<-merge(x=qseqids, y = lcasp[,c(2,4:10)], by = "qseqid",all = T)
  
  com_level[com_level=="unknown"]<-NA
  
  #info
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  
  write.table(x = com_level,file = out,sep="\t",quote = F,row.names = F)
  
  message(c("Complete. ",total_hits, " hits from ", total_queries," queries processed in ",t3," mins."))
  
  message("Note: if all hits for a particular OTU are removed due to filters, 
          the results will be NA for all taxon levels.
          If the lca for a particular OTU is above kingdom, e.g. cellular organisms or root, 
          the results will be unknown for all taxon levels.")
}

filter.blast2<-function(blastfile,headers="qseqid evalue staxid pident qcovs",ncbiTaxDir,out,min_qcovs=70,
                        max_evalue=0.001,rm.unclassified=T){
  
  message("Reminder, filter.blast2 does not include 'top' threshold")
  
  if(length(grep("qcovs",headers))<1) stop("qcovs not in headers")
  if(length(grep("evalue",headers))<1) stop("evalue not in headers")
  if(length(grep("qseqid",headers))<1) stop("qseqid not in headers")
  if(length(grep("pident",headers))<1) stop("pident not in headers")
  if(length(grep("staxid",headers))<1) stop("staxid not in headers")
  
  if(is.null(ncbiTaxDir)) stop("ncbiTaxDir not specified")
  if(is.null(out)) stop("out not specified")
  
  message("reading blast results")
  btab<-as.data.frame(data.table::fread(file = blastfile,sep = "\t"))
  colnames(btab)<-strsplit(headers,split = " ")[[1]]
  
  message("applying global min_qcovs threshold")
  start<-nrow(btab)
  btab<-btab[btab$qcovs>min_qcovs,]
  after_min_qcovs<-nrow(btab)
  message(start-after_min_qcovs, " hits removed")
  
  message("applying global max_evalue threshold")
  btab<-btab[btab$evalue<max_evalue,]
  after_max_eval<-nrow(btab)
  message(after_min_qcovs-after_max_eval, " hits removed")
  
  #add lineage to results
  message("adding taxonomic lineages")
  btab$taxids<-btab$staxid #add.lineage.df requires this colname
  btab<-add.lineage.df(btab,ncbiTaxDir)
  
  btab$K<-gsub(" ","_",btab$K)
  btab$P<-gsub(" ","_",btab$P)
  btab$C<-gsub(" ","_",btab$C)
  btab$O<-gsub(" ","_",btab$O)
  btab$F<-gsub(" ","_",btab$F)
  btab$G<-gsub(" ","_",btab$G)
  
  #remove crappy hits 
  if(rm.unclassified==T){
    #1. btab$S contains uncultured
    message("Removing species containing the terms: uncultured, environmental, 
            unidentified,fungal, eukaryote or unclassified")
    if(length(grep("uncultured",btab$S,ignore.case = T))>0) btab<-btab[-grep("uncultured",btab$S,ignore.case = T),]
    if(length(grep("environmental",btab$S,ignore.case = T))>0) btab<-btab[-grep("environmental",btab$S,ignore.case = T),]
    if(length(grep("unclassified",btab$S,ignore.case = T))>0) btab<-btab[-grep("unclassified",btab$S,ignore.case = T),]
    if(length(grep("unidentified",btab$S,ignore.case = T))>0) btab<-btab[-grep("unidentified",btab$S,ignore.case = T),]
    if(length(grep("fungal ",btab$S,ignore.case = T))>0) btab<-btab[-grep("fungal ",btab$S,ignore.case = T),]
    if(length(grep("eukaryote",btab$S,ignore.case = T))>0) btab<-btab[-grep("eukaryote",btab$S,ignore.case = T),]
    if(length(grep("synthetic",btab$S,ignore.case = T))>0) btab<-btab[-grep("synthetic",btab$S,ignore.case = T),]
  }
  
  write.table(x = btab,file = out,sep="\t",quote = F,row.names = F)
  
}

bin.blast3<-function(filtered_blastfile,ncbiTaxDir,
                     out,spident=98,gpident=95,fpident=92,abspident=80,topS=1,topG=1,topF=1,topAbs=1,
                     disabledTaxaFiles=NULL,disabledTaxaOut=NULL,
                     force=F,full.force=F,consider_sp.=F){
  t1<-Sys.time()
  
  if(is.null(out)) stop("out not specified")
  if(is.null(ncbiTaxDir)) stop("ncbiTaxDir not specified")

  ###################################################
  
  btab<-data.table::fread(filtered_blastfile,sep="\t",data.table = F)
  
  #preparing some things for final step
  total_hits<-length(btab$qseqid) #for info later
  total_queries<-length(unique(btab$qseqid))
  qseqids<-as.data.frame(unique(btab$qseqid))
  qseqids$qseqid<-qseqids$`unique(btab$qseqid)`
  qseqids$`unique(btab$qseqid)`=NULL
  
  #read and check disabled taxa file(s) 
  if(!is.null(disabledTaxaFiles)){
    
    require(treemap)
    
    disabledTaxaDf<-merge.and.check.disabled.taxa.files(disabledTaxaFiles,disabledTaxaOut,force = force,full.force = full.force)
    
    #get taxids at 7 levels
    disabledTaxaDf<-add.lineage.df(disabledTaxaDf,ncbiTaxDir,as.taxids = T)
    
    disabledSpecies<-disabledTaxaDf[disabledTaxaDf$disable_species==T,]
    disabledSpecies<-disabledSpecies[!is.na(disabledSpecies$taxids),]
    #get children of all disabled species
    disabledSpecies$taxids<-disabledSpecies$S
    if(nrow(disabledSpecies)>0)  {
      childrenS<-get.children.taxonkit(disabledSpecies) 
    } else {
      childrenS<-NULL
    } 
    
    disabledGenus<-disabledTaxaDf[disabledTaxaDf$disable_genus==T,]
    disabledGenus<-disabledGenus[!is.na(disabledGenus$taxids),]
    #get children of all disabled genera
    disabledGenus$taxids<-disabledGenus$G
    if(nrow(disabledGenus)>0)  {
      childrenG<-get.children.taxonkit(disabledGenus) 
    } else {
      childrenG<-NULL
    }
    
    
    disabledFamily<-disabledTaxaDf[disabledTaxaDf$disable_family==T,]
    disabledFamily<-disabledFamily[!is.na(disabledFamily$taxids),]
    #get children of all disabled families
    disabledFamily$taxids<-disabledFamily$F
    if(nrow(disabledFamily)>0)  {
      childrenF<-get.children.taxonkit(disabledFamily) 
    } else {
      childrenF<-NULL
    }
    
    message("The following taxa are disabled at species level")
    if(!is.null(childrenS)){
      childrenAlldf<-as.data.frame(unique(c(childrenS,childrenG,childrenF)))
      colnames(childrenAlldf)<-"taxids"
      childrenNames<-add.lineage.df(childrenAlldf,ncbiTaxDir)
      childrenNames$pathString<-paste("Family",childrenNames$F,childrenNames$G,childrenNames$S,sep = "/")
      childrenNames$pathString<-lapply(childrenNames$pathString, gsub, pattern = "unknown", replacement = "", fixed = TRUE)
      disabledtree <- data.tree::as.Node(childrenNames)
      print(disabledtree,limit = NULL)
    } else (message("No species disabled"))
    
    message("The following taxa are disabled at genus level")
    if(!is.null(childrenG)){
      childrenAlldf<-as.data.frame(unique(c(childrenG,childrenF)))
      colnames(childrenAlldf)<-"taxids"
      childrenNames<-add.lineage.df(childrenAlldf,ncbiTaxDir)
      childrenNames$pathString<-paste("Family",childrenNames$F,childrenNames$G,sep = "/")
      childrenNames$pathString<-lapply(childrenNames$pathString, gsub, pattern = "unknown", replacement = "", fixed = TRUE)
      disabledtree <- data.tree::as.Node(childrenNames)
      print(disabledtree,limit = NULL)
    } else (message("No genera disabled"))
    
    
    message("The following taxa are disabled at family level")
    if(!is.null(childrenF)){
      childrenAlldf<-as.data.frame(unique(c(childrenF)))
      colnames(childrenAlldf)<-"taxids"
      childrenNames<-add.lineage.df(childrenAlldf,ncbiTaxDir)
      childrenNames$pathString<-paste("Family",childrenNames$F,sep = "/")
      childrenNames$pathString<-lapply(childrenNames$pathString, gsub, pattern = "unknown", replacement = "", fixed = TRUE)
      disabledtree <- data.tree::as.Node(childrenNames)
      print(disabledtree,limit = NULL)
    } else (message("No families disabled"))
  }
  
  #species-level assignments
  message("binning at species level")
  
  btabsp<-btab[btab$S!="unknown",]
  
  message("applying species top threshold of ",topS)
  topdf<-aggregate(x = btabsp[,"pident"],by=list(btabsp$qseqid),FUN = max)
  colnames(topdf)<-c("qseqid","pident")
  topdf$min_pident<-topdf$pident-topS
  btabsp$min_pident<-topdf[match(btabsp$qseqid,topdf$qseqid),"min_pident"]
  btabsp<-btabsp[btabsp$pident>btabsp$min_pident,]
  
  if(!is.null(disabledTaxaFiles)){
    btabsp<-btabsp[!btabsp$taxids %in% unique(c(childrenS,childrenG,childrenF)),]
  }
  
  if(consider_sp.==F){
    message("Not considering species with 'sp.', numbers or more than one space")
    if(length(grep(" sp\\.",btabsp$S,ignore.case = T))>0) btabsp<-btabsp[-grep(" sp\\.",btabsp$S,ignore.case = T),]
    if(length(grep(" .* .*",btabsp$S,ignore.case = T))>0) btabsp<-btabsp[-grep(" .* .*",btabsp$S,ignore.case = T),]
    if(length(grep("[0-9]",btabsp$S))>0) btabsp<-btabsp[-grep("[0-9]",btabsp$S),]
  } else(message("Considering species with 'sp.', numbers or more than one space"))
  
  btabsp<-btabsp[btabsp$pident>spident,]
 
  if(nrow(btabsp)>0){
    
    btabsp$path<-paste(btabsp$K,btabsp$P,btabsp$C,btabsp$O,btabsp$F,btabsp$G,btabsp$S,sep = ";")
    lcasp = aggregate(btabsp$path, by=list(btabsp$qseqid),function(x) lca(x,sep=";"))
    colnames(lcasp)<-c("qseqid","binpath")
    lcasp<-add.unknown.lca(lcasp)
    mat<-do.call(rbind,stringr::str_split(lcasp$binpath,";"))
    lcasp<-as.data.frame(cbind(lcasp$qseqid,mat[,1],mat[,2],mat[,3],mat[,4],mat[,5],mat[,6],mat[,7]))
    colnames(lcasp)<-c("qseqid","K","P","C","O","F","G","S")
    
  } else {
    lcasp<-data.frame(matrix(nrow=1,ncol = 8))
    colnames(lcasp)<-c("qseqid","K","P","C","O","F","G","S")
  }
  rm(btabsp)
  
  #genus-level assignments
  message("binning at genus level") 
  
  btabg<-btab[btab$G!="unknown",]  
  #reason - can have g=unknown and s=known (e.g. Ranidae isolate), these should be removed
  #can have g=unknown and s=unknown (e.g. Ranidae), these should be removed
  #can have g=known and s=unknown (e.g. Leiopelma), these should be kept
  
  message("applying genus top threshold of ",topG)
  topdf<-aggregate(x = btabg[,"pident"],by=list(btabg$qseqid),FUN = max)
  colnames(topdf)<-c("qseqid","pident")
  topdf$min_pident<-topdf$pident-topG
  btabg$min_pident<-topdf[match(btabg$qseqid,topdf$qseqid),"min_pident"]
  btabg<-btabg[btabg$pident>btabg$min_pident,]
  
  if(!is.null(disabledTaxaFiles)){
    btabg<-btabg[!btabg$taxids %in% unique(c(childrenG,childrenF)),]
  }
  
  btabg<-btabg[btabg$pident>gpident,] 
  
  if(nrow(btabg)>0){
    
    btabg$path<-paste(btabg$K,btabg$P,btabg$C,btabg$O,btabg$F,btabg$G,btabg$S,sep = ";")
    lcag = aggregate(btabg$path, by=list(btabg$qseqid),function(x) lca(x,sep=";"))
    colnames(lcag)<-c("qseqid","binpath")
    lcag<-add.unknown.lca(lcag)
    mat<-do.call(rbind,stringr::str_split(lcag$binpath,";"))
    lcag<-as.data.frame(cbind(lcag$qseqid,mat[,1],mat[,2],mat[,3],mat[,4],mat[,5],mat[,6],mat[,7]))
    colnames(lcag)<-c("qseqid","K","P","C","O","F","G","S")
    
  } else {
    lcag<-data.frame(matrix(nrow=1,ncol = 8))
    colnames(lcag)<-c("qseqid","K","P","C","O","F","G","S")
  }
  
  rm(btabg)
  
  #family-level assignments
  message("binning at family level") 
  #can have f=known, g=unknown, s=unknown, these should be kept
  #can have f=unknown, g=known, s=known, these should be kept
  #can have f=unknown, g=known, s=unknown, these should be kept
  #can have f=known, g=known, s=unknown, these should be kept
  #can have f=known, g=known, s=known, these should be kept
  #can have f=unknown, g=known, s=unknown, these should be kept
  #can have f=known, g=unknown, s=known, these should be kept
  
  #can have f=unknown, g=unknown, s=known, these should be removed - 
  #assumes that this case would be a weird entry (e.g. Ranidae isolate)
  
  #can have f=unknown, g=unknown, s=unknown, these should be removed
  
  btabf<-btab[!(btab$F=="unknown" & btab$G=="unknown" & btab$S=="unknown"),]  ####line changed 
  btabf<-btabf[!(btabf$F=="unknown" & btabf$G=="unknown"),] ####line changed 
  
  message("applying family top threshold of ",topF)
  topdf<-aggregate(x = btabf[,"pident"],by=list(btabf$qseqid),FUN = max)
  colnames(topdf)<-c("qseqid","pident")
  topdf$min_pident<-topdf$pident-topF
  btabf$min_pident<-topdf[match(btabf$qseqid,topdf$qseqid),"min_pident"]
  btabf<-btabf[btabf$pident>btabf$min_pident,]
  
  if(!is.null(disabledTaxaFiles)){
    btabf<-btabf[!btabf$taxids %in% unique(c(childrenF)),]
  }
  
  btabf<-btabf[btabf$pident>fpident,]
  
  if(nrow(btabf)>0){
    
    btabf$path<-paste(btabf$K,btabf$P,btabf$C,btabf$O,btabf$F,btabf$G,btabf$S,sep = ";")
    lcaf = aggregate(btabf$path, by=list(btabf$qseqid),function(x) lca(x,sep=";"))
    colnames(lcaf)<-c("qseqid","binpath")
    lcaf<-add.unknown.lca(lcaf)
    mat<-do.call(rbind,stringr::str_split(lcaf$binpath,";"))
    lcaf<-as.data.frame(cbind(lcaf$qseqid,mat[,1],mat[,2],mat[,3],mat[,4],mat[,5],mat[,6],mat[,7]))
    colnames(lcaf)<-c("qseqid","K","P","C","O","F","G","S")
    
  } else {
    lcaf<-data.frame(matrix(nrow=1,ncol = 8))
    colnames(lcaf)<-c("qseqid","K","P","C","O","F","G","S")
  }
  
  rm(btabf)
  
  #higher-than-family-level assignments
  message("binning at higher-than-family level")
  btabhtf<-btab[btab$K!="unknown",]
  
  message("applying htf top threshold of ",topAbs)
  topdf<-aggregate(x = btabhtf[,"pident"],by=list(btabhtf$qseqid),FUN = max)
  colnames(topdf)<-c("qseqid","pident")
  topdf$min_pident<-topdf$pident-topAbs
  btabhtf$min_pident<-topdf[match(btabhtf$qseqid,topdf$qseqid),"min_pident"]
  btabhtf<-btabhtf[btabhtf$pident>btabhtf$min_pident,]
  
  btabhtf<-btabhtf[btabhtf$pident>abspident,]
  
  if(nrow(btabhtf)>0){
    
    btabhtf$path<-paste(btabhtf$K,btabhtf$P,btabhtf$C,btabhtf$O,btabhtf$F,btabhtf$G,btabhtf$S,sep = ";")
    lcahtf = aggregate(btabhtf$path, by=list(btabhtf$qseqid),function(x) lca(x,sep=";"))
    colnames(lcahtf)<-c("qseqid","binpath")
    lcahtf<-add.unknown.lca(lcahtf)
    mat<-do.call(rbind,stringr::str_split(lcahtf$binpath,";"))
    lcahtf<-as.data.frame(cbind(lcahtf$qseqid,mat[,1],mat[,2],mat[,3],mat[,4],mat[,5],mat[,6],mat[,7]))
    colnames(lcahtf)<-c("qseqid","K","P","C","O","F","G","S")
    
  } else {
    lcahtf<-data.frame(matrix(nrow=1,ncol = 8))
    colnames(lcahtf)<-c("qseqid","K","P","C","O","F","G","S")
  }
  
  rm(btabhtf)
  
  ###################################################
  #combine
  #combine
  sp_level<-lcasp[lcasp$S!="unknown",]
  g_level<-lcag[lcag$G!="unknown",]
  if(nrow(g_level)>0) g_level$S<-NA
  g_level<-g_level[!g_level$qseqid %in% sp_level$qseqid,]
  f_level<-lcaf[lcaf$F!="unknown",]
  if(nrow(f_level)>0) f_level$G<-NA
  if(nrow(f_level)>0) f_level$S<-NA
  f_level<-f_level[!f_level$qseqid %in% sp_level$qseqid,]
  f_level<-f_level[!f_level$qseqid %in% g_level$qseqid,]
  
  abs_level<-lcahtf
  if(nrow(abs_level)>0) abs_level$G<-NA
  if(nrow(abs_level)>0) abs_level$S<-NA
  if(nrow(abs_level)>0) abs_level$F<-NA
  abs_level<-abs_level[!abs_level$qseqid %in% sp_level$qseqid,]
  abs_level<-abs_level[!abs_level$qseqid %in% g_level$qseqid,]
  abs_level<-abs_level[!abs_level$qseqid %in% f_level$qseqid,]
  
  com_level<-rbind(sp_level,g_level,f_level,abs_level)
  com_level<-merge(x=qseqids, y = com_level, by = "qseqid",all.x = T)
  
  #info
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  
  write.table(x = com_level,file = out,sep="\t",quote = F,row.names = F)
  
  message(c("Complete. ",total_hits, " hits from ", total_queries," queries processed in ",t3," mins."))
  
  message("Note: If none of the hits for a BLAST query pass the binning thesholds, the results will be NA for all levels.
                 If the LCA for a query is above kingdom, e.g. cellular organisms or root, the results will be 'unknown' for all levels.
                 Queries that had no BLAST hits, or did not pass the filter.blast step will not appear in results.  ")
}

# =========================================================
# Copyright 2019-2020,  Nuno A. Fonseca (nuno dot fonseca at gmail dot com)
#
#
# This is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# if not, see <http://www.gnu.org/licenses/>.
#
#
# =========================================================
lca <- function(paths, threshold=1.0, sep=":",                
                remove.dups=FALSE,
                normalize.entries=FALSE,
                remove.trailing.nas=FALSE) {
  
  if (is.null(paths)) return(NA)
  # remove dups...  or not 
  if (remove.dups) paths<-unique(paths)
  ## Normalize netries
  if ( normalize.entries ) {
    paths <- sapply(paths,FUN=sub,pattern=" [^:]+$",replacement="",ignore.case=TRUE)
    paths <- sapply(paths,FUN=sub,pattern="([^0-9])[0-9]+$",replacement="\\1",ignore.case=TRUE)
    
  }
  
  ## remove trailing NAs
  if ( remove.trailing.nas ) {
    paths <- sapply(paths,FUN=sub,pattern="(:NA)+$",replacement="",ignore.case=TRUE)   
  }
  ## workaround to handle paths with toplevel entries only
  paths<-paste(paths,":NA-",sep="")
  v<-sapply(paths,function(l) strsplit(l,sep)[[1]])
  ##
  if ( typeof(v) == "list" ) {
    ## ensure that length of all elements is the same
    lens<-sapply(v,length)
    target<-max(lens,na.rm=TRUE)
    unilength<-function(l,tlen) {
      nnewels<-tlen-length(l)
      if (nnewels==0) return(l)
      return(append(l,rep(x="NA",nnewels)))
    }
    v<-lapply(v,unilength,tlen=target)
  } 
  
  v<-data.frame(v,check.names=FALSE,stringsAsFactors=FALSE)    
  nt<-apply(v,1,table,useNA="no")
  if (typeof(nt) == "integer") {
    ## single entry        
    return(sub(":NA-$","",x=paths[1]))
  }
  vals<-lapply(nt,function(x) x/sum(x,na.rm=TRUE))#,simplify=FALSE)
  
  w<-NA
  for ( e in 1:length(vals)) {
    w<-vals[[e]]
    if (max(w)<threshold) {
      e<-e-1
      break;
    }
  }
  if (e==0) return(NA)
  w<-vals[[e]]
  # get the full paths
  l<-names(which(w>=threshold))
  cols<-as.character(v[e,])%in%l
  v2<-v[1:e,cols,drop=FALSE]
  
  r<-unique(apply(v2,c(2),paste,sep=sep,collapse=sep))
  return(sub(":NA-$","",x=r))
}

add.unknown.lca<-function(lca.out){
  
  lca.out$binpath[is.na(lca.out$binpath)]<-"unknown;unknown;unknown;unknown;unknown;unknown;unknown"
  lca.out$binpath[stringr::str_count(lca.out$binpath,";")==5]<-paste0(lca.out$binpath[stringr::str_count(lca.out$binpath,";")==5],";unknown")
  lca.out$binpath[stringr::str_count(lca.out$binpath,";")==4]<-paste0(lca.out$binpath[stringr::str_count(lca.out$binpath,";")==4],";unknown;unknown")
  lca.out$binpath[stringr::str_count(lca.out$binpath,";")==3]<-paste0(lca.out$binpath[stringr::str_count(lca.out$binpath,";")==3],";unknown;unknown;unknown")
  lca.out$binpath[stringr::str_count(lca.out$binpath,";")==2]<-paste0(lca.out$binpath[stringr::str_count(lca.out$binpath,";")==2],";unknown;unknown;unknown;unknown")
  lca.out$binpath[stringr::str_count(lca.out$binpath,";")==1]<-paste0(lca.out$binpath[stringr::str_count(lca.out$binpath,";")==1],";unknown;unknown;unknown;unknown;unknown")
  lca.out$binpath[stringr::str_count(lca.out$binpath,";")==0]<-paste0(lca.out$binpath[stringr::str_count(lca.out$binpath,";")==0],";unknown;unknown;unknown;unknown;unknown;unknown")
  
  return(lca.out)
}

report.blast.maxmin<-function(blastfile,pident_col="V4",qseqid_col="V1"){
  
  btab<-data.table::fread(blastfile,sep="\t",data.table = F)
  
  tophits<-aggregate(btab[,pident_col],by=list(btab[,qseqid_col]),FUN=max)
  colnames(tophits)<-c("qseqid","pident")
  tophits$minmax<-"max"
  
  minhits<-aggregate(btab[,pident_col],by=list(btab[,qseqid_col]),FUN=min)
  colnames(minhits)<-c("qseqid","pident")
  minhits$minmax<-"min"
  
  tophits<-rbind(tophits,minhits)
  
  return(tophits)
}

compare.sample.reads<-function(taxatab1,taxatab2,rm.nohits=T,rm.NAs=T){
  # plot and compare reads / sample
  
  if(rm.nohits){
    message("Removing no_hits")
    a<-sum(taxatab1[,-1])
    b<-sum(taxatab2[,-1])
    taxatab1<-taxatab1[taxatab1$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]
    taxatab2<-taxatab2[taxatab2$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]
    message(a-sum(taxatab1[,-1]), " reads removed from taxatab1")
    message(b-sum(taxatab2[,-1]), " reads removed from taxatab2")
  }
  
  if(rm.NAs){
    a<-sum(taxatab1[,-1])
    b<-sum(taxatab2[,-1])
    message("Removing NAs - sequences that had blast hits but were not assigned to any taxon")
    taxatab1<-taxatab1[taxatab1$taxon!="NA;NA;NA;NA;NA;NA;NA",]
    taxatab1<-taxatab1[taxatab1$taxon!="unknown;unknown;unknown;unknown;NA;NA;NA",]
    
    taxatab2<-taxatab2[taxatab2$taxon!="NA;NA;NA;NA;NA;NA;NA",]
    taxatab2<-taxatab2[taxatab2$taxon!="unknown;unknown;unknown;unknown;NA;NA;NA",]
    
    message(a-sum(taxatab1[,-1]), " reads removed from taxatab1")
    message(b-sum(taxatab2[,-1]), " reads removed from taxatab2")
  }
  
  results<-data.frame(matrix(ncol = 3, nrow = 6))
  colnames(results)<-c("stats",name.taxatab1,name.taxatab2)
  results$stats<-c("min.sample.reads","max.sample.reads","mean.sample.reads","SD.sample.reads","n.samples","total.reads")
  
  reads1<-data.frame(reads=colSums(rm.0readtaxSam(taxatab1)[,-1]),taxatab=name.taxatab1)
  reads2<-data.frame(reads=colSums(rm.0readtaxSam(taxatab2)[,-1]),taxatab=name.taxatab2)
  
  results[1,name.taxatab1]<-range(reads1$reads)[1]
  results[2,name.taxatab1]<-range(reads1$reads)[2]
  results[1,name.taxatab2]<-range(reads2$reads)[1]
  results[2,name.taxatab2]<-range(reads2$reads)[2]
  results[3,name.taxatab1]<-mean(reads1$reads)
  results[3,name.taxatab2]<-mean(reads2$reads)
  results[4,name.taxatab1]<-sd(reads1$reads)
  results[4,name.taxatab2]<-sd(reads2$reads)
  results[5,name.taxatab1]<-nrow(reads1)
  results[5,name.taxatab2]<-nrow(reads2)
  results[6,name.taxatab1]<-sum(reads1$reads)
  results[6,name.taxatab2]<-sum(reads2$reads)
  
  reads1and2<-rbind(reads1,reads2)
  
  reads.plot<-ggplot(reads1and2,aes(x=taxatab,y=reads))+geom_boxplot() + scale_y_continuous(labels = scales::comma) + ylab("reads/sample") +
    xlab("") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  m1<-lm(formula = reads~taxatab,data = reads1and2)
  mout<-capture.output(interp.lm(m1))
  
  out<-list(results,reads.plot,mout)
  
  return(out)
}

chisq.taxatab.ranks<-function(taxatab1,taxatab2,name.taxatab1,name.taxatab2){
  
  message("Reminder: comparison includes levels below, e.g. genus means reads attaining genus or lower")
  
  ranks<-c("htf","family","genus" ,"species")
  
  rank.taxatab1<-stats.by.rank(taxatab1)
  rank.taxatab1$taxatab<-name.taxatab1
  rank.taxatab2<-stats.by.rank(taxatab2)
  rank.taxatab2$taxatab<-name.taxatab2
  
  rank.taxatab1and2<-rbind(rank.taxatab1,rank.taxatab2)
  rank.taxatab1and2$rank<-factor(rank.taxatab1and2$rank,levels = ranks)
  
  #compare proportions
  totalx<-aggregate(rank.taxatab1and2$reads,by=list(rank.taxatab1and2$taxatab), FUN=sum)
  
  output<-data.frame(matrix(nrow = length(ranks),ncol=6))
  colnames(output)<-c("rank",name.taxatab1,name.taxatab2,"p-value","chisq","df")
  output$rank<-ranks
  output$df=1
  
  for(i in 1:length(ranks)){
    message(ranks[i])
    temptable<-rank.taxatab1and2
    if(ranks[i]=="species")  {
      propx<-temptable[temptable$rank==ranks[i],c("reads","taxatab")]
      colnames(propx)<-c("x","Group.1")
    }
    if(ranks[i]=="genus") {
      temptable[temptable$rank=="species","rank"]<-"genus"
      temptable<-temptable[temptable$rank=="genus",c("reads","taxatab")]
      propx<-aggregate(temptable$reads,by=list(temptable$taxatab),FUN=sum)
    }
    if(ranks[i]=="family") {
      temptable[temptable$rank=="species","rank"]<-"family"
      temptable[temptable$rank=="genus","rank"]<-"family"
      temptable<-temptable[temptable$rank=="family",c("reads","taxatab")]
      propx<-aggregate(temptable$reads,by=list(temptable$taxatab),FUN=sum)
    }
    if(ranks[i]=="htf") {
      temptable[temptable$rank=="species","rank"]<-"htf"
      temptable[temptable$rank=="genus","rank"]<-"htf"
      temptable[temptable$rank=="family","rank"]<-"htf"
      temptable<-temptable[temptable$rank=="htf",c("reads","taxatab")]
      propx<-aggregate(temptable$reads,by=list(temptable$taxatab),FUN=sum)
    }
    
    merged<-merge(propx,totalx,by = "Group.1")
    message("prop1=",merged[1,1],"; prop2=",merged[2,1])
    a<-prop.test(merged[,2],merged[,3])
    print(a)
    
    output[output$rank==ranks[i],"chisq"]<-a$statistic
    output[output$rank==ranks[i],"p-value"]<-a$p.value
    output[output$rank==ranks[i],merged[1,1]]<-round(a$estimate[1],digits=3)
    output[output$rank==ranks[i],merged[2,1]]<-round(a$estimate[2],digits=3)
    
  }
  return(output)
}

blast.maxhits.2.files.compare<-function(blastfile1,blastfile2,name.blastfile1,name.blastfile2,pident_col="V5",qseqid_col="V1",cutoff){
  
  message("Comparing proportion of hits above ", cutoff, "% identity")
  
  btab1<-data.table::fread(blastfile1,sep="\t",data.table = F)
  btab2<-data.table::fread(blastfile2,sep="\t",data.table = F)
  
  tophits1<-aggregate(btab1[,pident_col],by=list(btab1[,qseqid_col]),FUN=max)
  colnames(tophits1)<-c("qseqid","pident")
  tophits1$blastfile<-name.blastfile1
  tophits2<-aggregate(btab2[,pident_col],by=list(btab2[,qseqid_col]),FUN=max)
  colnames(tophits2)<-c("qseqid","pident")
  tophits2$blastfile<-name.blastfile2
  
  a1<-as.data.frame(as.factor(quantile(tophits1$pident, c(.01,.1, .2, .3, .4, .5, .6, .7, .8,.9, .95, .99, 1)))) 
  a2<-as.data.frame(as.factor(quantile(tophits2$pident, c(.01,.1, .2, .3, .4, .5, .6, .7, .8,.9, .95, .99, 1)))) 
  
  colnames(a1)<-"top_hits"
  a1$quantile<-rownames(a1)
  a1$quantile <- factor(a1$quantile, levels = a1$quantile)
  a1$top_hits<-round(as.numeric(as.character(a1$top_hits)),digits = 1)
  a1$blastfile<-name.blastfile1
  
  colnames(a2)<-"top_hits"
  a2$quantile<-rownames(a2)
  a2$quantile <- factor(a2$quantile, levels = a2$quantile)
  a2$top_hits<-round(as.numeric(as.character(a2$top_hits)),digits = 1)
  a2$blastfile<-name.blastfile2
  
  a1anda2<-rbind(a1,a2)
  
  plot1<-ggplot(data = a1anda2,aes(x=quantile,y=top_hits,color=blastfile,group=blastfile)) + geom_line() +  
    xlab("% of queries with top hit above y") + ylab("% identity") +
    theme(axis.text.x=element_text(size=8,angle=45, hjust=1),legend.title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  allhits<-rbind(tophits1,tophits2)
  
  totalx<-aggregate(allhits$pident,by=list(allhits$blastfile), FUN=sum)
  
  propx<-aggregate(allhits[allhits$pident>cutoff,]$pident,by=list(allhits[allhits$pident>cutoff,]$blastfile), FUN=sum)
  
  merged<-merge(propx,totalx,by = "Group.1")
  message("prop1=",merged[1,1],"; prop2=",merged[2,1])
  a<-prop.test(merged[,2],merged[,3])
  
  plot2<-ggplot(allhits,aes(x=blastfile,y=pident)) + geom_violin() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    ylab("% identity of top hit") + xlab("")
  
  output<-list(plot1,plot2,a)
}

plot.taxatab.rank.props<-function(taxatab.list,y="reads",grayscale=F){
  
  #y can be "reads, "dxns"
  
  message("Reminder: must be a named list, see split.taxatab.by.facets")
  
  message("Not including no_hits")
  
  taxatab.list2<-list()
  for(i in 1:length(taxatab.list)){
    taxatab.list2[[i]]<-taxatab.list[[i]][taxatab.list[[i]]$taxon!="no_hits;no_hits;no_hits;no_hits;no_hits;no_hits;no_hits",]
    taxatab.list2[[i]]<-stats.by.rank(taxatab.list[[i]],grouphtf = F)
    taxatab.list2[[i]]$taxatab<-names(taxatab.list)[i]
  }
  
  rank.taxatabs<-do.call(rbind,taxatab.list2)
  
  colnames(rank.taxatabs)<-gsub("rank","Rank",colnames(rank.taxatabs))
  
  #rank.taxatabs$Rank<-gsub("htf","above family",rank.taxatabs$Rank)
  
  rank.taxatabs$Rank<-factor(rank.taxatabs$Rank,levels = c("above_kingdom","kingdom","phylum", "class","order","family","genus","species" ))
  
  if(y=="reads"){
    a<-ggplot(rank.taxatabs,aes(x=taxatab,y=reads,fill=Rank))+geom_bar(stat="identity",position="fill")+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x=element_text(size=8,angle=45, hjust=1))+
      ylab("Proportion of reads") + xlab("") + scale_fill_manual(values = MyCols)
  }
  
  if(y=="dxns"){
    a<-ggplot(rank.taxatabs,aes(x=taxatab,y=dxns,fill=Rank))+geom_bar(stat="identity",position="fill")+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text.x=element_text(size=8,angle=45, hjust=1))+
      ylab("Proportion of detections") + xlab("") + scale_fill_manual(values = MyCols)
  }
  if(grayscale==T) a<-a+ scale_fill_grey()
  
  a
  
}

taxatab.heatmap<-function(taxatab,master_sheet,group.by="ss_sample_id",values="nreads",current.grouping = "ss_sample_id",colour.bar=NULL,
                          facetcol=NULL,inc.values=T,tidy.taxon.names="order",taxafontsize=10,colfontsize=10){
  
  if(!is.data.frame(master_sheet)) stop("master_sheet should be a data frame")
  
  if(values!="ndxns" & values!="nreads" & values!="dxn") stop("values must be one of nreads, ndxns, dxn")
  
  require("ComplexHeatmap")
  
  taxatab2<-taxatab
  
  if(values=="ndxns" | values=="dxn") taxatab2<-binarise.taxatab(taxatab2,t=T)
  
  if(!is.null(facetcol)){
    #if facetcol is not null, then need to split taxatabs by facetcol
    #first split mastersheet
    split.list<-list()
    for(i in 1:length(facetcol)){
      split.list[[i]]<-as.factor(master_sheet[,facetcol[i]])
    }
    
    ms.split<-split(master_sheet,split.list)
    
    #then split taxatab, removing any that have no samples for that facet
    taxatab.list<-list()
    for(i in 1:length(ms.split)){
      taxatabtemp<-cbind(taxon=taxatab2$taxon,taxatab2[colnames(taxatab2) %in% ms.split[[i]][,current.grouping]])
      if(ncol(taxatabtemp)>1) taxatab.list[[i]]<-taxatabtemp else taxatab.list[[i]]<-NULL
    }
    
    #remove taxatables and split.list if NULL (accountign that taxtable.list may be shorter than ms.list)
    ms.split<-ms.split[c(1:length(taxatab.list))]
    if(length(which(sapply(taxatab.list, is.null)))>0) {
        ms.split<-ms.split[-which(sapply(taxatab.list, is.null))]
        taxatab.list<-taxatab.list[-which(sapply(taxatab.list, is.null))]
    }
    
    } else {
    taxatab.list<-list()
    taxatab.list[[1]]<-taxatab2
    ms.split<-list()
    ms.split[[1]]<-master_sheet
  }
  
  
  
  #then group by group.by
  for(i in 1:length(taxatab.list)){
    taxatab.list[[i]]<-sumreps(taxatab = taxatab.list[[i]],ms_ss = ms.split[[i]],grouping = group.by,discard = F,
                               current.grouping = current.grouping)
  }
  
  #then make all taxa present in all taxatables, in same order
  all.taxa<-unique(unlist(lapply(taxatab.list,function(x) x[,1])))
  for(i in 1:length(taxatab.list)){
    missing.taxa<-all.taxa[!all.taxa %in% taxatab.list[[i]]$taxon]
    if(length(missing.taxa)>0){
      fake.taxatab<-taxatab.list[[i]][1:length(missing.taxa),]
      fake.taxatab$taxon<-missing.taxa
      fake.taxatab[,-1]=0
      taxatab.list[[i]]<-rbind(taxatab.list[[i]],fake.taxatab)
      taxatab.list[[i]]<-taxatab.list[[i]][order(taxatab.list[[i]]$taxon),]
    } else taxatab.list[[i]]<-taxatab.list[[i]][order(taxatab.list[[i]]$taxon),]
  }
  
  #sort by rank
  for(i in 1:length(taxatab.list)){
    ranks<-factor(bas.get.ranks(taxatab.list[[i]]),levels = c("species","genus","family","htf"))
    taxatab.list[[i]]<-taxatab.list[[i]][order(ranks),]
    ranks<-factor(bas.get.ranks(taxatab.list[[i]]),levels = c("species","genus","family","htf"))
  }
  
  #ditto for samples
  all.samples<-unique(unlist(lapply(taxatab.list,function(x) colnames(x[,-1]))))
  for(i in 1:length(taxatab.list)){
    missing.samples<-all.samples[!all.samples %in% colnames(taxatab.list[[i]][,-1])]
    if(length(missing.samples)>0){
      a<-data.frame(matrix(ncol=length(missing.samples),nrow = nrow(taxatab.list[[i]])))
      colnames(a)<-missing.samples
      a[,]=0
      taxatab.list[[i]]<-cbind(taxatab.list[[i]],a)
      taxatab.list[[i]]<-taxatab.list[[i]][,order(colnames(taxatab.list[[i]]))]
      taxatab.list[[i]]<-taxatab.list[[i]]  %>% select(taxon,everything())
    } else {
      taxatab.list[[i]]<-taxatab.list[[i]][,order(colnames(taxatab.list[[i]]))]
      taxatab.list[[i]]<-taxatab.list[[i]]  %>% select(taxon,everything())
    }
  }
  
  #colour bar
  if(!is.null(colour.bar)) {
    
    if(length(colour.bar)>3) stop("More than 3 colour.bars requested, limit is 3")
    
    outgroups<-list()
    outcols<-list()
    
    for(i in 1:length(colour.bar)){
    
      colour.bar.groups<-master_sheet[match(colnames(taxatab.list[[1]][,-1,drop=F]),master_sheet[,group.by]),colour.bar[i]]
      
      if(sum(is.na(colour.bar.groups))>0) {
        message("some colour bar groups were NA, changing to unknown")
        colour.bar.groups[is.na(colour.bar.groups)]<-"unknown"
      }
      
      if(i==1){
        if(length(unique(colour.bar.groups))<29) { 
          MyColsx<-MyCols
        } else MyColsx<-rainbow(length(unique(colour.bar.groups)))
        colnum1<-length(unique(colour.bar.groups))
      }
      
      if(i==2){
        if(length(unique(colour.bar.groups))+colnum1<29) { 
          MyColsx<-MyCols
          MyColsx<-MyColsx[(colnum1+1):length(MyColsx)]
        } else {
          MyColsx<-rainbow(colnum1+length(unique(colour.bar.groups)))
          MyColsx<-MyColsx[(colnum1+1):length(MyColsx)]
        }
        colnum2<-length(unique(colour.bar.groups))
      }
      
      if(i==3){
        if(length(unique(colour.bar.groups))+colnum1+colnum2<29) { 
          MyColsx<-MyCols
          MyColsx<-MyColsx[(colnum1+colnum2+1):length(MyColsx)]
        } else {
          MyColsx<-rainbow(colnum1+colnum2+length(unique(colour.bar.groups)))
          MyColsx<-MyColsx[(colnum1+colnum2+1):length(MyColsx)]
        }
      }
      
      colours.list<-MyColsx[as.numeric(as.factor(colour.bar.groups))]
      names(colours.list)<-as.factor(colour.bar.groups)
      colours.list<-list(colours.list)
      
      if(i==1) names(colours.list)<-"bar"
      if(i==2) names(colours.list)<-"foo"
      if(i==3) names(colours.list)<-"tang"
      
      outgroups[[i]]<-colour.bar.groups
      outcols[[i]]<-colours.list
    }
    
    if(length(colour.bar)==1) {
    
      ha = ComplexHeatmap::HeatmapAnnotation(bar=colour.bar.groups,show_annotation_name = F,col = colours.list,
                                           annotation_legend_param = list(bar = list(title = colour.bar)))
                                           }
    
    if(length(colour.bar)==2) {
      
      ha = ComplexHeatmap::HeatmapAnnotation(bar=outgroups[[1]],foo=outgroups[[2]],show_annotation_name = F,
                                             col = c(outcols[[1]],outcols[[2]]),
                                             annotation_legend_param = list(bar = list(title = colour.bar[1])
                                             ,foo=list(title=colour.bar[2])))
    }
    
    if(length(colour.bar)==3) {
      
      ha = ComplexHeatmap::HeatmapAnnotation(bar=outgroups[[1]],foo=outgroups[[2]],tang=outgroups[[3]],show_annotation_name = F,
                                             col = c(outcols[[1]],outcols[[2]],outcols[[3]]),
                                             annotation_legend_param = list(bar = list(title = colour.bar[1])
                                                                            ,foo=list(title=colour.bar[2])
                                                                            ,tang=list(title=colour.bar[3])))
    }
    
    
  } else ha=NULL
  
  
  taxatab.list2<-list()
  hm.list<-list()
  for(i in 1:length(taxatab.list)){
    
    if(values=="dxn") taxatab.list2[[i]]<-binarise.taxatab(taxatab.list[[i]],t=T) else taxatab.list2[[i]]<-taxatab.list[[i]]
    
    if(!is.null(tidy.taxon.names)) taxatab.list2[[i]]<-tidy.taxon(taxatab=taxatab.list2[[i]],rm.trailing.NA=T,
                                                                  rm.preceeding.above=tidy.taxon.names)
    
    #convert to matrix
    rownames(taxatab.list2[[i]])<-make.unique(taxatab.list2[[i]]$taxon)
    taxatab.list2[[i]]$taxon=NULL
    taxatab.list2[[i]]<-as.matrix(taxatab.list2[[i]])
  }
   
  
    
    if(values=="dxn"){
      
      for(i in 1:length(taxatab.list2)){
      
      #if detected in all samples, need to change a little. 
      if(mean(taxatab.list2[[i]])==1) {
        coldxns="black"
        levels=1
        labels="Detected" 
      } else {
        coldxns=c("grey90","black")
        levels=c(1,0)
        labels=c("Detected","Not Detected") 
      }
      
      CurrentData<-round(taxatab.list2[[i]],digits = 0)
      
      if(inc.values) {values_func<-local({
        CurrentData = CurrentData
        function(a, b, x, y, width, height, fill) {
          if(CurrentData[b, a] > 0) grid.text(sprintf("%.0f", CurrentData[b,a]), x, y, gp = gpar(fontsize = 8,col="red"))}})
      } else values_func<-NULL
      
      if(!is.null(colour.bar)) {
        ordercols<-with(as.data.frame(colnames(taxatab.list2[[i]])),order(colour.bar.groups,as.factor(colour.bar.groups)))
      } else ordercols<-NULL
      
      hm.list[[i]]<-ComplexHeatmap::Heatmap(taxatab.list2[[i]],
                                            bottom_annotation=ha, 
                                            col=coldxns,name = " ",
                                            #column_title = plot_title,
                                            rect_gp = grid::gpar(col = "white", lwd = 2),
                                            border = FALSE,
                                            cluster_rows = FALSE,
                                            cluster_columns = FALSE,
                                            column_names_rot = 45,
                                            row_names_side = "left",
                                            row_names_max_width = max_text_width(rownames(taxatab.list2[[i]]), gp = grid::gpar(fontsize = 12)),
                                            row_names_gp = gpar(fontsize = taxafontsize),
                                            column_names_gp = gpar(fontsize = colfontsize),
                                            column_order = ordercols,
                                            heatmap_legend_param = list(at=levels,labels = labels),
                                            column_title = names(ms.split)[i],
                                            column_title_gp = gpar(fontsize = 10),
                                            cell_fun = values_func,
                                            row_split = ranks,
                                            #column_split=newfacets2
      )
      }
    }
    
    if(values=="ndxns") {
      
      for(i in 1:length(taxatab.list2)){
      
      CurrentData<-round(taxatab.list2[[i]],digits = 0)
      
      if(inc.values) {values_func<-local({
        CurrentData = CurrentData
        function(a, b, x, y, width, height, fill) {
          if(CurrentData[b, a] > 0) grid.text(sprintf("%.0f", CurrentData[b,a]), x, y, gp = gpar(fontsize = 8))}})
      } else values_func<-NULL
      
      if(!is.null(colour.bar)) {
        ordercols<-with(as.data.frame(colnames(taxatab.list2[[i]])),order(colour.bar.groups,as.factor(colour.bar.groups)))
      } else ordercols<-NULL
      
      hm.list[[i]]<-ComplexHeatmap::Heatmap(taxatab.list2[[i]],
                                            bottom_annotation=ha, 
                                            col=circlize::colorRamp2(c(0, max(unlist(lapply(taxatab.list2,FUN=max)))), c("grey90", "red")),
                                            name = "Number of detections",
                                            #column_title = plot_title,
                                            rect_gp = grid::gpar(col = "white", lwd = 2),
                                            border = FALSE,
                                            cluster_rows = FALSE,
                                            cluster_columns = FALSE,
                                            column_names_rot = 45,
                                            row_names_side = "left",
                                            row_names_max_width = max_text_width(rownames(taxatab.list2[[i]]), gp = grid::gpar(fontsize = 12)),
                                            row_names_gp = gpar(fontsize = taxafontsize),
                                            column_names_gp = gpar(fontsize = colfontsize),
                                            column_order = ordercols,
                                            column_title = names(ms.split)[i],
                                            column_title_gp = gpar(fontsize = 10),
                                            cell_fun = values_func,
                                            row_split = ranks
                                            #column_split=newfacets2
      )
      }
    }
    
    if(values=="nreads")  {
      
      for(i in 1:length(taxatab.list2)){
        
      #log scale
      taxatab.list2[[i]][taxatab.list2[[i]]==0] = NaN
      taxatab.list2[[i]][taxatab.list2[[i]]==1] = 1.01
      taxatab.list2[[i]]<-log(taxatab.list2[[i]],10)  
      taxatab.list2[[i]][is.nan(taxatab.list2[[i]])] = 0
      
      }
      
      for(i in 1:length(taxatab.list2)){
      
      CurrentData<-round(taxatab.list2[[i]],digits = 7)
      
      if(inc.values) {values_func<-local({
        CurrentData = CurrentData
        function(a, b, x, y, width, height, fill) {
          if(CurrentData[b, a] > 0) grid.text(sprintf("%.1f", CurrentData[b,a]), x, y, gp = gpar(fontsize = 8),rot = 90)}})
      } else values_func<-NULL
      
      if(!is.null(colour.bar)) {
        ordercols<-with(as.data.frame(colnames(taxatab.list2[[i]])),order(colour.bar.groups,as.factor(colour.bar.groups)))
      } else ordercols<-NULL
      
      hm.list[[i]]<-ComplexHeatmap::Heatmap(taxatab.list2[[i]],
                                            bottom_annotation=ha, 
                                            col=circlize::colorRamp2(c(0, max(unlist(lapply(taxatab.list2,FUN=max)))), 
                                                                     c("grey90", "red")),
                                            name = "Reads (log10)",
                                            #column_title = plot_title,
                                            rect_gp = grid::gpar(col = "white", lwd = 2),
                                            border = FALSE,
                                            cluster_rows = FALSE,
                                            cluster_columns = FALSE,
                                            column_names_rot = 45,
                                            row_names_side = "left",
                                            row_names_max_width = max_text_width(rownames(taxatab.list2[[i]]), 
                                                                                 gp = grid::gpar(fontsize = 12)),
                                            row_names_gp = gpar(fontsize = taxafontsize),
                                            column_names_gp = gpar(fontsize = colfontsize),
                                            column_order = ordercols,
                                            column_title = names(ms.split)[i],
                                            column_title_gp = gpar(fontsize = 10),
                                            cell_fun = values_func,
                                            row_split = ranks
                                            #column_split=newfacets2
                                            )
    }
    }
  
  for(i in 1:length(hm.list)){
    if(i==1) hm.out<-hm.list[[i]]
    if(i>1) hm.out<-hm.out + hm.list[[i]]
  }
  
  plot_title<-paste("x axis=",group.by,"; values=",values,";colour.bar=",paste(colour.bar,collapse=","),"; facetcol=",
                    paste(facetcol,collapse = ","))
  
  draw(hm.out,column_title=plot_title)
  
}


#print colours in a vector
print.cols<-function(vector){
  cols <- function(a) image(1:length(vector), 1, as.matrix(1:length(vector)), col=vector, axes=FALSE , xlab="", ylab="")
  a<-(1:length(vector))
  cols(a)
}

split.taxatab.by.facets<-function(taxatab,master_sheet,facetcol,current.grouping="ss_sample_id"){
  
  split.list<-list()
  for(i in 1:length(facetcol)){
    split.list[[i]]<-as.factor(master_sheet[,facetcol[i]])
  }
  
  ms.split<-split(master_sheet,split.list)
  
  #then split taxatab
  taxatab.list<-list()
  for(i in 1:length(ms.split)){
    taxatab.list[[i]]<-cbind(taxon=taxatab$taxon,taxatab[colnames(taxatab) %in% ms.split[[i]][,current.grouping]])
  }
  
  names(taxatab.list)<-names(ms.split)
  
  return(taxatab.list)
}


tidy.taxon<-function(taxatab,rm.trailing.NA=T,rm.preceeding.above="family"){
  
  if(!is.null(rm.preceeding.above)) {
    taxasplit<-do.call(rbind, stringr::str_split(taxatab$taxon,";"))
    if(rm.preceeding.above=="phylum") taxasplit<-taxasplit[,c(2:7)]
    if(rm.preceeding.above=="class") taxasplit<-taxasplit[,c(3:7)]
    if(rm.preceeding.above=="order") taxasplit<-taxasplit[,c(4:7)]
    if(rm.preceeding.above=="family") taxasplit<-taxasplit[,c(5:7)]
    if(rm.preceeding.above=="genus") taxasplit<-taxasplit[,c(6:7)]
    if(rm.preceeding.above=="species") taxasplit<-taxasplit[,7]
    
    taxatab$taxon<-do.call(paste,c(as.data.frame(taxasplit),sep=";"))
  }
  
  if(rm.trailing.NA) for(i in 1:7) taxatab$taxon<-gsub(";NA$","",taxatab$taxon)
  
  return(taxatab)
  
}
#' Perform BLAST 
#' @param infasta Name of fasta file to BLAST. Must be in current directory.
#' @param refdb The BLAST database to query.
#' @param blast_exec Path to blastn executable, if not using the blastn in current environment.
#' @param wait logical indicating whether R should wait for BLASTs to finish. Default=T. Recommended.
#' @param taxidlimit A vector of NCBI taxids to restrict the BLAST to. 
#' @param inverse logical to indicate whether the search should be restricted to the inverse of taxidlimit. Only required if using taxidlimit.
#' @param ncbiTaxDir Path to the folder containing the NCBI taxonomy dump files (see getNCBITaxonomyDump). Only required if using taxidlimit.
#' @param overWrite logical indicating whether to overwrite any existing BLAST files of the same name, if they exist.
#' @param out Name of the file for the BLAST results, Default = paste0(gsub(".fasta", ".blast.txt",infasta))
#' @param options A vector with flags (including dashes) and a single item for each
#'   Default=c("-outfmt", "6 qseqid evalue staxid pident qcovs","-num_threads", "64", "-max_target_seqs", 100, 
#'   "-max_hsps","1""-word_size", 6, "-perc_identity", 50, "-qcov_hsp_perc", 90,
#'   "-gapopen", 0, "-gapextend", 2, "-reward", 1, "-penalty", -1)
#' @return 
#' @examples
blast.min.bas2<-function(infasta,refdb,blast_exec="blastn",wait=T,taxidlimit=NULL,inverse=F,ncbiTaxDir=NULL,overWrite=F,out=NULL,
                         opts=c("-task","blastn","-outfmt", "6 qseqid evalue pident qcovs saccver staxid ssciname sseq","-num_threads", 64,
                                "-max_target_seqs", 100, "-max_hsps",1,"-word_size", 11, "-perc_identity", 50,
                                "-qcov_hsp_perc", 98, "-gapopen", 0, "-gapextend", 2, "-reward", 1, "-penalty", -1)){
  
  t1<-Sys.time()
  
  require(processx)
  
  if(!is.null(taxidlimit)) if(is.null(ncbiTaxDir)) stop("to use taxidlimit, ncbiTaxDir must be supplied")
  if(is.null(out)) out<-paste0(gsub(".fasta", ".blast.txt",infasta))
  
  outdircheck<-
  
  if(overWrite==F) if(file.exists(out)) stop("The following file already exists ", out, "Use overWrite=T to overwrite")
  
  if(!is.null(taxidlimit)){
    
    h<-list()
    
    #generate children of taxids and store in file
    taxid.list<-list()
    
    taxids_fileA<-paste0("taxids",as.numeric(Sys.time()),".txt")
    
    for(i in 1:length(taxidlimit)){
      system2(command = "taxonkit",args = c("list", "--ids", taxidlimit[i], "--indent", '""',"--data-dir",ncbiTaxDir)
              ,wait=T,stdout = taxids_fileA)
      taxid.list[i]<-read.table(taxids_fileA)
    }
    
    write.table(unlist(taxid.list),taxids_fileA,row.names = F,quote = F,col.names = F,)
    
    #change options to include taxidlimit
    opts<-c(opts,"-taxidlist",taxids_fileA)
    
    if(inverse) opts<-gsub("-taxidlimit","-negative_taxids",opts)
  }
  
  #run BLAST  
  
  error.log.file<-paste0("blast.error.temp.processx.file",as.numeric(Sys.time()),".txt")
  
  h<-process$new(command = blast_exec, args=c("-query", infasta, "-db",refdb,opts, "-out", out),echo_cmd = T,
                 stderr = error.log.file)
  
  Sys.sleep(time = 2)
  
  #report PID
  message(paste("PID:",h$get_pid()))
  
  #check immediate exit status
  exits<-h$get_exit_status()
  
  if(1 %in% exits){
    message("************
             There was a problem with ", infasta[match(1,exits)], ", aborting blast
             ************")
    print(grep("Error",readLines(error.log.file),value = T))
    
    h$kill()
  }
  
  if(wait==T){
    h$wait()
    exits<-h$get_exit_status()
    message(readLines(error.log.file))
    message("exit_status=",exits)
    file.remove(error.log.file)
    if(!is.null(taxidlimit)) file.remove(taxids_fileA)
  }
  
  headers<-paste0(paste("'1i",paste(unlist(strsplit(opts[match("-outfmt",opts)+1]," "))[-1],collapse = "\t"),collapse = "\t"),"'")
  
  system2("sed",c("-i", headers, out),wait = T)
  
  t2<-Sys.time()
  t3<-round(difftime(t2,t1,units = "mins"),digits = 2)
  
  message(c("All blasts complete in ",t3," mins."))
  
  return(h)
}

filter.blast3<-function(blastfile,ncbiTaxDir,out,rm.unclassified=T){
  
  message("Reminder, filter.blast3 does not include 'top' threshold")
  
  if(is.null(ncbiTaxDir)) stop("ncbiTaxDir not specified")
  if(is.null(out)) stop("out not specified")
  
  message("reading blast results")
  btab<-data.table::fread(file = blastfile,sep = "\t",header = T,data.table = F)
  
  if(!"staxid" %in% colnames(btab)) stop("No column called staxid")

  #add lineage to results
  message("adding taxonomic lineages")
  btab$taxids<-btab$staxid #add.lineage.df requires this colname
  btab<-add.lineage.df(btab,ncbiTaxDir)
  
  btab$K<-gsub(" ","_",btab$K)
  btab$P<-gsub(" ","_",btab$P)
  btab$C<-gsub(" ","_",btab$C)
  btab$O<-gsub(" ","_",btab$O)
  btab$F<-gsub(" ","_",btab$F)
  btab$G<-gsub(" ","_",btab$G)
  
  #remove crappy hits 
  if(rm.unclassified==T){
    #1. btab$S contains uncultured
    message("Removing species containing the terms: uncultured, environmental, 
            unidentified,fungal, eukaryote or unclassified")
    if(length(grep("uncultured",btab$S,ignore.case = T))>0) btab<-btab[-grep("uncultured",btab$S,ignore.case = T),]
    if(length(grep("environmental",btab$S,ignore.case = T))>0) btab<-btab[-grep("environmental",btab$S,ignore.case = T),]
    if(length(grep("unclassified",btab$S,ignore.case = T))>0) btab<-btab[-grep("unclassified",btab$S,ignore.case = T),]
    if(length(grep("unidentified",btab$S,ignore.case = T))>0) btab<-btab[-grep("unidentified",btab$S,ignore.case = T),]
    if(length(grep("fungal ",btab$S,ignore.case = T))>0) btab<-btab[-grep("fungal ",btab$S,ignore.case = T),]
    if(length(grep("eukaryote",btab$S,ignore.case = T))>0) btab<-btab[-grep("eukaryote",btab$S,ignore.case = T),]
    if(length(grep("synthetic",btab$S,ignore.case = T))>0) btab<-btab[-grep("synthetic",btab$S,ignore.case = T),]
  }
  
  write.table(x = btab,file = out,sep="\t",quote = F,row.names = F)
}

bas.get.flagged.ids.google<-function(out){
  
  sheeturl<-"https://docs.google.com/spreadsheets/d/1wgFzeM3YLyT-77mz0vANHJY6PCBIRrMeLFbrg690uso/edit#gid=0"
  
  url2<-stringr::str_split(sheeturl,"/d/")[[1]][2]
  url2<-stringr::str_split(url2,"/")[[1]][1]
  ss_info<-googlesheets4::read_sheet(ss = url2)
  write.table(ss_info$UNIQUE,file = out,row.names = F,quote=F,col.names = F)
}

threshold.bin.blast2<-function(df,qseqidcol="origseqid",qseqcol="origseq",TaxlevelTest="G",taxidcol="origtaxid",qpath="origpath"){
  #TaxlevelTest can be G or F
  if(!qseqidcol %in% colnames(df)) stop(qseqidcol," missing")
  if(!qseqcol %in% colnames(df)) stop(qseqcol," missing")
  if(!taxidcol %in% colnames(df)) stop(taxidcol," missing")
  if(!qpath %in% colnames(df)) stop(qpath," missing")
  
  if(TaxlevelTest!="F" & TaxlevelTest!="G") stop("TaxlevelTest must be F or G")
  
  colnames(df)<-gsub(qseqidcol,"qseqid",colnames(df))
  colnames(df)<-gsub(qseqcol,"seq.text",colnames(df))
  colnames(df)<-gsub(taxidcol,"taxids",colnames(df))
  colnames(df)<-gsub(qpath,"qpath",colnames(df))
  
  df$seq.name<-paste0(df$qseqid," taxid=",df$taxids,";") 
  
  lineage<-as.data.frame(stringr::str_split(df$qpath,";",simplify = T))
  colnames(lineage)<-c("K","P","C","O","F","G","S")
  
  df<-cbind(df,lineage)
  
  # if(TaxlevelTest=="S") {
  #   #ex.seqid.group<-"S"
  #   out<-paste0("tempBLASTDB.",TaxlevelTest,".tsv")
  #   if(file.exists(out)) file.remove(out)
  # }
  
  if(TaxlevelTest=="G") {
    ex.seqid.group<-"S"
    out<-paste0("tempBLASTDB.",TaxlevelTest,".tsv")
    if(file.exists(out)) file.remove(out)
  }
  
  if(TaxlevelTest=="F") {
    ex.seqid.group<-"G"
    out<-paste0("tempBLASTDB.",TaxlevelTest,".tsv")
    file.remove(out)
    if(file.exists(out)) file.remove(out)
  }
  
  #write blastdb
  db<-df[!duplicated(df$qseqid),]
  phylotools::dat2fasta(db[,c("seq.name","seq.text")],"temp.db.fasta")
  #make blastdb using metabinkit
  system2("metabinkit_blastgendb",c("-f","temp.db.fasta","-o", "temp.db","-c"),wait = T)
  a<-Sys.time()
  for(i in 1:length(unique(df$qseqid))){
    
    #message("loop ",i)
    
    #make query fasta
    b<-df[df$qseqid==unique(df$qseqid)[i],]
    b<-b[1,]
    phylotools::dat2fasta(b[,c("seq.name","seq.text")],"temp.seq.fasta")
    
    #get seqids of query group
    ex.seqids<-unique(df[df[,ex.seqid.group]==b[,ex.seqid.group],"qseqid"])  
    write.table(ex.seqids,"ex.seqids.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = F)
    
    #blast (cant use metabin_blast because of -negative_seqidlist)
    # system2(command = "blastdb_aliastool",
    #         args=c("-seqid_file_in", "ex.seqids.txt","-seqid_file_out","ex.seqids.out.txt"), wait = T)
    # 
    system2(command = "blastn",
            args=c("-task", "megablast", "-query", "temp.seq.fasta", "-db","temp.db","-outfmt",
                   "'6 qseqid saccver ssciname evalue staxid pident qcovs qseq'","-evalue",1,"-num_threads", 16, "-max_target_seqs", 
                   100, "-max_hsps",1,"-word_size", 20,"-perc_identity", 70,"-qcov_hsp_perc",98,
                   "-gapopen", 0, "-gapextend", 2, "-reward", 1, "-penalty", -1, "-dust","no", 
                   #  "-negative_seqidlist", "ex.seqids.out.txt", 
                   "-out",
                   "temp.seq.blast.txt"), wait = T)
    #could exclude evalue and qcovs, might save a bit of time
    
    #store blast results
    lblast<-system2("wc",c("-l","temp.seq.blast.txt"),wait=T)
    
    #if(lblast>0) {
    blastResults<-data.table::fread("temp.seq.blast.txt",sep = "\t",data.table = F)
    write.table(blastResults,file = out,append = T,quote = F,row.names = F,sep = "\t",col.names = F)
  }
  
  b<-Sys.time()
  
  #hard coded for now, note the "taxids" rather than "staxid", which metabin does not accept
  headers<-paste0("'1i",paste(c("qseqid", "saccver", "ssciname","evalue", "taxids", "pident", "qcovs","qseq"),collapse = "\t"),"'")
  
  system2("sed",c("-i", headers, out),wait = T)
  
  message("Ouput saved to ",out)
}

plot.thresh<-function(thresher.final.table,limit.plot.to.taxon=NULL,plot.at.level="O"){
  require(ggplot2)
  final.table<-data.table::fread(thresher.final.table,data.table = F)
  
  allcounts<-list()
  
  for(j in 1:length(unique(final.table$level))){
    
    current.level<-unique(final.table$level)[j]
    
    final.tableS<-final.table[final.table$level==current.level,]
    
    nsettings<-length(unique(final.tableS$settings))
    
    countsS<-data.frame(settings=rep("none",nsettings)
                        ,no_hits=rep(0,nsettings),
                        correct=rep(0,nsettings)
                        ,above=rep(0,nsettings),
                        incorrect=rep(0,nsettings),
                        failed=rep(0,nsettings)
    )
    
    countsS$settings<-as.character(countsS$settings)
    
    countsS$level<-current.level
    
    for(i in 1:length(unique(final.tableS$settings))){
      current.setting<-unique(final.tableS$settings)[i]
      countsS$settings[i]<-current.setting
      countsS$no_hits[i]<-sum(final.tableS[final.tableS$settings==current.setting,"no.hits"])
      
      final.tableSx<-final.tableS[final.tableS$no.hits==FALSE,]
      
      countsS$correct[i]<-sum(final.tableSx[final.tableSx$settings==current.setting,paste0("correct",current.level)],na.rm = T)
      countsS$incorrect[i]<-sum(final.tableSx[final.tableSx$settings==current.setting,paste0("incorrect",current.level)],na.rm = T)
      countsS$above[i]<-sum(final.tableSx[final.tableSx$settings==current.setting,paste0("above",current.level)],na.rm = T)
      countsS$failed[i]<-sum(final.tableSx[final.tableSx$settings==current.setting,paste0("failed",current.level)],na.rm = T)
    }
    
    allcounts[[j]]<-countsS
    
  }
  
  allcounts<-do.call(rbind,allcounts)
  
  allcounts$sum<-rowSums(allcounts[,c("no_hits","correct","above","incorrect","failed")])
  
  #write.table(allcounts,paste0(outDir,counts.out),quote = F,sep = "\t",row.names = F)
  
  #######################################################
  #PLOTTING 
  allcounts$settings<-paste0(allcounts$settings,allcounts$level)
  
  longcount<-reshape2::melt(allcounts[,c("no_hits","correct","above","incorrect","failed","settings")],id.vars="settings")
  longcount$level<-substr(longcount$settings,nchar(longcount$settings),nchar(longcount$settings))
  longcount$settings<-substr(longcount$settings,1,nchar(longcount$settings)-1)
  longcount$settings<-gsub("top_","",longcount$settings)
  longcount$settings<-gsub("pidents_","",longcount$settings)
  longcount$settings<-as.factor(longcount$settings)
  
  plot.cols<-c("gray70","yellow4","khaki2","#E31A1C","darkturquoise","green1")
  count.plot.A<-ggplot(data=longcount , aes(y=value, x=settings, fill=variable))+geom_bar(stat = "identity")+
    theme(legend.title = element_text(size=10), legend.text=element_text(size=10),
          axis.text.x=element_text(size=8,angle=45, hjust=1),legend.position="right",legend.direction="vertical") +
    scale_fill_manual(values = plot.cols) + ggtitle("Overview")+
    ylab("Number of queries") + xlab("Settings (top, S, G, F, AF)")+
    facet_wrap(~level,scales = "free_x",)
  
  #print(count.plot)
  
  #############################################  
  #plot by taxon, optional filter by taxa  
  
  out.plot.list<-list()
  
  final.table.list<-split(final.table,final.table$level)
  
  if(!is.null(limit.plot.to.taxon))   message("Limiting plots to ",limit.plot.to.taxon[1],": taxonomic level ",limit.plot.to.taxon[2],
                                              ". Overview plot will not be affected")
  message("Plotting at level ",plot.at.level)
  
  for(i in 1:length(final.table.list)){
    
    xx<-final.table.list[[i]]
    
    levelxx<-xx[1,c("level")]
    longcount<-reshape2::melt(xx[,c("no.hits",paste0("correct",levelxx),paste0("above",levelxx),paste0("incorrect",levelxx),
                                    paste0("failed",levelxx)
                                    ,"settings",paste0("origpath",levelxx))],
                              id.vars=c("settings",paste0("origpath",levelxx)))
    
    colnames(longcount)[2]<-"origpathS"
    
    #extract table to limit by taxon
    if(!is.null(limit.plot.to.taxon)){
      indexTax<-match(limit.plot.to.taxon[2],table = c("K","P","C","O","F","G","S"))
      longcount<-longcount[do.call(rbind,stringr::str_split(longcount$origpathS,";"))[,indexTax]==limit.plot.to.taxon[1],]
    }
    
    #plot at level
    indexn<-match(plot.at.level,table = c("K","P","C","O","F","G","S"))
    if(grep(plot.at.level, c("K","P","C","O","F","G","S"))>grep(levelxx, c("K","P","C","O","F","G","S"))){
      message("Cannot plot at ", plot.at.level," level for ",levelxx,", reverting to ",levelxx)
      indexn<-match(levelxx,table = c("K","P","C","O","F","G","S"))
    }
    longcount$plotpath<-do.call(rbind,stringr::str_split(longcount$origpathS,";"))[,indexn]
    
    #sublabel
    if(levelxx=="S") sublabel<-"Binning outcomes at species level if species is in database"
    if(levelxx=="G") sublabel<-"Binning outcomes at genus level if species is not in database"
    if(levelxx=="F") sublabel<-"Binning outcomes at family level if genus is not in database"
    
    longcount$settings<-gsub("top_","",longcount$settings)
    longcount$settings<-gsub("pidents_","",longcount$settings)
    
    #plot
    plot.cols<-c("gray70","yellow4","khaki2","#E31A1C","darkturquoise","green1")
    count.plot<-ggplot(data=longcount , aes(y=as.numeric(as.logical(value)), x=settings, fill=variable))+geom_bar(stat = "identity")+
      theme(legend.title = element_text(size=10), legend.text=element_text(size=10),
            axis.text.x=element_text(size=8,angle=45, hjust=1),legend.position="right",legend.direction="vertical") +
      scale_fill_manual(values = plot.cols) +
      facet_wrap(~plotpath,scales = "free") + 
      ggtitle(paste0("TestLevel=",levelxx,"; limit=",limit.plot.to.taxon[1],"; PlotLevel=",plot.at.level), subtitle = sublabel) +
      ylab("Number of queries") + xlab("Settings (top, S, G, F, AF)")
    
    
    out.plot.list[[i]]<-count.plot
    names(out.plot.list)[i]<-levelxx
  }
  
  out.plot.list<-list(count.plot.A,out.plot.list)
  
  print(out.plot.list[[1]])
  print(out.plot.list[[2]])
}

bin.thresh<-function(blast.thresh.input,tops=c(0,1,100),
                     pidents.list=list(one=c(99,97,95,90),two=c(98,94,92,88),three=c(93,85,75,60)),
                     known_flags=NULL,final.table.out,SpeciesBL=NULL,GenusBL=NULL,FamilyBL=NULL,ncbiTaxDir= ncbiTaxDir){
  
  sb2<-data.table::fread(blast.thresh.input,data.table = F)
  
  #need to add back origtaxids so we can do taxa disabling
  sb2$origtaxids<-sb2[match(substr(sb2$qseqid,1,nchar(sb2$qseqid)-2),sb2$saccver),"taxids"]
  
  remove.taxa.from.list<-function(BL,sb2,ncbiTaxDir){
    species_bl<-data.table::fread(BL,data.table = F, header = F)
    exclude<-get.children.taxonkit(df = species_bl,column = "V1",ncbiTaxDir)
    sb2<-sb2[!sb2$taxids %in% exclude,]
    sb2<-sb2[!sb2$origtaxids %in% exclude,]
  }
  
  if(!is.null(SpeciesBL)) {
    message("Disabling species")
    sb2<-remove.taxa.from.list(BL = SpeciesBL,sb2 = sb2,ncbiTaxDir)
  }
  
  if(!is.null(GenusBL)) {
    message("Disabling genera")
    sb2<-remove.taxa.from.list(BL = GenusBL,sb2 = sb2,ncbiTaxDir)
  }
  
  if(!is.null(FamilyBL)) {
    message("Disabling families")
    sb2<-remove.taxa.from.list(BL = FamilyBL,sb2 = sb2,ncbiTaxDir)
  }  
  
  write.table(sb2,paste0(blast.thresh.input,"temp"),append = F,sep = "\t",quote = F,row.names = F)
  
  #loop bin at all levels
  binfile.list<-list()
  countloop<-0
  for(j in 1:length(tops)){
    for(k in 1:length(pidents.list)){
      binfile<-paste0("top_",tops[j],".pidents_",paste(pidents.list[[k]],collapse = "."))
      countloop<-countloop+1
      binfile.list[[countloop]]<-binfile
      argsmbk<-c("-i",paste0(blast.thresh.input,"temp"), "-o",binfile,"--no_mbk")
      argsmbk<-c(argsmbk,"-S", pidents.list[[k]][1],"-G", pidents.list[[k]][2],"-F", pidents.list[[k]][3],"-A", pidents.list[[k]][4],
                 "--TopSpecies", tops[j],"--TopGenus",
                 tops[j],"--TopFamily", tops[j]
                 ,"--TopAF", tops[j])
      
      #in the end I dont think this makes sense for threshing, but maybe if removed prior to threshing
      
      # if(!is.null(SpeciesBL)) argsmbk<-c(argsmbk,"--SpeciesBL",SpeciesBL)
      # if(!is.null(GenusBL)) argsmbk<-c(argsmbk,"--GenusBL",GenusBL)
      # if(!is.null(FamilyBL)) argsmbk<-c(argsmbk,"--FamilyBL",FamilyBL)
      
      
      if(!is.null(known_flags)) argsmbk<-c(argsmbk,"--FilterFile",known_flags,"--FilterCol","saccver")
      system2("metabin",argsmbk, wait=T)
    }
  }
  
  file.remove(paste0(blast.thresh.input,"temp"))
  
  outfiles<-do.call(c,binfile.list)
  
  #set metabin args
  
  #read metabin results and merge with sb
  results<-list()
  sb3<-sb2[!duplicated(sb2$qseqid),]
  
  for(i in 1:length(outfiles)){
    a<-data.table::fread(paste0(outfiles[i],".tsv"),data.table = F)
    a$binpath<-paste(a$K,a$P,a$C,a$O,a$F,a$G,a$S,sep = ";")
    sb4<-merge(sb3[,c("qseqid","origpath")],a[,c("qseqid","binpath")],all.x = T,all.y = F,by = "qseqid")
    sb4$settings<-outfiles[i]
    results[[i]]<-sb4
    file.remove(paste0(outfiles[i],".tsv"))
    file.remove(paste0(outfiles[i],".info.tsv"))
    file.remove(paste0(outfiles[i],".versions.txt"))
  }
  
  results<-do.call(rbind, results)
  
  results$level<-substr(results$qseqid,nchar(results$qseqid),nchar(results$qseqid))
  results$binpathS<-results$binpath
  results$binpathG<-path.at.level(results$binpathS,level = "G")
  results$binpathF<-path.at.level(results$binpathS,level = "F")
  results$origpathS<-results$origpath
  results$origpathG<-path.at.level(results$origpathS,level = "G")
  results$origpathF<-path.at.level(results$origpathS,level = "F")
  
  final.table<-results
  
  final.table$no.hits<-FALSE
  
  #correct 
  final.table$correctS<-final.table$binpathS==final.table$origpathS
  final.table$correctG<-final.table$binpathG==final.table$origpathG
  final.table$correctF<-final.table$binpathF==final.table$origpathF
  
  #above desired rank (doesnt check if above rank is correct...)
  final.table$aboveS<-bas.get.ranks(data.frame(taxon=final.table$binpathS))
  final.table$aboveG<-bas.get.ranks(data.frame(taxon=paste0(final.table$binpathG,";NA")))
  final.table$aboveF<-bas.get.ranks(data.frame(taxon=paste0(final.table$binpathF,";NA;NA")))
  
  final.table[(final.table$aboveF=="htf") & final.table$level=="F","aboveF"]<-"above"
  final.table[(final.table$aboveG=="htf" | final.table$aboveG=="family") & final.table$level=="G","aboveG"]<-"above"
  final.table[(final.table$aboveS=="htf" | final.table$aboveS=="family" | final.table$aboveS=="genus") & 
                final.table$level=="S","aboveS"]<-"above"
  
  #incorrect
  final.table$incorrectS<-(!final.table$correctS) & final.table$aboveS=="species"
  final.table$incorrectG<-(!final.table$correctG) & final.table$aboveG=="genus"
  final.table$incorrectF<-(!final.table$correctF) & final.table$aboveF=="family"
  
  #change above from "above" to T/F
  final.table$aboveS<-final.table$aboveS!="species"
  final.table$aboveG<-final.table$aboveG!="genus"
  final.table$aboveF<-final.table$aboveF!="family"
  
  #failed
  final.table$failedS<-final.table$binpathS=="NA;NA;NA;NA;NA;NA;NA"
  final.table$failedG<-final.table$binpathG=="NA;NA;NA;NA;NA;NA"
  final.table$failedF<-final.table$binpathF=="NA;NA;NA;NA;NA"
  
  #if failed, other outcomes are NA
  final.table[final.table$failedS==TRUE,c("correctS","aboveS","incorrectS")]<-NA
  final.table[final.table$failedG==TRUE,c("correctG","aboveG","incorrectG")]<-NA
  final.table[final.table$failedF==TRUE,c("correctF","aboveF","incorrectF")]<-NA
  
  #final.table$settings<-paste0(final.table$Spident,"_",final.table$Gpident,"_",final.table$Fpident,"_top",final.table$top)
  
  #print incorrects
  incorrects<-list()
  incorrects2<-list()
  require(tidyverse)
  for(i in 1:length(unique(final.table$level))){
    levelincor<-unique(final.table$level)[i]
    finalincor<-final.table[final.table$level==levelincor,]
    #remove failed
    finalincor<-finalincor[!is.na(finalincor$correctS),]

    for(j in 1:length(unique(finalincor$settings))){
      finalsets<-unique(finalincor$settings)[j]
      finalincorset<-finalincor[(finalincor$settings==finalsets & finalincor[,paste0("incorrect",levelincor)]==TRUE),]
      if(nrow(finalincorset)>0) incorrects[[j]]<-finalincorset[,c(paste0("origpath",levelincor),paste0("binpath",levelincor),"qseqid","settings")]
      incorrects<-incorrects %>% discard(is.null)
    }

    if(length(incorrects)>0) {
      incorrects2[[i]]<-do.call(rbind,incorrects)
      colnames(incorrects2[[i]])<-c("origpath","binpath","qseqid","settings")
    } else incorrects2[[i]]<-NULL

  }

  incorrects2<-incorrects2 %>% discard(is.null)
  incorrectsdf<-do.call(rbind,incorrects2)
  #incorrectsdf<-incorrectsdf[!duplicated(incorrectsdf$qseqid),]
  message("Results incorrectly binned:")
  print(incorrectsdf)
  
  write.table(final.table,final.table.out,quote = F,sep = "\t",row.names = F)
}

