
add.counts.to.biasfile<-function(ncbiTaxDir,download.fasta,after.minL.fasta,after.checks.fasta,first.ecopcr.hit.table,mapped.fasta,out_bias_file){

#catted DLS
cattedDLS<-count.fams.in.fasta(download.fasta,ncbiTaxDir)

#after minL
afterminL<-count.fams.in.fasta(after.minL.fasta,ncbiTaxDir)

#after rm fams and checks
afterchecks<-count.fams.in.fasta(after.checks.fasta,ncbiTaxDir)

#after first ecopcr
firstecopcr<-data.table::fread(first.ecopcr.hit.table,sep = "\t")
firstecopcr$taxids<-firstecopcr$taxid
firstecopcr<-add.lineage.df(firstecopcr,ncbiTaxDir)
firstecopcr$count<-1
firstecopcr$path<-paste(firstecopcr$K,firstecopcr$P,firstecopcr$C,firstecopcr$O,firstecopcr$F,sep = ";")
a<-aggregate(firstecopcr$count,by=list(firstecopcr$path),FUN=sum)
colnames(a)<-c("Family","nseqs")
firstecopcr<-a

#after mapping back
aftermapping<-count.fams.in.fasta(mapped.fasta,ncbiTaxDir)

#merge all
merged<-merge(cattedDLS,afterminL,by = "Family",all = T)
colnames(merged)<-gsub("nseqs.x","downloaded",colnames(merged))
colnames(merged)<-gsub("nseqs.y","after_min_length",colnames(merged))
merged<-merge(merged,afterchecks,by = "Family",all = T)
colnames(merged)<-gsub("nseqs","after_checks",colnames(merged))
merged<-merge(merged,firstecopcr,by = "Family",all = T)
colnames(merged)<-gsub("nseqs","first_ecopcr",colnames(merged))
merged<-merge(merged,aftermapping,by = "Family",all = T)
colnames(merged)<-gsub("nseqs","after_mapping_back",colnames(merged))

sum(merged$downloaded,na.rm = T)
sum(merged$after_min_length,na.rm = T)
sum(merged$after_checks,na.rm = T)
sum(merged$firstecopcr,na.rm = T)
sum(merged$after_mapping,na.rm = T)

biastemp<-data.table::fread(out_bias_file,sep = "\t")
biastemp$path<-paste(biastemp$K,biastemp$P,biastemp$C,biastemp$O,biastemp$in.odb,sep = ";")

mergedbias<-merge(merged,biastemp,by.x ="Family",by.y = "path",all = T)
mergedbias$in.odb=NULL
mergedbias$nseqs.odb=NULL

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

