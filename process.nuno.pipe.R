library(processx)
library(dplyr)
library(mgsub)

setwd("/home/bastian.egeter/git_bastools/bastools/")
file.sources<-c("add.taxids.fasta.BAS.R","bin.blast.R","merge_MBC_otutab_with_bin_blast.R","blast.min.bas.R",
                "taxatab.filter.R","do.count.R","googlesheet.foos.R")
sapply(file.sources,source)
googlesheets::gs_auth()
obitaxdb<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/obitaxdump/obitaxdb"
ncbiTaxDir<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/taxdump/October-2019/"
#setwd("/mnt/Disk1/Minion_data/2019_August_002_Mussels/test_incr1/")
setwd("/mnt/Disk1/Minion_data/2019_September_001/res/")
filterpc<-0.1
#out.taxatable<-"2019_August_002_Mussels_ALL_PRIMERS_nuno.then.bas.taxatable.tf.txt"
out.taxatable<-"2019_September_001_ALL_PRIMERS_nuno.then.bas.taxatable.tf.txt"


#############################################################################
#FILTER BLAST
nuno.res<-data.table::fread("results.tsv",header = T,sep = "\t",fill = T,data.table = F)
filter.blast.nuno(blastfile = nuno.res,ncbiTaxDir = ncbiTaxDir,
                  out = "results.blast.filt.txt")

#############################################################################
#BIN READS
files<-list.files(pattern = ".blast.filt.txt")
for(i in 1:length(files)){
  message(paste("binning filtered blast results for",files[i]))
  filtered_blastfile<-files[i]
  binfile<-gsub(".blast.filt.txt",".bins.txt",files[i])
  bin.blast2(filtered_blastfile = filtered_blastfile,ncbiTaxDir = ncbiTaxDir,
             obitaxdb = obitaxdb,out = binfile,spident = 95,gpident = 92,fpident = 80,abspident = 70)
}

############################################################################
#make taxon table

#deal with multiple hits - 
  #1st order by pident
  nuno.res2<-nuno.res[order(-nuno.res$pident),]
  #then rm dup reads
  nuno.res2<-nuno.res2[!duplicated(nuno.res2$read),]

#next extract counts
  nuno.res2$count<-as.numeric(do.call(rbind,stringr::str_split(nuno.res2$read,":size="))[,3])

#merge with bins
  taxon_input<-data.table::fread(file = binfile, sep = "\t")
  taxon_input$path<-paste0(taxon_input$K,";",taxon_input$P,";",taxon_input$C,";",taxon_input$O,";",
                             taxon_input$F,";",taxon_input$G,";",taxon_input$S)
  merged.table<-merge(taxon_input[,c("qseqid","path")],nuno.res2[,c("read","count")],
                        by.x = "qseqid",by.y = "read",all = TRUE)
  merged.table<-merged.table[!is.na(merged.table$count),]
  
  
  barcodes<-do.call(rbind,stringr::str_split(merged.table$qseqid,"barcode="))[,2]
  merged.table$barcode<-do.call(rbind,stringr::str_split(barcodes,":"))[,1]
    
  taxatable<-reshape2::dcast(merged.table[,c("path","barcode","count")],path~barcode,value.var = "count",
                               fun.aggregate = sum)
  
  taxatable$path[is.na(taxatable$path)]<-"No_hits"
  
  #Apply taxa filter
  rownames(taxatable)<-taxatable$path
  taxatable$path=NULL
  taxatable<-taxatable[rowSums(taxatable)!=0,]
  taxatab.PCS<-sweep(taxatable, MARGIN = 1, STATS = rowSums(taxatable), FUN = "/")*100
  taxatab.PCS[taxatab.PCS<filterpc]<-0 
  taxatable[taxatab.PCS==0]<-0
  taxatable<-taxatable[rowSums(taxatable)!=0,]
  
  taxatable$taxon<-rownames(taxatable)
  taxatable<-taxatable[,c(length(colnames(taxatable)),1:(length(colnames(taxatable))-1))]
  write.table(taxatable,out.taxatable,row.names = F,quote = F,sep = "\t")
  
  
filter.blast.nuno<-function(blastfile,ncbiTaxDir,out, max_evalue=0.001,top=1){
  
  message("reading blast results")
  btab<-as.data.frame(blastfile)
  
  message("applying global max_evalue threshold")
  btab<-btab[btab$evalue<max_evalue,]
  message("applying global top threshold")
  if(length(btab[,1])==0) stop("No hits passing min_qcovs and max_evalue thresholds")
  topdf<-aggregate(x = btab[,colnames(btab) %in% c("read","pident")],by=list(btab$read),FUN = max)
  topdf$min_pident<-topdf$pident-topdf$pident*top/100
  btab<-merge(btab,topdf[,c(2,4)],by = "read", all = T) #can definitely do this differently and faster
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
  
  colnames(btab)<-gsub("read","qseqid",colnames(btab))
  
  write.table(x = btab,file = out,sep="\t",quote = F,row.names = F)
  
}








