library(processx)
library(dplyr)
library(mgsub)

setwd("/home/bastian.egeter/git_bastools/bastools/")
file.sources<-c("add.taxids.fasta.BAS.R","bin.blast.R","merge_MBC_otutab_with_bin_blast.R","blast.min.bas.R",
                "taxatab.filter.R","do.count.R","googlesheet.foos.R")
sapply(file.sources,source)
obitaxdb<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/obitaxdump/obitaxdb"
ncbiTaxDir<-"/mnt/Disk1/Tools/BLAST+/DBs/nt_taxonomy/taxdump/October-2019/"
setwd("/mnt/Disk1/Minion_data/2019_August_002_Mussels/test_incr1/")
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


############################################################################
#make otu table





