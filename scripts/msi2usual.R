##msi.2.single.fasta.and.single.otutab

#transform the ouptput fastas from msi (fastas in multiple folders) into usual for illscript3,
#i.e. the same as the output from mbc

##start
require(tidyverse)
setwd(root_folder)
dmplxd_dirs<-grep("barcode.*",list.dirs(recursive = F),value = T)
master<-data.table::fread(master_sheet,data.table = F)
if(grep(".fasta$",fasta.out)==0) stop("fasta.out must end in '.fasta'")
if(length(grep(".tsv$",otutab.out))==0) stop("otutab.out must end in '.tsv'")

##read fastas
  fasta_files<-paste0(dmplxd_dirs,gsub(pattern = "./","/",dmplxd_dirs),".centroids.fasta")
  fasta.list<-list()
  for(i in 1:length(fasta_files)){
    fasta.list[[i]]<-phylotools::read.fasta(fasta_files[i])
  }
  
  all.fasta<-do.call(rbind,fasta.list)
  message(nrow(all.fasta), " sequences,", length(unique(all.fasta$seq.text))," unique")
  
##metadata
  all.fasta$metadata<-stringr::str_split(all.fasta$seq.name,"size=",simplify = T)[,3]
  all.fasta$nreads<-as.numeric(stringr::str_split(all.fasta$metadata,":",simplify=T)[,1])
  all.fasta$primer_found<-stringr::str_split(all.fasta$metadata,":",simplify=T)[,2]
  all.fasta$barcode<-stringr::str_split(all.fasta$seq.name,"barcode=",simplify = T)[,2]
  all.fasta$barcode<-stringr::str_split(all.fasta$barcode,":",simplify = T)[,1]
  all.fasta$primer_used<-stringr::str_split(all.fasta$metadata,"adapter=",simplify=T)[,2]
  all.fasta$primer_used<-stringr::str_split(all.fasta$primer_used,":",simplify=T)[,1]
  
  all.fasta$otu<-paste0("OTU.",1:nrow(all.fasta))
  all.fasta$seq.name<-paste0(all.fasta$otu,";size=",all.fasta$nreads)
  
  total_reads<-sum(all.fasta$nreads)
  message("Total reads=",total_reads)
  
##rm no_adapter
  before<-nrow(all.fasta)
  all.fasta<-all.fasta[all.fasta$primer_found!="adapter=no_adapter",]
  message(before-nrow(all.fasta)," otus (",total_reads-sum(all.fasta$nreads)," reads) removed, no adapters")
  
##split by primer
  all.fasta.spl<-split(all.fasta,f=all.fasta$primer_used)

##otu tabs
  otutab.list<-list()
  for(i in 1:length(all.fasta.spl)){
    otutab.list[[i]]<-tidyr::pivot_wider(names_from = "barcode",values_from="nreads",
                             data=all.fasta.spl[[i]][,c("otu","nreads","barcode")])
    otutab.list[[i]][is.na(otutab.list[[i]])]<-0
    
    colnames(otutab.list[[i]])<-c("#OTU ID",
    master[match(colnames(otutab.list[[i]][,-1]),master$barcode_name),"ss_sample_id"])
  }
  
##write files
  for(i in 1:length(all.fasta.spl)){
    phylotools::dat2fasta(all.fasta.spl[[i]][,c("seq.name","seq.text")],
                          outfile = gsub(".fasta$",paste0("_",all.fasta.spl[[i]]$primer_used[1],".fasta"),fasta.out))
    write.table(otutab.list[[i]],quote = F,sep = "\t",row.names = F,
                file=gsub(".tsv$",paste0("_",all.fasta.spl[[i]]$primer_used[1],".tsv"),otutab.out))
  }
  