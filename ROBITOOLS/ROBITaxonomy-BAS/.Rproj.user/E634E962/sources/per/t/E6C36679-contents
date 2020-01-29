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
  merged.table<-merge(taxon_input[,c("qseqid","path")],otutab_input[,c("id","barcode","count")],
                      by.x = "qseqid",by.y = "id",all = TRUE)
    
    taxatable<-reshape2::dcast(merged.table[,c("path","barcode","count")],path~barcode,value.var = "count",
                    fun.aggregate = sum)
    taxatable$path[is.na(taxatable$path)]<-"No_hits"
    
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

#merging taxatables 
bas.merge.taxatables<-function(taxatabs){
  
  taxatabs.list<-list()
  counts<-data.frame(file=taxatabs,reads=0)
  for(i in 1:length(taxatabs)){
    taxatabs.list[[i]]<-data.table::fread(taxatabs[i],sep = "\t",data.table = F)
    counts[i,2]<-sum(taxatabs.list[[i]][,2:length(colnames(taxatabs.list[[i]]))])
    message(counts[i,1])
    message(counts[i,2])
  }
  
  library(tidyverse)
  all.taxatabs<-taxatabs.list %>% reduce(full_join, by = "taxon")
  
  #remove NAs
  all.taxatabs[is.na(all.taxatabs)]<-0
  
  #check sums
  a<-sum(all.taxatabs[,2:length(colnames(all.taxatabs))])
  message(paste("merged taxatable",a))
  if(a==sum(counts$reads)) message("Read counts all good") else stop("Read counts do not match")
  
  return(all.taxatabs)
  
}





