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
obitab_bin_blast_merge_minion<-function(obitabfile,binfile,mastersheetfile=NA,experiment_id,used.obiuniq=F,out){
  
  if(used.obiuniq==F){
    obitab_input<-data.table::fread(obitabfile, sep = "\t")
    taxon_input<-data.table::fread(file = binfile, sep = "\t")
    #obitab_input$barcode<-gsub("barcode=","",do.call(rbind,stringr::str_split(obitab_input$definition," "))[,3])
    taxon_input$path<-paste0(taxon_input$K,";",taxon_input$P,";",taxon_input$C,";",taxon_input$O,";",
                             taxon_input$F,";",taxon_input$G,";",taxon_input$S)
    
    merged.table<-merge(taxon_input[,c("qseqid","path")],obitab_input[,c("id","barcode")],
                        by.x = "qseqid",by.y = "id",all = TRUE)
    merged.table$count<-1
    
    taxatable<-reshape2::dcast(merged.table[,c("path","barcode","count")],path~barcode,value.var = "count",
                    fun.aggregate = sum)
    taxatable$path[is.na(taxatable$path)]<-"No_hits"
    
    if(!is.null(mastersheetfile)){
    #replace barcodes with sample_names
    final.barcodes<-as.data.frame(colnames(taxatable[,2:length(colnames(taxatable))]))
    colnames(final.barcodes)<-"barcode_id"
    master<-data.table::fread(file = mastersheetfile, sep = "\t",data.table = F)
    mapping.samples<-master[grep(experiment_id, master$experiment_id),c("barcode_id","ss_sample_id")]
    final.samples<-merge(final.barcodes,mapping.samples,by = "barcode_id",all.y = F,all.x = T)
    colnames(taxatable)<-c("taxon",final.samples$ss_sample_id)
    }
    
    write.table(x = taxatable,file = out,sep="\t",quote = F,row.names = F)
  }
  
  if(used.obiuniq==T) message("Havent writtent his yet!")
}

######################################################################
#SPLIT TAXATABLES
splice.taxatables<-function(files,mastersheet){
#read files
taxatables<-list()
for(i in 1:length(files)){
  taxatables[[i]]<-data.table::fread(files[i],data.table = F,sep = "\t")
}

#split by project
ssdf<-data.table::fread(file = mastersheet,sep = "\t")
projectnames<-unique(do.call(rbind,stringr::str_split(string = ssdf$ss_sample_id,pattern = "-"))[,1])
projectnames<-paste0(projectnames,"-")

taxatablesplit<-list()
taxatablesplit2<-list()

for(i in 1:length(taxatables)){
  for(j in 1:length(projectnames)){
    taxatablesplit[[j]]<-taxatables[[i]][,grep(projectnames[j],colnames(taxatables[[i]]))]
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





