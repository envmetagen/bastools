message("settings:
        ")
print(ls.str())

source(paste0(bastoolsDir,"master_functions.R"))
source(paste0(bastoolsDir,"bin.blast.R"))
library(processx)
library(dplyr)
library(ggplot2)
library(tidyverse)
obitaxoR<-ROBITaxonomy::read.taxonomy(dbname = obitaxo)
setwd(outDir)
#######################################################
#taxonlimit is to subset the starting ecopcr
if(!is.null(taxonlimit)) if(is.null(taxonlimitlevel)) stop("taxonlimitlevel cannot be NULL if taxonlimit specified")
if(!is.null(taxonlimitlevel)) {
  if(!taxonlimitlevel %in% c("K","P","C","O","F","G","S")) stop("taxonlimitlevel must be one of K,P,C,O,F,G,S")
}

#######################################################
#BUILD ECOPCRDB
if("step1" %in% stepstotake){
  
  message("STEP1 - building ecoPCR database")
  
  obiconvert.Bas(infile = starting_fasta,taxo = obitaxo,in_type = "fasta",out_type = "--ecopcrdb-output",
                 out = gsub(".fasta",".ecopcrdb",starting_fasta))
  
  message("STEP1 COMPLETE")
}

#######################################################
#RUN ECOPCRDB
if("step2" %in% stepstotake){
  
  message("STEP2 - run ecoPCR")

  #run ecopcr with buffer option
  message("Creating references. maxL not used, min_length=10, buffer = F, max_error=",max_error_buildrefs)
  
  ecoPCR.Bas(Pf = Pf,Pr = Pr,ecopcrdb = gsub(".fasta",".ecopcrdb",starting_fasta),max_error = max_error_buildrefs,
           min_length = 10,out = gsub(".fasta",".ecopcrResults.txt",starting_fasta))

  #clean
  clean.ecopcroutput(ecopcrfile = gsub(".fasta",".ecopcrResults.txt",starting_fasta),rm.buffer = F,buffer.used = F)
}

#######################################################
#LOOP BLAST
if("step3" %in% stepstotake){
  
  message("STEP3 - loop blast")

  #read ecopcroutput and add taxonomy
  a<-data.table::fread(gsub(".fasta",".ecopcrResults.txt.clean",starting_fasta),sep = "\t",data.table = F)
  colnames(a)<-gsub("taxid","taxids",colnames(a))
  a<-add.lineage.df(a,ncbiTaxDir)
  a$path<-paste(a$K,a$P,a$C,a$O,a$F,a$G,a$S,sep = ";")
  
  #subset by taxonlimit
  message("ecoPCR table length:",nrow(a))
  if(!is.null(taxonlimit)) {
    a2<-a[a[,taxonlimitlevel]==taxonlimit,]
    message("ecoPCR table length after taxonlimit subset:",nrow(a2))
  } else a2<-a
  
  if(nrow(a2)==0) stop("No sequences left...")
  
  #loop blast for each threshold
  outloopblast.files<-list()
  for(i in 1:length(TaxlevelTest)){
    
    #remove seqs that dont have taxonomy at desired level
    a3<-a2[a2[,TaxlevelTest[i]]!="unknown",]
    
    outloopblast.files[[i]]<-paste0(gsub(".fasta","",starting_fasta),".",primer_set_name,".",TaxlevelTest[i],"_loopBlast.txt")
    
    threshold.blast(ecopcr.clean.df = a3,ncbiTaxDir = ncbiTaxDir,
                  out=outloopblast.files[[i]],TaxlevelTest=TaxlevelTest[i],
                  blast_exec=blast_exec,makeblastdb_exec=makeblastdb_exec,task = task)
  }

  #read blast results, add query paths and plot
  outloopblast.results<-list()
  outloopblast.plots<-list()
  for(i in 1:length(outloopblast.files)){
    outloopblast.results[[i]]<-data.table::fread(outloopblast.files[[i]],sep="\t",data.table = F)
    tophits<-aggregate(V5~V1,outloopblast.results[[i]],max)
    colnames(tophits)<-c("AC","maxpident")
    tophits.m<-merge(tophits,a[,c("AC","path")],by = "AC",all.x = T,all.y = F)
    tophits.m$lev<-path.at.level(tophits.m$path,TaxlevelTest[i])
    outloopblast.plots[[i]]<-ggplot(tophits.m,aes(y=lev,x=maxpident)) + geom_point() + ggtitle(outloopblast.files[[i]])
    ggsave(filename = gsub(".txt",".tophits.pdf",outloopblast.files[[i]]),outloopblast.plots[[i]],device = "pdf", width = 15,height = 10)
  }

message("STEP 3 COMPLETE") 
}

#######################################################
#LOOP BINNING 
if("step4" %in% stepstotake){
  
  message("STEP4 - loop binning")

  results<-list()
  for(i in 1:length(outloopblast.files)){
    #function takes a single file and TaxlevelTest, but loops through tops and pidents
    #so here getting it to loop through TaxLevels, by giving it two blast files (one for each TaxlevelTest)
    results[[i]]<-threshold.bin.blast(blastfile = outloopblast.files[[i]],ecopcr.clean.df = a3,
                                      headers = "qseqid sseqid evalue staxid pident qcovs",
                                      ncbiTaxDir = ncbiTaxDir,max_evalue = 10,min_qcovs = 70,top = tops,
                                      TaxlevelTest=TaxlevelTest[i], pidents=pidents.list[[i]]) 
    
    write.table(results[[i]][[2]],file = gsub(".txt",".bin.results.txt",outloopblast.files[[i]]),append = F,quote = F,sep = "\t",row.names = F)
  }
  
  combocounts<-do.call(rbind,lapply(results, `[[`, 1))
  write.table(combocounts,file = paste0(gsub(".fasta","",starting_fasta),".",primer_set_name,"_summary.counts.txt"),
              append = F,quote = F,sep = "\t",row.names = F)
  
  message("add code to delete leftover files")
  message("STEP 4 COMPLETE") 
  
}
  
#######################################################
#LOOP BINNING 
if("step5" %in% stepstotake){
  
  message("STEP5 - plotting") 
  
  #overview plot
  combocounts<-data.table::fread(paste0(gsub(".fasta","",starting_fasta),".",primer_set_name,"_summary.counts.txt"),data.table = F)
  
  longcount<-data.table::melt(combocounts[,-7],id.vars="file")
  longcount$pident<-rep(combocounts$pident,6)
  longcount$variable <- factor(longcount$variable, levels = c("no.hits","fail.filt","fail.bin","incorrect","above","correct"))
  longcount$pident <- factor(longcount$pident, levels = unique(c(pidents.list[[1]],pidents.list[[2]])))
  
  total.reads<-sum(combocounts[1,-8])
  
  count.plot<-ggplot2::ggplot(data=longcount , aes(y=value, x=pident, fill=variable))+geom_bar(stat = "identity")+
              theme(legend.title = element_text(size=10), legend.text=element_text(size=10),
              axis.text.x=element_text(size=8,angle=45, hjust=1),legend.position="right",legend.direction="vertical")+
              scale_fill_manual(values = MyCols) +
              geom_text(aes(label = paste(round(value/total.reads*100,digits = 0),"%")), 
              position = position_stack(vjust = 0.5), size = 2) +
              facet_wrap(~file,scales = "free") 
  
  ggsave(filename = paste0(gsub(".fasta","",starting_fasta),".",primer_set_name,"_summary.counts.plot.pdf"),
         count.plot,device = "pdf", width = 15,height = 10)
  
#plot by taxon, optional filter by taxa  
  results<-list()
  outloopblast.files<-list()
  for(i in 1:length(TaxlevelTest)){
    outloopblast.files[[i]]<-paste0(gsub(".fasta","",starting_fasta),".",primer_set_name,".",TaxlevelTest[i],"_loopBlast.txt")
    results[[i]]<-data.table::fread(gsub(".txt",".bin.results.txt",outloopblast.files[[i]]),data.table = F)
  }
  
  facet.plots.taxon<-list()
  facet.plots.taxon.list<-list() 
  filelist<-list()
  file.list.list<-list()

  for(j in 1:length(results)) { 
    
    df<-results[[j]]
    
    df$plotpath<-path.at.level(df$origpath,plot.level)
    
    for(i in 1:length(unique(df$file))){
      
      df2<-df[df$file==unique(df$file)[i],]
    
      final.table2<-data.table::melt(df2[,c(length(colnames(df2)),
                                                  grep("rank",colnames(df2)))],id.vars="plotpath")
    
      final.table2$pident<-do.call(rbind,stringr::str_split(final.table2$variable,"_"))[,2]
      final.table2$count<-1 
      colnames(final.table2)<-gsub("value","outcome",colnames(final.table2))
    
      final.table3<-aggregate(count~plotpath+pident+outcome,final.table2,sum)
      final.table3$outcome <- factor(final.table3$outcome, levels = c("no.hits","fail.filt","fail.bin","incorrect","above","correct"))
      final.table3$pident <- factor(final.table3$pident, levels = pidents.list[[j]])
    
      if(!is.null(plot.taxa.limit)) {
        final.table3.list<-list()
        for(i in 1:length(plot.taxa.limit)){
          final.table3.list[[i]]<-final.table3[final.table3$plotpath==plot.taxa.limit[i],]
        }
        final.table3<-do.call(rbind,final.table3.list)
      }
      
      taxa<-do.call(rbind,stringr::str_split(final.table3$plotpath,";"))
      
      final.table3$plotpath<-taxa[,ncol(taxa)]
      
      facet.plots.taxon[[i]]<-ggplot2::ggplot(data=final.table3 , aes(x=pident,y=count,fill=outcome))+geom_bar(stat = "identity")+
                            theme(legend.title = element_text(size=10), legend.text=element_text(size=10),
                            axis.text.x=element_text(size=8,angle=45, hjust=1),legend.position="right",legend.direction="vertical")+
                            scale_fill_manual(values = MyCols) +
                            facet_wrap(~plotpath,scales = "free") +
                            geom_text(aes(label = round(count,digits = 0)), 
                                position = position_stack(vjust = 0.5), size = 2)
      
      filelist[[i]]<-unique(df2$file)
    
      ggsave(filename = paste0(unique(df2$file),".pdf"),plot = facet.plots.taxon[[i]],
             device = "pdf", width = 15,height = 10)
    }
    
    file.list.list[[j]]<-filelist
    facet.plots.taxon.list[[j]]<-facet.plots.taxon
  }
  
  rmarkdown::render(input = paste0(bastoolsDir,"scripts/Plot_threshold_output.Rmd"),output_file = plot.file.html)
  
  message("STEP5 - COMPLETE")
}

