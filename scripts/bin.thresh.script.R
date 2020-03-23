message("settings:
        ")
print(ls.str())

source(paste0(bastoolsDir,"master_functions.R"))
library(processx)
library(dplyr)
library(ggplot2)
library(tidyverse)
setwd(outDir)
#######################################################
#check taxonlimit ok 
if(!is.null(taxonlimit)) if(is.null(taxonlimitlevel)) stop("taxonlimitlevel cannot be NULL if taxonlimit specified")

#set some file names
outloopblast.files<-list()
ecopcr.files<-list()
for(i in 1:length(TaxlevelTest)){
  outloopblast.files[[i]]<-paste0(gsub(".fasta","",starting_fasta),".",primer_set_name,".",TaxlevelTest[i],"_loopBlast.txt")
  ecopcr.files[[i]]<-paste0(gsub(".fasta","",starting_fasta),".",primer_set_name,".",TaxlevelTest[i],"_ecopcrDB.txt")
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
           min_length = 10,out = gsub(".fasta",paste0(".",primer_set_name,".ecopcrResults.txt"),starting_fasta))

  #clean
  clean.ecopcroutput(ecopcrfile = gsub(".fasta",paste0(".",primer_set_name,".ecopcrResults.txt"),starting_fasta),rm.buffer = F,buffer.used = F)
  
  #read ecopcroutput and add taxonomy
  a<-data.table::fread(paste0(gsub(".fasta",paste0(".",primer_set_name,".ecopcrResults.txt"),starting_fasta),".clean"),sep = "\t",data.table = F)
  colnames(a)<-gsub("taxid","taxids",colnames(a))
  a<-add.lineage.df(a,ncbiTaxDir)
  a$path<-paste(a$K,a$P,a$C,a$O,a$F,a$G,a$S,sep = ";")
  write.table(a,paste0(gsub(".fasta",paste0(".",primer_set_name,".ecopcrResults.txt"),starting_fasta),".clean"),append = F,quote = F,row.names = F,sep = "\t")
  
  #krona plot
  a$count<-1
  a<-a[,c("path","count")]
  write.table(a,paste0(gsub(".fasta",paste0(".",primer_set_name,".ecopcrResults.txt"),starting_fasta),".clean.Forkrona.txt")
              ,append = F,quote = F,row.names = F,sep = "\t")
              
  bas.krona.plot(taxatable = paste0(gsub(".fasta",paste0(".",primer_set_name,".ecopcrResults.txt"),starting_fasta),".clean.Forkrona.txt"))
}


#######################################################
#LOOP BLAST
if("step3" %in% stepstotake){
  
  message("STEP3 - loop blast")
  a<-data.table::fread(paste0(gsub(".fasta",paste0(".",primer_set_name,".ecopcrResults.txt"),starting_fasta),".clean"),sep = "\t",data.table = F)
  
  #subset by taxonlimit
  message("ecoPCR table length:",nrow(a))
  a.list<-list()
  if(!is.null(taxonlimit)) {
    for(i in 1:length(taxonlimit)){
      a.list[[i]]<-a[a[,taxonlimitlevel[i]]==taxonlimit[i],]
    }
    a2<-do.call(rbind,a.list)
    message("ecoPCR table length after taxonlimit subset:",nrow(a2))
  } else a2<-a
  
  if(nrow(a2)==0) stop("No sequences left...")
  
  #if(assume.no.other.taxa==T) a2<-a2[a2$path %in% plot.taxa.limit,]
  
  #loop blast for each taxonlevel
  
  for(i in 1:length(TaxlevelTest)){
    
    #remove seqs that dont have taxonomy at desired level
    a3<-a2[a2[,TaxlevelTest[i]]!="unknown",]
    
    #if TaxlevelTest[i]="S" remove sp.-type entries
    if(TaxlevelTest[[i]]=="S"){
      message("Removing species with 'sp.', numbers or more than one space")
      count1<-nrow(a3)
      if(length(grep(" sp\\.",a3$S,ignore.case = T))>0) a3<-a3[-grep(" sp\\.",a3$S,ignore.case = T),]
      if(length(grep(" .* .*",a3$S,ignore.case = T))>0) a3<-a3[-grep(" .* .*",a3$S,ignore.case = T),]
      if(length(grep("[0-9]",a3$S))>0) a3<-a3[-grep("[0-9]",a3$S),]
      message(count1-nrow(a3), " sequences removed")
    }
    
    #save ecopcroutput used for later
    write.table(a3,ecopcr.files[[i]],append = F,quote = F,sep = "\t",row.names = F)
    
    threshold.blast(ecopcr.clean.df = a3,ncbiTaxDir = ncbiTaxDir,
                  out=outloopblast.files[[i]],TaxlevelTest=TaxlevelTest[i],
                  blast_exec=blast_exec,makeblastdb_exec=makeblastdb_exec,task = task,all.sp.in.db = all.sp.in.db)
    }
  

  #read blast results, add query paths and plot
  outloopblast.results<-list()
  outloopblast.plots<-list()
  for(i in 1:length(outloopblast.files)){
    outloopblast.results[[i]]<-data.table::fread(outloopblast.files[[i]],sep="\t",data.table = F)
    tophits<-aggregate(V5~V1,outloopblast.results[[i]],max)
    colnames(tophits)<-c("AC","pident")
    tophits$minmax<-"max"
    minhits<-aggregate(V5~V1,outloopblast.results[[i]],min)
    colnames(minhits)<-c("AC","pident")
    minhits$minmax<-"min"
    tophits<-rbind(tophits,minhits)
    tophits.m<-merge(tophits,a[,c("AC","path")],by = "AC",all.x = T,all.y = F)
    if(TaxlevelTest[i]!="S") {
      tophits.m$lev<-path.at.level(tophits.m$path,TaxlevelTest[i])
    } else {
      tophits.m$lev<-tophits.m$path
    }
    outloopblast.plots[[i]]<-ggplot(tophits.m,aes(y=lev,x=pident,colour=minmax)) + geom_point() + ggtitle(outloopblast.files[[i]])
    ggsave(filename = gsub(".txt",".tophits.pdf",outloopblast.files[[i]]),outloopblast.plots[[i]],device = "pdf", width = 15,height = 10)
  }

message("STEP 3 COMPLETE") 
}

#######################################################
#LOOP BINNING 
if("step4" %in% stepstotake){
  
  message("STEP4 - loop binning")
  
  results<-list()
  for(i in 1:length(TaxlevelTest)){
    #function takes a single file and TaxlevelTest, but loops through tops and pidents
    #so here getting it to loop through TaxLevels, by giving it one blast file for each TaxlevelTest
    
    ecopcr<-data.table::fread(ecopcr.files[[i]],data.table = F)
    
    results[[i]]<-threshold.bin.blast(blastfile = outloopblast.files[[i]],ecopcr.clean.df = ecopcr,
                                      headers = "qseqid sseqid evalue staxid pident qcovs",
                                      ncbiTaxDir = ncbiTaxDir,max_evalue = 20,min_qcovs = 70,top = tops,
                                      TaxlevelTest=TaxlevelTest[i], pidents=pidents.list[[i]]) 
    
    write.table(results[[i]][[2]],file = gsub(".txt",".bin.results.txt",outloopblast.files[[i]]),append = F,quote = F,sep = "\t",row.names = F)
  }
  
  #finally, do a check if species not in db, do we get to family level
  if(all.sp.in.db==F) if("G" %in% TaxlevelTest) {
    results[[i+1]]<-threshold.bin.blast(blastfile = grep(".G_",outloopblast.files,value = T),ecopcr.clean.df = ecopcr,
                                       headers = "qseqid sseqid evalue staxid pident qcovs",
                                       ncbiTaxDir = ncbiTaxDir,max_evalue = 20,min_qcovs = 70,top = tops,
                                       TaxlevelTest="F", pidents=pidents.list[[grep("F",TaxlevelTest)]])
    
    results[[i+1]][[1]]$file<-gsub("loopBlast.txt_","loopBlast.txt_SbyF_",results[[i+1]][[1]]$file)
    results[[i+1]][[2]]$file<-gsub("loopBlast.txt_","loopBlast.txt_SbyF_",results[[i+1]][[2]]$file)
    
   write.table(results[[i+1]][[2]],file = gsub(".txt",".SbyF.bin.results.txt",grep(".G_",outloopblast.files,value = T)),
              append = F,quote = F,sep = "\t",row.names = F)
  }
  
  
  combocounts<-do.call(rbind,lapply(results, `[[`, 1))
  write.table(combocounts,file = paste0(gsub(".fasta","",starting_fasta),".",primer_set_name,"_summary.counts.txt"),
              append = F,quote = F,sep = "\t",row.names = F)
  
  message("add code to delete leftover files")
  message("STEP 4 COMPLETE") 
  
}
  
#######################################################
#PLOTTING 
if("step5" %in% stepstotake){
  
  message("STEP5 - plotting") 
  
  #overview plot
  combocounts<-data.table::fread(paste0(gsub(".fasta","",starting_fasta),".",primer_set_name,"_summary.counts.txt"),data.table = F)
  
  longcount<-reshape2::melt(combocounts[,-7],id.vars="file")
  longcount$pident<-rep(combocounts$pident,6)
  longcount$variable <- factor(longcount$variable, levels = c("no.hits","fail.filt","fail.bin","incorrect","above","correct"))
  longcount$pident <- factor(longcount$pident, levels = sort(unique(unlist(pidents.list)),decreasing = T))
  longcount$settings<-do.call(rbind,stringr::str_split(longcount$file,"loopBlast.txt_"))[,2]
  longcount$taxlevel<-do.call(rbind,stringr::str_split(longcount$file,"loopBlast.txt_"))[,1]
  longcount$taxlevel<-stringr::str_sub(longcount$taxlevel,start = -2)
  longcount$settings<-paste0(longcount$taxlevel,longcount$settings)
  
  b<-list()
  b2<-list()
  for(j in 1:length(TaxlevelTest)){
    for(i in 1:length(tops)){
      greps<-paste0("^",TaxlevelTest[j],"(.)*top_",sort(tops)[i],"$")
      b[[i]]<-grep(greps,unique(longcount$settings))
    }
    b2[[j]]<-b
  }
  
  longcount$settings<-factor(longcount$settings, levels = unique(longcount$settings))
  
  
  plot.cols<-c("gray70","yellow4","khaki2","#E31A1C","darkturquoise","green1")
  names(plot.cols)<-c("no.hits","fail.filt","fail.bin","incorrect","above","correct")
  
  total.reads<-sum(combocounts[1,-8])
  
  count.plot<-ggplot2::ggplot(data=longcount , aes(y=value, x=pident, fill=variable))+geom_bar(stat = "identity")+
              theme(legend.title = element_text(size=10), legend.text=element_text(size=10),
              axis.text.x=element_text(size=8,angle=45, hjust=1),legend.position="right",legend.direction="vertical")+
              scale_fill_manual(values = plot.cols) +
              geom_text(aes(label = paste(round(value/total.reads*100,digits = 0),"%")), 
              position = position_stack(vjust = 0.5), size = 2) +
              facet_wrap(~settings,scales = "free") 
  
  ggsave(filename = paste0(gsub(".fasta","",starting_fasta),".",primer_set_name,"_summary.counts.plot.pdf"),
         count.plot,device = "pdf", width = 15,height = 10)
  
  
#############################################  
#plot by taxon, optional filter by taxa  
  results<-list()
  for(i in 1:length(outloopblast.files)){
    results[[i]]<-data.table::fread(gsub(".txt",".bin.results.txt",outloopblast.files[[i]]),data.table = F)
  }
  if(all.sp.in.db==F) if("F" %in% TaxlevelTest) results[[i+1]]<-data.table::fread(gsub("F_loopBlast.bin.results.txt","G_loopBlast.SbyF.bin.results.txt",
                                                                   gsub(".txt",".bin.results.txt",grep(".F_",outloopblast.files,value = T)))
                                                                   ,data.table = F)
  
  facet.plots.taxon<-list()
  facet.plots.taxon.list<-list() 
  filelist<-list()
  file.list.list<-list()

  for(j in 1:length(results)) { 
    
    df<-results[[j]]
    
    for(l in 1:length(plot.levels)){
      
      df$plotpath<-path.at.level(pathvector = df$origpath,level = plot.levels[l])
      
      if(!is.null(plot.taxa.limit)) {
        plot.taxa.limit2<-path.at.level(pathvector = plot.taxa.limit,level = plot.levels[l]) 
        }else plot.taxa.limit2<-NULL
      
      for(i in 1:length(unique(df$file))){
        
        df2<-df[df$file==unique(df$file)[i],]
        
        if(!is.null(plot.taxa.limit2)) {
          df2.list<-list()
          for(n in 1:length(plot.taxa.limit2)){
            df2.list[[n]]<-df2[df2$plotpath==plot.taxa.limit2[n],]
          }
          df2<-do.call(rbind,df2.list)
        
        
          #counts of taxa in final plots
          #report some counts
          taxa.counts.list<-list()
          taxlevels<-c("K","P","C","O","F","G")
          plot.title.list<-list()
          for(k in 1:length(plot.taxa.limit2)){
            for(m in 1:length(taxlevels)) {
              taxa.counts<-path.at.level(df2.list[[k]]$path,taxlevels[m])
              taxa.counts.list[[m]]<-paste0(taxlevels[m],":",length(unique(taxa.counts)))
            }
            plot.title.list[[k]]<-paste(plot.taxa.limit2[k],toString(unlist(taxa.counts.list)),",S:",length(unique(df2.list[[k]]$path)))
          }
          plot.title<-toString(unlist(plot.title.list))
        }
      
        final.table2<-reshape2::melt(df2[,c(length(colnames(df2)),
                                                    grep("rank",colnames(df2)))],id.vars="plotpath")
      
        final.table2$pident<-do.call(rbind,stringr::str_split(final.table2$variable,"_"))[,2]
        final.table2$count<-1 
        colnames(final.table2)<-gsub("value","outcome",colnames(final.table2))
      
        final.table3<-aggregate(count~plotpath+pident+outcome,final.table2,sum)
        final.table3$outcome <- factor(final.table3$outcome, levels = c("no.hits","fail.filt","fail.bin","incorrect","above","correct"))
        
        if(length(grep("SbyF", unique(df$file)[i]))>0) pident.levs<-pidents.list[[grep("F",TaxlevelTest)]] else pident.levs<-pidents.list[[j]]

        final.table3$pident <- factor(final.table3$pident, levels = pident.levs)
      
        taxa<-do.call(rbind,stringr::str_split(final.table3$plotpath,";"))
        
        final.table3$plotpath<-taxa[,ncol(taxa)]
        
        facet.plots.taxon[[i]]<-ggplot2::ggplot(data=final.table3 , aes(x=pident,y=count,fill=outcome))+geom_bar(stat = "identity")+
                              theme(legend.title = element_text(size=10), legend.text=element_text(size=10),
                              axis.text.x=element_text(size=8,angle=45, hjust=1),legend.position="right",legend.direction="vertical")+
                              scale_fill_manual(values = plot.cols) +
                              facet_wrap(~plotpath,scales = "free") +
                              geom_text(aes(label = round(count,digits = 0)), 
                                  position = position_stack(vjust = 0.5), size = 2)
        
        if(!is.null(plot.taxa.limit)) {
          facet.plots.taxon[[i]]<-facet.plots.taxon[[i]]+labs(caption=toString(plot.title))+
            theme(plot.caption=element_text(size=8, hjust=0, margin=margin(15,0,0,0)))
        }
        
        filelist[[i]]<-paste0(unique(df2$file),"_plotlev",plot.levels[l],".pdf")
      
        ggsave(filename = paste0(unique(df2$file),"_plotlev",plot.levels[l],".pdf"),plot = facet.plots.taxon[[i]],
               device = "pdf", width = 15,height = 10)
      }
    }
      
    file.list.list[[j]]<-filelist
    facet.plots.taxon.list[[j]]<-facet.plots.taxon
  }
  
  #rmarkdown::render(input = paste0(bastoolsDir,"scripts/Plot_threshold_output.Rmd"),output_file = plot.file.html)
  #should probably fix this
  message("STEP5 - COMPLETE")
}

