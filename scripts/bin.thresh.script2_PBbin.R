#steps
#1. do thorough blast with sseq in results (must have sseq, staxid, saccver)
#2. make blastdb from results (metabinkit)
#3. do thresh blast (blastn)
#4. do thresh binning (metabinkit)

#blast reads e.g.
###settings used for this example
# '/home/bastian.egeter/Tools/ncbi-blast-2.9.0+/bin/blastn' -query \
# 2019_August_002.VENE.lenFilt.trimmed.ids.SC4.pol.fasta -db \
# '/mnt/Disk1/Tools/BLAST+/DBs/nt_v5/nt' -task blastn -outfmt \
# '6 qseqid pident qcovs saccver staxid ssciname sseq' -num_threads 64 \
# -max_target_seqs 500 -max_hsps 1 -word_size 11 -perc_identity 50 \
# -qcov_hsp_perc 98 -gapopen 0 -gapextend 2 -reward 1 -penalty -1 -dust no \
# -out 2019_August_002.VENE.lenFilt.trimmed.ids.SC4.pol.blast.txt


outDir<-"/home/bastian.egeter/PB_tests/"
input<-"/home/bastian.egeter/git_bastools/bastools/beta/PB_classifier/2019_August_002.VENE.lenFilt.trimmed.ids.SC4.pol.blast.txt" #example
bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/" #needed only for add.lineage.df
ncbiTaxDir<-"/home/bastian.egeter/metabinkit.install/db/"
=======
##Note: was getting quite a few "incorrects" with max_target_seqs 50 and not choosing BEST staxid 
#(was previously choosing one, potentially random, ssacver)

#TODO do a check if species not in db, do we get to family level
#TODO add code to delete leftover files
#TODO divide into chunks: 1. loop blast 2. loop bin & results 3. plotting 
#TODO dont run top*pident combos that dont make sense
#TODO count "failed binning thresholds" rather than lumping it into "above"

outDir<-"/home/tutorial/SCRIPTS/PB_tests/"
input<-"/home/tutorial/TOOLS/bastools/beta/PB_classifier/2019_August_002.VENE.lenFilt.trimmed.ids.SC4.pol.blast.txt" #example
bastoolsDir<-"/home/tutorial/TOOLS/bastools/" #needed only for add.lineage.df
ncbiTaxDir<-"/home/tutorial/TOOLS/metabinkit.install/db/"
TaxlevelTest=c("S","G","F")

#binning settings to loop through
tops<-c(0.1,1,10,100)
pidents.list<-list(Spident=c(100,98,95,93),Gpident=c(99,97,95,93),Fpident=c(99,90,80,70)) 

#plotting
plot.at.level<-"O"
limit.plot.to.taxon<-c("Bivalvia","C")

#saving main results (can be fed to plotting chunk)
counts.out<-"PBthresh_counts.tsv"
final.table.out<-"PBthresh_final.table.tsv"

create.extended.db=T
nt="/mnt/Disk1/Tools/BLAST+/DBs/nt_v5/nt"

source(paste0(bastoolsDir,"master_functions.R"))
library(ggplot2)
library(tidyverse)


setwd(outDir)
#######################################################
#makeblastdb

#read full blast results
fullblast<-data.table::fread(input,data.table = F)

#remove hyphens sseq
fullblast$sseq<-gsub("-","",fullblast$sseq)

#need to reduce db size. decided best seq per taxid, but maybe other ways!
#keep best hit per saccver 22050
#keep best hit per sseq=29796
#keep best hit per staxid=9419, if then only keeping unique sseqs->7953 (could then map results back later,but I have not implemented this)
#keep best hit per species=8765 (think this is not ideal)
fullblast2<-fullblast[order(fullblast$pident,decreasing = T),]
#fullblast2<-fullblast2[!duplicated(fullblast2$saccver),] #testing with one per saccaver instead of taxid, not much difference
fullblast2<-fullblast2[!duplicated(fullblast2$taxids),] 

#make new header names with taxids, for building db
fullblast2$seq.name<-paste0(fullblast2$saccver," taxid=", fullblast2$taxids,"; ",fullblast2$ssciname)

#change sseq to seq.text, for phylotools
colnames(fullblast2)<-gsub("sseq","seq.text",colnames(fullblast2)) 

#add lineages for loop blasting later
fullblast2<-add.lineage.df(fullblast2,ncbiTaxDir) 

#just for testing
#fullblast2<-fullblast2[1:100,]

#make path for comparison later
fullblast2$origpathS<-paste(fullblast2$K,fullblast2$P,fullblast2$C,fullblast2$O,fullblast2$F,fullblast2$G,fullblast2$S,sep=";")
fullblast2$origpathF<-path.at.level(fullblast2$origpathS,level = "F")
fullblast2$origpathG<-path.at.level(fullblast2$origpathS,level = "G")

if(create.extended.db) {
  ##TODO
  
  phylotools::dat2fasta(fullblast2[,c("seq.name","seq.text")],"for.extended.blast.fasta")
  
  system2(command = "blastn",
          args=c("-query", "for.extended.blast.fasta", "-db",nt,"-outfmt",
                 "'6 qseqid saccver ssciname evalue staxid pident qcovs'","-evalue",1,"-num_threads", 2, "-max_target_seqs", 
                 100, "-max_hsps",1,"-word_size", 11,"-perc_identity", 50,"-qcov_hsp_perc",98,
                 "-gapopen", 0, "-gapextend", 2, "-reward", 1, "-penalty", -1, "-dust","no", "-out", "for.extended.blast.tsv"), wait = T)
  #hard coded for now, note the "taxids" rather than "staxid", which metabin does not accept
  headers<-paste0("'1i",paste(c("qseqid", "saccver", "ssciname","evalue", "taxids", "pident", "qcovs"),collapse = "\t"),"'")
  
  system2("sed",c("-i", headers, "for.extended.blast.tsv"),wait = T)
}


#export as fasta, seqs that have taxonomy at desired level
fullblast3<-list()

for(i in 1:length(TaxlevelTest)){
  
  count1<-nrow(fullblast2)
  
  fullblast3[[i]]<-fullblast2[fullblast2[,TaxlevelTest[i]]!="unknown",]
  names(fullblast3)[i]<-TaxlevelTest[i]

  #if TaxlevelTest[i]="S" remove sp.-type entries
  if(TaxlevelTest[[i]]=="S"){
    message("Removing species with 'sp.', numbers or more than one space")
    if(length(grep(" sp\\.",fullblast3[[i]]$S,ignore.case = T))>0) fullblast3[[i]]<-fullblast3[[i]][-grep(" sp\\.",fullblast3[[i]]$S,ignore.case = T),]
    if(length(grep(" .* .*",fullblast3[[i]]$S,ignore.case = T))>0) fullblast3[[i]]<-fullblast3[[i]][-grep(" .* .*",fullblast3[[i]]$S,ignore.case = T),]
    if(length(grep("[0-9]",fullblast3[[i]]$S))>0) fullblast3[[i]]<-fullblast3[[i]][-grep("[0-9]",fullblast3[[i]]$S),]
  }
  
  message(count1-nrow(fullblast3[[i]]), " sequences removed")
  
  x<-data.frame(taxon=unique(fullblast3[[1]]$origpathS),count=1)
  bas.krona.plot(taxatable = x,out = paste0("db",TaxlevelTest[[i]],".html"))
  
  phylotools::dat2fasta(fullblast3[[i]][,c("seq.name","seq.text")],paste0("tempBLASTDB.",TaxlevelTest[i],".fasta"))
  
  #make blastdb using metabinkit
  system2("metabinkit_blastgendb",c("-f",paste0("tempBLASTDB.",TaxlevelTest[i],".fasta"),"-o",
                                    paste0("tempBLASTDB.",TaxlevelTest[i]),"-c"),wait = T)
}

#######################################################
#LOOP BLAST

#loop blast for each taxonlevel
  
for(j in 1:length(TaxlevelTest)){
    
    if(TaxlevelTest[j]=="S"){
      message("species level, straight BLAST running...")
      
      out<-paste0("tempBLASTDB.",TaxlevelTest[j],".tsv")
      
      file.remove(out)
      
      system2(command = "blastn",
              args=c("-query", paste0("tempBLASTDB.",TaxlevelTest[j],".fasta"), "-db",paste0("tempBLASTDB.",TaxlevelTest[j]),"-outfmt",
                     "'6 qseqid saccver ssciname evalue staxid pident qcovs'","-evalue",1,"-num_threads", 16, "-max_target_seqs", 
                     100, "-max_hsps",1,"-word_size", 11,"-perc_identity", 50,"-qcov_hsp_perc",98,
                     "-gapopen", 0, "-gapextend", 2, "-reward", 1, "-penalty", -1, "-dust","no", "-out", out), wait = T)
      #hard coded for now, note the "taxids" rather than "staxid", which metabin does not accept
      headers<-paste0("'1i",paste(c("qseqid", "saccver", "ssciname","evalue", "taxids", "pident", "qcovs"),collapse = "\t"),"'")
      
      system2("sed",c("-i", headers, out),wait = T)
      
      
    } else {
      
      if(TaxlevelTest[j]=="G") {
        ex.seqid.group<-"S"
        out<-paste0("tempBLASTDB.",TaxlevelTest[j],".tsv")
        
        file.remove(out)
      }
      
      if(TaxlevelTest[j]=="F") {
        ex.seqid.group<-"G"
        out<-paste0("tempBLASTDB.",TaxlevelTest[j],".tsv")
        
        file.remove(out)
        }

      for(i in 1:nrow(fullblast3[[j]])){
      
        message("loop ",i)
        
        #make query fasta
        b<-fullblast3[[j]][fullblast3[[j]]$saccver==fullblast3[[j]]$saccver[i],]
        phylotools::dat2fasta(b[,c("seq.name","seq.text")],"temp.seq.fasta")
        
        #get saccvers of query group
        ex.seqids<-fullblast3[[j]][fullblast3[[j]][,ex.seqid.group]==b[,ex.seqid.group],"saccver"]  
        write.table(ex.seqids,"ex.seqids.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = F)
        
        #blast (cant use metabin_blast because of -negative_seqidlist)
        system2(command = "blastdb_aliastool",
                args=c("-seqid_file_in", "ex.seqids.txt","-seqid_file_out","ex.seqids.out.txt"), wait = T)
        
        system2(command = "blastn",
                args=c("-query", "temp.seq.fasta", "-db",paste0("tempBLASTDB.",TaxlevelTest[j]),"-outfmt",
                       "'6 qseqid saccver ssciname evalue staxid pident qcovs'","-evalue",1,"-num_threads", 16, "-max_target_seqs", 
                       100, "-max_hsps",1,"-word_size", 11,"-perc_identity", 50,"-qcov_hsp_perc",98,
                       "-gapopen", 0, "-gapextend", 2, "-reward", 1, "-penalty", -1, "-dust","no", 
                       "-negative_seqidlist", "ex.seqids.out.txt", 
                       "-out",
                       "temp.seq.blast.txt"), wait = T)
        #could exclude evalue and qcovs, might save a bit of time
        
        #store blast results
        lblast<-system2("wc",c("-l","temp.seq.blast.txt"),wait=T)
        
        #if(lblast>0) {
          blastResults<-data.table::fread("temp.seq.blast.txt",sep = "\t",data.table = F)
          write.table(blastResults,file = out,append = T,quote = F,row.names = F,sep = "\t",col.names = F)
         # }
        }
      #hard coded for now, note the "taxids" rather than "staxid", which metabin does not accept
      headers<-paste0("'1i",paste(c("qseqid", "saccver", "ssciname","evalue", "taxids", "pident", "qcovs"),collapse = "\t"),"'")
      
      system2("sed",c("-i", headers, out),wait = T)
     }
  
  #merge lineage to blast results
  a<-data.table::fread(out,data.table=F)
  a2<-merge(a,fullblast2[,c("saccver","K","P","C","O","F","G","S")],by = "saccver",all.x=T,all.y=F)
  write.table(a2,out,sep="\t",quote=F,row.names = F)
    
  message("need to add species vs family thresh blast")
}

#keep as objects for inspection
files<-c("tempBLASTDB.S.tsv","tempBLASTDB.G.tsv","tempBLASTDB.F.tsv")
bin.results<-list()
for(i in 1:length(files)){
  bin.results[[i]]<-data.table::fread(files[i],data.table = F)
}
bin.resultsS<-bin.results[[1]]
bin.resultsG<-bin.results[[2]]
bin.resultsF<-bin.results[[3]]


#######################################################

#LOOP BINNING 
outfiles<-"outfiles"
l=0

for(i in 1:length(files)){
  pidents.list2<-pidents.list
  
  ##TODO carefully think about these, I am not convinced of the zeros, maybe should just be all the same...
  ##and htf should be the same also
  ##then, bin at species e.g. 99, then G99, F99, A99...if using anything less, and some queries only have a couple of (poor) hits
  #AF will still bin at species level and results are wrong!
  
  if(files[i]=="tempBLASTDB.S.tsv") {
    pidents.list2$Gpident<-pidents.list2$Spident
    pidents.list2$Fpident<-pidents.list2$Spident
    pidents.list2$HTFpident<-pidents.list2$Spident
    
    settings<-expand.grid(pidents.list2)
    settings<-settings[settings$Gpident==settings$Spident & settings$Fpident==settings$Spident & settings$HTFpident==settings$Spident,]
  }
  
  if(files[i]=="tempBLASTDB.G.tsv") {
    pidents.list2$Spident<-pidents.list$Gpident 
    pidents.list2$Fpident<-pidents.list$Gpident
    pidents.list2$HTFpident<-pidents.list2$Gpident
    
    settings<-expand.grid(pidents.list2)
    settings<-settings[settings$Spident==settings$Gpident & settings$Fpident==settings$Gpident & settings$HTFpident==settings$Gpident,]
  }
  
  if(files[i]=="tempBLASTDB.F.tsv") {
    pidents.list2$Spident<-pidents.list$Fpident
    pidents.list2$Gpident<-pidents.list$Fpident
    pidents.list2$HTFpident<-pidents.list2$Fpident
    
    settings<-expand.grid(pidents.list2)
    settings<-settings[settings$Spident==settings$Fpident & settings$Gpident==settings$Fpident & settings$HTFpident==settings$Fpident,]
    #TODO remove extraneous
  }
  
  for(j in 1:length(tops)){
    for(k in 1:nrow(settings)){
      l=l+1
      
      outfiles[l]<-paste0(gsub("tsv","mbk-",files[i]),"top",tops[j], "-S", settings$Spident[k],"-G", settings$Gpident[k]
                        ,"-F", settings$Fpident[k])
      
      system2("metabin",c("-i",files[i], "-o",outfiles[l], "-S", settings$Spident[k],"-G", settings$Gpident[k]
                       ,"-F", settings$Fpident[k],"-A", settings$HTFpident[k], "--TopSpecies", tops[j],"--TopGenus",
                       tops[j],"--TopFamily", tops[j]
                       ,"--TopAF", tops[j],"--no_mbk"),wait=T)
    }
  }
}

#read metabin results
results<-list()

for(i in 1:length(outfiles)){
  
  results[[i]]<-data.table::fread(paste0(outfiles[i],".tsv"),data.table = F)
  
  level<-gsub("\\.","",stringr::str_extract(outfiles[i],"\\..*\\.mbk"))
  level<-gsub("mbk","",level)
  
  if(level=="S") fullblast4<-fullblast3$S
  if(level=="G") fullblast4<-fullblast3$G
  if(level=="F") fullblast4<-fullblast3$F
  
  #add no hit to each table
  results[[i]]<-merge(results[[i]], fullblast4[,c("saccver","origpathS","origpathG","origpathF")],
                      by.x = "qseqid",by.y = "saccver",all=T)
  results[[i]]$no.hits<-is.na(results[[i]]$pident)
  results[[i]]$level<-level
  results[[i]]$Spident<-gsub("S","",unlist(stringr::str_split(outfiles[i],"-"))[3])
  results[[i]]$Gpident<-gsub("G","",unlist(stringr::str_split(outfiles[i],"-"))[4])
  results[[i]]$Fpident<-gsub("F","",unlist(stringr::str_split(outfiles[i],"-"))[5])
  results[[i]]$top<-gsub("top","",unlist(stringr::str_split(outfiles[i],"-"))[2])
  
  #make path for comparison
  results[[i]]$binpathS<-paste(results[[i]]$K,results[[i]]$P,results[[i]]$C,results[[i]]$O,results[[i]]$F,
                                results[[i]]$G,results[[i]]$S,sep=";")
  results[[i]]$binpathG<-path.at.level(results[[i]]$binpathS,level = "G")
  results[[i]]$binpathF<-path.at.level(results[[i]]$binpathS,level = "F")
  
  names(results)[i]<-outfiles[i]
}

final.table<-do.call(rbind,results)

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

#change above to "above" to T/F
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

final.table$settings<-paste0(final.table$Spident,"_",final.table$Gpident,"_",final.table$Fpident,"_top",final.table$top)

write.table(final.table,final.table.out,quote = F,sep = "\t",row.names = T)

#print incorrects
incorrects<-list()
incorrects2<-list()
require(tidyverse)
for(i in 1:length(unique(final.table$level))){
  levelincor<-unique(final.table$level)[i]
  finalincor<-final.table[final.table$level==levelincor,]
  
  for(j in 1:length(unique(finalincor$settings))){
    finalsets<-unique(finalincor$settings)[j]
    finalincorset<-finalincor[(finalincor$settings==finalsets & finalincor[,paste0("incorrect",levelincor)]==TRUE),]
    if(nrow(finalincorset)>0) incorrects[[j]]<-finalincorset[,c(paste0("origpath",levelincor),paste0("binpath",levelincor))]
    incorrects<-incorrects %>% discard(is.null)
  }
  
  if(length(incorrects)>0) {
    incorrects2[[i]]<-do.call(rbind,incorrects)
    colnames(incorrects2[[i]])<-c("origpath","binpath")
  } else incorrects2[[i]]<-NULL
   
}

incorrects2<-incorrects2 %>% discard(is.null)
incorrectsdf<-do.call(rbind,incorrects2)

#summary counts

allcounts<-list()

for(j in 1:length(unique(final.table$level))){
  
  final.tableS<-final.table[final.table$level==unique(final.table$level)[j],]
  
  nsettings<-length(unique(final.tableS$settings))

  countsS<-data.frame(settings=rep("none",nsettings)
                    ,no_hits=rep(0,nsettings),
                    correct=rep(0,nsettings)
                    ,above=rep(0,nsettings),
                    incorrect=rep(0,nsettings),
                    failed=rep(0,nsettings)
                    )
  
  countsS$settings<-as.character(countsS$settings)
  
  countsS$level<-unique(final.table$level)[j]
  
  for(i in 1:length(unique(final.tableS$settings))){
    current.setting<-unique(final.tableS$settings)[i]
    countsS$settings[i]<-current.setting
    countsS$no_hits[i]<-sum(final.tableS[final.tableS$settings==current.setting,"no.hits"])
    
    final.tableSx<-final.tableS[final.tableS$no.hits==FALSE,]
                                            
    countsS$correct[i]<-sum(final.tableSx[final.tableSx$settings==current.setting,paste0("correct",unique(final.table$level)[j])],na.rm = T)
    countsS$incorrect[i]<-sum(final.tableSx[final.tableSx$settings==current.setting,paste0("incorrect",unique(final.table$level)[j])],na.rm = T)
    countsS$above[i]<-sum(final.tableSx[final.tableSx$settings==current.setting,paste0("above",unique(final.table$level)[j])],na.rm = T)
    countsS$failed[i]<-sum(final.tableSx[final.tableSx$settings==current.setting,paste0("incorrect",unique(final.table$level)[j])],na.rm = T)
  }
  
  allcounts[[j]]<-countsS
  
}

allcounts<-do.call(rbind,allcounts)

allcounts$sum<-rowSums(allcounts[,c("no_hits","correct","above","incorrect","failed")])

write.table(allcounts,counts.out,quote = F,sep = "\t",row.names = F)

#######################################################
#PLOTTING 
allcounts$settings<-paste0(allcounts$settings,allcounts$level)

longcount<-reshape2::melt(allcounts[,c("no_hits","correct","above","incorrect","failed","settings")],id.vars="settings")

plot.cols<-c("gray70","yellow4","khaki2","#E31A1C","darkturquoise","green1")
count.plot<-ggplot(data=longcount , aes(y=value, x=settings, fill=variable))+geom_bar(stat = "identity")+
  theme(legend.title = element_text(size=10), legend.text=element_text(size=10),
        axis.text.x=element_text(size=8,angle=45, hjust=1),legend.position="right",legend.direction="vertical") +
  scale_fill_manual(values = plot.cols) 
  #geom_text(aes(label = paste(round(value/total.reads*100,digits = 0),"%")), 
   #         position = position_stack(vjust = 0.5), size = 2)

count.plot

#############################################  
#plot by taxon, optional filter by taxa  

out.plot.list<-list()

final.table.list<-split(final.table,final.table$level)

for(i in 1:length(final.table.list)){

  xx<-final.table.list[[i]]
  
  levelxx<-xx[1,c("level")]
  longcount<-reshape2::melt(xx[,c("no.hits",paste0("correct",levelxx),paste0("above",levelxx),paste0("incorrect",levelxx),
                                  paste0("failed",levelxx)
                                  ,"settings","origpathS")],
                            id.vars=c("settings","origpathS"))
  
  #extract table to limit by taxon
  if(!is.null(limit.plot.to.taxon)){
    indexTax<-match(limit.plot.to.taxon[2],table = c("K","P","C","O","F","G","S"))
    longcount<-longcount[do.call(rbind,stringr::str_split(longcount$origpathS,";"))[,indexTax]==limit.plot.to.taxon[1],]
  }
  
  #plot at level
  indexn<-match(plot.at.level,table = c("K","P","C","O","F","G","S"))
  longcount$plotpath<-do.call(rbind,stringr::str_split(longcount$origpathS,";"))[,indexn]
  
  #sublabel
  if(levelxx=="S") sublabel<-"Binning outcomes if species is in database"
  if(levelxx=="G") sublabel<-"Binning outcomes if species is not in database"
  if(levelxx=="F") sublabel<-"Binning outcomes if genus is not in database"
  
  #plot
  plot.cols<-c("gray70","yellow4","khaki2","#E31A1C","darkturquoise","green1")
  count.plot<-ggplot(data=longcount , aes(y=as.numeric(as.logical(value)), x=settings, fill=variable))+geom_bar(stat = "identity")+
    theme(legend.title = element_text(size=10), legend.text=element_text(size=10),
          axis.text.x=element_text(size=8,angle=45, hjust=1),legend.position="right",legend.direction="vertical") +
  scale_fill_manual(values = plot.cols) +
    facet_wrap(~plotpath,scales = "free") + 
    ggtitle(paste0("TestLevel=",levelxx,"; limit=",limit.plot.to.taxon[1],"; PlotLevel=",plot.at.level), subtitle = sublabel)
  
  
  out.plot.list[[i]]<-count.plot
  names(out.plot.list)[i]<-levelxx
}

out.plot.list[[1]]
out.plot.list[[2]]
out.plot.list[[3]]

#also make new incorrectsdf based on plot settings
indexTax<-match(limit.plot.to.taxon[2],table = c("K","P","C","O","F","G","S"))
incorrectlimited<-incorrectsdf[as.data.frame(stringr::str_split_fixed(incorrectsdf$origpath,";",n = 7))[,indexTax]==limit.plot.to.taxon[1],]

#inspect incorrect reasons
#KM055377.1	
#100_100_100_top1
#tempBLASTDB.S.mbk-top0.1-S100-G100-F100 #rownames final.table
xz<-results$'tempBLASTDB.S.mbk-top0.1-S100-G100-F100'
xz[xz$qseqid=="KM055377.1",]  
bin.resultsG[bin.resultsG$qseqid=="KM055377.1",]
