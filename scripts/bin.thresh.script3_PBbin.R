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

bastoolsDir<-"/home/bastian.egeter/git_bastools/bastools/" #needed only for add.lineage.df

#TODO do a check if species not in db, do we get to family level
#TODO dont run top*pident combos that dont make sense
#TODO count "failed binning thresholds" rather than lumping it into "above"

outDir<-"/home/tutorial/temp/"
selfblastout<-"16S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.blast.filt.tempBLASTDB.tsv" #example
#using output from barcode gap report

bastoolsDir<-"/home/tutorial/TOOLS/bastools/" #needed only for add.lineage.df
ncbiTaxDir<-"/home/tutorial/TOOLS/metabinkit.install/db/"
TaxlevelTest=c("S","G","F")
known_flags<-"/media/sf_Documents/WORK/CIBIO/temp/known_flags_downloaded_22-7-20.txt"

#binning settings to loop through
tops<-c(0.1,1,10,100)
pidents.list<-list(Spident=c(100,99),Gpident=c(98,97),Fpident=c(96,70)) 

#plotting
plot.at.level<-"F"
limit.plot.to.taxon<-c("Mammalia","C")

#saving main results (can be fed to plotting chunk)
counts.out<-"PBthresh_counts.tsv"
final.table.out<-"PBthresh_final.table.tsv"

source(paste0(bastoolsDir,"master_functions.R"))
library(ggplot2)
library(tidyverse)


setwd(outDir)
#######################################################
sb<-data.table::fread(selfblastout,data.table = F)
hierarchy<-c("K","P","C","O","F","G","S")
#ok, so this is a more or less complete comparison of everything in blast database,
#can I run "fake loop blast" by just excluding results taxon by taxon?
#make fake loop blast files
#how about just remove all the non-species level entries and go from there, easier
sb$origtax<-stringr::str_split(sb$origpath,";",simplify = T)[,7]
sb$hittax<-stringr::str_split(sb$hitpath,";",simplify = T)[,7]

sb<-sb[sb$origtax!="unknown",]
sb<-sb[sb$hittax!="unknown",]

message("Removing species with 'sp.', numbers or more than one space")

if(length(grep(" sp\\.",sb$origtax,ignore.case = T))>0) sb<-sb[-grep(" sp\\.",sb$origtax,ignore.case = T),]
if(length(grep(" .* .*",sb$origtax,ignore.case = T))>0) sb<-sb[-grep(" .* .*",sb$origtax,ignore.case = T),]
if(length(grep("[0-9]",sb$origtax))>0) sb<-sb[-grep("[0-9]",sb$origtax),]

if(length(grep(" sp\\.",sb$hittax,ignore.case = T))>0) sb<-sb[-grep(" sp\\.",sb$hittax,ignore.case = T),]
if(length(grep(" .* .*",sb$hittax,ignore.case = T))>0) sb<-sb[-grep(" .* .*",sb$hittax,ignore.case = T),]
if(length(grep("[0-9]",sb$hittax))>0) sb<-sb[-grep("[0-9]",sb$hittax),]

#remove crappy hits 

message("Removing species containing the terms: uncultured, environmental, 
            unidentified,fungal, eukaryote, unclassified, synthetic")

crap<-c("uncultured","environmental","unclassified","unidentified","fungal ","eukaryote","synthetic")

for(i in 1:length(crap)){
  if(length(grep(crap[i],sb$origtax,ignore.case = T))>0) sb<-sb[-grep(crap[i],sb$origtax,ignore.case = T),]
  if(length(grep(crap[i],sb$hittaxid,ignore.case = T))>0) sb<-sb[-grep(crap[i],sb$hittaxid,ignore.case = T),]
}

lineage<-as.data.frame(stringr::str_split(sb$hitpath,";",simplify = T))
colnames(lineage)<-c("K","P","C","O","F","G","S")

###using this point as bin blast, could be earlier

sb<-sb[,c("origseqid","hitseqid", "pident","hittaxid","origpath","origseq")]
colnames(sb)<-c("qseqid","saccver", "pident","taxids","origpath","origseq")
sb<-cbind(sb,lineage)
sb$origgenus<-stringr::str_split(sb$origpath,";",simplify = T)[,6]
sb$origspecies<-stringr::str_split(sb$origpath,";",simplify = T)[,7]
sb$qseqid<-paste0(sb$qseqid,"_S")

sbforpaths<-sb[!duplicated(sb$qseqid),c("qseqid","origpath")]

for(h in 1:length(TaxlevelTest)){
  #genus should be, for each query remove the hits to the same species
  #family should be, for each query remove hits to same genus
  
  #problems with this approach, how many queries disappear because there are no hits that did not belong to that genus, looks like a lot
  #so these are not being tested at all
  
  #soooo, back to thresh blast, can I make identical output to sb2
  
  if(TaxlevelTest[h]=="G") {
    threshold.bin.blast2(df = sb,TaxlevelTest = "G",qseqidcol = "qseqid",qseqcol = "origseq",taxidcol = "taxids")
    sb.tbbG<-data.table::fread("tempBLASTDB.G.tsv",data.table = F)
    sb.tbbG<-add.lineage.df(sb.tbbG,ncbiTaxDir = ncbiTaxDir) #for metabin
    sb.tbbG<-merge(sb.tbbG,sbforpaths,all.x=T,all.y = F,by = "qseqid")
    sb.tbbG$origgenus<-stringr::str_split(sb.tbbG$origpath,";",simplify = T)[,6]
    sb.tbbG$origspecies<-stringr::str_split(sb.tbbG$origpath,";",simplify = T)[,7]
    sb.tbbG$origseq<-sb.tbbG$qseq
    sb.tbbG<-sb.tbbG[,c("qseqid","saccver","pident","taxids", "origpath","origseq", "K","P","C" ,"O","F","G","S","origgenus","origspecies")]
    sb.tbbG$qseqid<-gsub("_S$","_G",sb.tbbG$qseqid)
  }
  
  if(TaxlevelTest[h]=="F") {
    threshold.bin.blast2(df = sb,TaxlevelTest = "F",qseqidcol = "qseqid",qseqcol = "origseq",taxidcol = "taxids")
    sb.tbbF<-data.table::fread("tempBLASTDB.F.tsv",data.table = F)
    sb.tbbF<-add.lineage.df(sb.tbbF,ncbiTaxDir = ncbiTaxDir)
    sb.tbbF<-merge(sb.tbbF,sbforpaths,all.x=T,all.y = F,by = "qseqid")
    sb.tbbF$origgenus<-stringr::str_split(sb.tbbF$origpath,";",simplify = T)[,6]
    sb.tbbF$origspecies<-stringr::str_split(sb.tbbF$origpath,";",simplify = T)[,7]
    sb.tbbF<-sb.tbbF[,c("qseqid","saccver","pident","taxids", "origpath","K","P","C" ,"O","F","G","S","origgenus","origspecies")]
    sb.tbbF$origseq<-sb.tbbF$qseq
    sb.tbbF<-sb.tbbF[,c("qseqid","saccver","pident","taxids", "origpath","origseq", "K","P","C" ,"O","F","G","S","origgenus","origspecies")]
    sb.tbbF$qseqid<-gsub("_S$","_F",sb.tbbF$qseqid)
  }
  
  if(TaxlevelTest[h]=="S") sb.splS<-sb #else sb.spl<-split(sb,sb$qseqid)
  
  # if(TaxlevelTest[h]=="G") {
  #   sb.splG<-lapply(sb.spl,function(x) x<-x[x$origspecies!=x$S,])
  #   sb.splG<-do.call(rbind,sb.splG)
  #   sb.splG$qseqid<-gsub("_S$","_G",sb.splG$qseqid)
  # }
  # 
  # if(TaxlevelTest[h]=="F") {
  #   sb.splF<-lapply(sb.spl,function(x) x<-x[x$origgenus!=x$G,])
  #   sb.splF<-do.call(rbind,sb.splF)
  #   sb.splF$qseqid<-gsub("_S$","_F",sb.splF$qseqid)
  # }
}

sb2<-rbind(sb.splS,sb.tbbG,sb.tbbF)

write.table(sb2,"tempBLASTDB.ALL.tsv",quote = F,row.names = F,sep = "\t")

#loop bin at all levels
#binning settings to loop through
tops<-c(0,1,2,10,100)
#order=S,G,F,AF
pidents.list<-list(strict=c(99,97,95,90),medium=c(98,96,93,88),relaxed=c(93,85,81,80)) 

binfile.list<-list()
countloop<-0
for(j in 1:length(tops)){
  for(k in 1:length(pidents.list)){
    binfile<-paste0("top_",tops[j],".pidents_",paste(pidents.list[[k]],collapse = "."))
    countloop<-countloop+1
    binfile.list[[countloop]]<-binfile
    argsmbk<-c("-i","tempBLASTDB.ALL.tsv", "-o",binfile,"--no_mbk")
    argsmbk<-c(argsmbk,"-S", pidents.list[[k]][1],"-G", pidents.list[[k]][2],"-F", pidents.list[[k]][3],"-A", pidents.list[[k]][4],
               "--TopSpecies", tops[j],"--TopGenus",
               tops[j],"--TopFamily", tops[j]
               ,"--TopAF", tops[j])
    if(!is.null(known_flags)) argsmbk<-c(argsmbk,"--FilterFile",known_flags,"--FilterCol","saccver")
    system2("metabin",argsmbk, wait=T)
  }
}

outfiles<-do.call(c,binfile.list)

#set metabin args
# if(!is.null(SpeciesBL)) argsmbk<-c(argsmbk,"--SpeciesBL",SpeciesBL)
# if(!is.null(GenusBL)) argsmbk<-c(argsmbk,"--GenusBL",GenusBL)
# if(!is.null(FamilyBL)) argsmbk<-c(argsmbk,"--FamilyBL",FamilyBL)
#######################################################

#read metabin results and merge with sb
results<-list()
sb3<-sb2[!duplicated(sb2$qseqid),]

for(i in 1:length(outfiles)){
  a<-data.table::fread(paste0(outfiles[i],".tsv"),data.table = F)
  a$binpath<-paste(a$K,a$P,a$C,a$O,a$F,a$G,a$S,sep = ";")
  sb4<-merge(sb3[,c("qseqid","origpath","saccver")],a[,c("qseqid","binpath")],all = T,by = "qseqid")
  sb4$settings<-outfiles[i]
  results[[i]]<-sb4
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

#no hits?

  # #add no hit to each table
  # results[[i]]<-merge(results[[i]], fullblast4[,c("saccver","origpathS","origpathG","origpathF")],
  #                     by.x = "qseqid",by.y = "saccver",all=T)
  # results[[i]]$no.hits<-is.na(results[[i]]$pident)
  # results[[i]]$level<-level
  # results[[i]]$Spident<-gsub("S","",unlist(stringr::str_split(outfiles[i],"-"))[3])
  # results[[i]]$Gpident<-gsub("G","",unlist(stringr::str_split(outfiles[i],"-"))[4])
  # results[[i]]$Fpident<-gsub("F","",unlist(stringr::str_split(outfiles[i],"-"))[5])
  # results[[i]]$top<-gsub("top","",unlist(stringr::str_split(outfiles[i],"-"))[2])
  # 
 

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

write.table(final.table,final.table.out,quote = F,sep = "\t",row.names = F)

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
    if(nrow(finalincorset)>0) incorrects[[j]]<-finalincorset[,c(paste0("origpath",levelincor),paste0("binpath",levelincor),"qseqid","saccver")]
    incorrects<-incorrects %>% discard(is.null)
  }
  
  if(length(incorrects)>0) {
    incorrects2[[i]]<-do.call(rbind,incorrects)
    colnames(incorrects2[[i]])<-c("origpath","binpath","qseqid","saccver")
  } else incorrects2[[i]]<-NULL
  
}

incorrects2<-incorrects2 %>% discard(is.null)
incorrectsdf<-do.call(rbind,incorrects2)

#summary counts

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
