####
message("#############
        Please check the settings below are correct. Please check the end of the log for any warnings. Immediately prior to warnings is 
        a printout of the script used
        ##############")

message("settings:
        ")
print(ls.str())

source(paste0(bastoolsDir,"master_functions.R"))
source(paste0(bastoolsDir,"bin.blast.R"))
library(processx)
library(dplyr)
obitaxoR<-ROBITaxonomy::read.taxonomy(dbname = obitaxo)
setwd(outDir)
#######################################################
#step 1 Convert Nunos table into fasta 
if("step1" %in% stepstotake){
  
  message("RUNNING STEP1 - Converting genbank s2 table into fasta")

  a<-data.table::fread(cattedDLS.s2, data.table = F)
  
  if(!is.null(taxonlimit)){
    message("Selecting only those in ",taxonlimit)
    a<-a[grep(taxonlimit,a$V7),]
  }
  
  message("Removing seqs with a lot (>27) of NNNs")
  if(length(grep("NNNNNNNNNNNNNNNNNNNNNNNNNNNN",a$V16))>0) a<-a[-grep("NNNNNNNNNNNNNNNNNNNNNNNNNNNN",a$V16),] 
  
  #make unique names, hmm lots, I dont think there should be - tell Nuno
  a$ids<-make.unique(a$V2, sep = "_")
  
  #output as fasta
  a$seq.name<-paste0(a$ids," taxid=",a$V8,";")
  a$seq.text<-a$V16
  phylotools::dat2fasta(a[,c("seq.name","seq.text")],outfile = catted_DLS)
  
  message("STEP1 COMPLETE")
}

#######################################################
#step 2 Add taxonomy and remove family=unknown
if("step2" %in% stepstotake){
  
  message("RUNNING STEP2 - Add taxonomy and remove family=unknown")
  
  #add lineage
  a<-phylotools::read.fasta(catted_DLS)
  a$taxids<-stringr::str_match(a$seq.name, "taxid=(.*?);")[,2]
  a<-add.lineage.df(a,ncbiTaxDir)
  
  #remove family=unknown sequences
  message("Removing sequences with family=unknown")
  before<-nrow(a)
  if(length(grep("unknown",a$F))>0){
    a<-a[-grep("unknown",a$F),]}
  after<-nrow(a)
  
  message(before-after," reads removed")
  
  #write to fasta
  a$seq.name<-paste0(do.call(rbind,stringr::str_split(a$seq.name," "))[,1]," taxid=",a$taxids,";")
  
  phylotools::dat2fasta(a[,c("seq.name","seq.text")],gsub(".fasta",".checked.fasta",catted_DLS))
  
  message("STEP2 COMPLETE")
}
#######################################################
#step 3 find targets 
if("step3" %in% stepstotake){
  
  message("RUNNING STEP3 - find targets")
  
  #convert to ecopcrdb
  obiconvert.Bas(infile = gsub(".fasta",".checked.fasta",catted_DLS),in_type = "fasta",
                 out = gsub(".fasta",".checked.ecopcrdb",catted_DLS),taxo = obitaxo,
                 out_type = "--ecopcrdb-output",add2existing = F)
  
  #run ecopcr with buffer option
  message("Creating references. maxL not used, min_length=10, buffer = T, max_error=",max_error_buildrefs)
  ecoPCR.Bas(Pf,Pr,ecopcrdb = gsub(".fasta",".checked.ecopcrdb",catted_DLS),max_error = max_error_buildrefs,
             min_length = 10,out = gsub(".fasta",".ecopcrResults.txt",catted_DLS),  buffer = buffer)
  
  #clean
  clean.ecopcroutput(ecopcrfile = gsub(".fasta",".ecopcrResults.txt",catted_DLS),rm.buffer = F,buffer.used = T)
  
  #build refs fasta from results, assumes buffer was used in ecopcr, assumes clean was done (with rm.buffer=F). Uses 'śequence' (buff-p-insert-p-buff)
  ecopcr2refs(ecopcrfile=gsub(".fasta",".ecopcrResults.txt.clean",catted_DLS),
               outfile = gsub(".fasta",".refs.fasta",catted_DLS),bufferecopcr = buffer,Pf = Pf,Pr = Pr,selection = NULL)
    
  message("STEP3 COMPLETE")
  
}
#######################################################
#step 4 map to refdb 
if("step4" %in% stepstotake){
  
  message("RUNNING STEP4 - map sequences that did not have primers found against sequences that had primer found")
  
  #find positions of refs on sequences that had primers found
  add.target.positions(cattedDLS.checked=gsub(".fasta",".checked.fasta",catted_DLS),refs=gsub(".fasta",".refs.fasta",catted_DLS)
                       ,out=gsub(".fasta",".checked.wPos.fasta",catted_DLS))
  
  #map sequences that did not have primers found against sequences that had primers found
  ##first make faste of sequences that did not have primers found
  a<-phylotools::read.fasta(gsub(".fasta",".checked.wPos.fasta",catted_DLS))
  ref.ids<-do.call(rbind,stringr::str_split(a$seq.name," "))[,1]
  write.table(ref.ids,quote = F,row.names = F,append = F,file = "ref.ids.tmp")
  system2(command = "seqkit", args=c("grep","-v", "-f","ref.ids.tmp",gsub(".fasta",".checked.fasta",catted_DLS))
          ,stdout = gsub(".fasta",".checked.NoPri.fasta",catted_DLS),wait = T)
  unlink("ref.ids.tmp")
  
  #make blast db of sequences that had primers found
  system2(command = "makeblastdb", args=c("-in", gsub(".fasta",".checked.wPos.fasta",catted_DLS)
                                          , "-dbtype", "nucl", "-parse_seqids","-out",gsub(".fasta",".checked.wPos.blastdb",catted_DLS)),wait=T)
  
  #map sequences that did not have primers found against sequences that had primer found
  system2(command = "blastn", args=c("-query", gsub(".fasta",".checked.NoPri.fasta",catted_DLS), "-task", "megablast","-db",
                                     gsub(".fasta",".checked.wPos.blastdb",catted_DLS),
                                     "-outfmt",'"6 qseqid sseqid stitle qlen qstart qend slen sstart send length pident qcovs sstrand"',
                                     "-num_threads", "16","-max_target_seqs", "1","-subject_besthit"),
          stdout=gsub(".fasta",".checked.NoPri.blast",catted_DLS),wait = T)
  
  message("STEP4 COMPLETE")
}

#######################################################
#step 5 check that sequences which mapped overlap target region
if("step5" %in% stepstotake){
  
  message("RUNNING STEP5 - check that sequences which mapped overlap target region, and make final ecopcrdb")
  
  #check that sequences which mapped overlap target region  
  ##read in blast results
  a<-data.table::fread(gsub(".fasta",".checked.NoPri.blast",catted_DLS),data.table = F)
  colnames(a)<-c("qseqid","sseqid", "stitle", "qlen", "qstart", "qend", "slen", "sstart", "send", "length", "pident", "qcovs", "sstrand")
  ###note that there are multiple hsps using blast settings above
  
  ##correct minus strands, change sstart to send and send to sstart
  kplus<-a[a$sstrand=="plus",]
  kminus<-a[a$sstrand=="minus",]
  kminus$sstart.tmp<-kminus$send
  kminus$send<-kminus$sstart
  kminus$sstart<-kminus$sstart.tmp
  kminus$sstart.tmp=NULL
  kboth<-rbind(kplus,kminus)
  
  #extract target start/stop to new cols
  kboth$target.start<-as.numeric(gsub(";","",gsub("target.start=","",stringr::str_extract(kboth$stitle,"target.start=(.*?);"))))
  kboth$target.end<-as.numeric(gsub(";","",gsub("target.end=","",stringr::str_extract(kboth$stitle,"target.end=(.*?);"))))
  
  #check that alignments cover whole target region
  kboth$aln.left.of.target.start<-kboth$sstart<kboth$target.start+1
  kboth$aln.right.of.target.end<-kboth$send>kboth$target.end-1
  kboth.contains.frag<-kboth[kboth$aln.left.of.target.start==T & kboth$aln.right.of.target.end==T,]
  
  #what about if sstart is left of target.start, and sstart+remaining length is right of send
  #kboth$aln.left.qend.full<-kboth$sstart-1+(kboth$qlen-kboth$qstart) - come back to this, not essential
  #what about if send is right of target.end, and send-remaining length is left of sstart - come back to this, not essential
  
  #join kboth.contains.frag and gsub(".fasta",".checked.wPos.fasta",catted_DLS) into single fasta
  ##extract those seqs first
  write.table(kboth.contains.frag$qseqid,quote = F,row.names = F,append = F,file = "ref.ids.tmp")
  system2(command = "seqkit", args=c("grep","-f","ref.ids.tmp",gsub(".fasta",".checked.fasta",catted_DLS))
          ,stdout = gsub(".fasta",".checked.mapcaught.fasta",catted_DLS),wait = T)
  unlink("ref.ids.tmp")
  
  #join 
  system2(command = "cat", args=c(gsub(".fasta",".checked.wPos.fasta",catted_DLS),gsub(".fasta",".checked.mapcaught.fasta",catted_DLS))
          ,stdout = gsub(".fasta",".checked.wPos.plus.mapcaught.fasta",catted_DLS),wait = T)
  
  
  #create final ecopcrdb
  obiconvert.Bas(infile = gsub(".fasta",".checked.wPos.plus.mapcaught.fasta",catted_DLS),
                 in_type = "fasta",out = gsub(".fasta",".checked.wPos.plus.mapcaught.ecopcrdb",catted_DLS),taxo = obitaxo,
                 out_type = "--ecopcrdb-output",add2existing = F)
  
  message("STEP5 COMPLETE")
}

#######################################################
#step 6 final ecopcr
if("step6" %in% stepstotake){
   
  message("RUNNING STEP6 - final ecopcr and add stats")
    
  #run final ecopcr without buffer option
  ecoPCR.Bas(Pf,Pr,ecopcrdb = gsub(".fasta",".checked.wPos.plus.mapcaught.ecopcrdb",catted_DLS),max_error = max_error_ecopcr,
             min_length=min_length,max_length=max_length,out = gsub(".fasta",".checked.wPos.plus.mapcaught.final.hits.txt",catted_DLS))
  
  #clean ecopcrfile
  clean.ecopcroutput(ecopcrfile = gsub(".fasta",".checked.wPos.plus.mapcaught.final.hits.txt",catted_DLS),
                     buffer.used = F,rm.buffer.insert = F,min_length = min_length,max_length = max_length)
  
  #Add stats to ecopcrfile
  message("Adding stats to ecopcroutput")
  ecopcroutput<-data.table::fread(gsub(".fasta",".checked.wPos.plus.mapcaught.final.hits.txt.clean",catted_DLS),data.table = F)
  ecopcrout.wstats<-add.stats.ecopcroutput(ecopcroutput = ecopcroutput,ncbiTaxDir = ncbiTaxDir,Ta = Ta,Pf = Pf,Pr = Pr)
  
  #add resolution using simple method
  message("Adding resolution")
  lca = aggregate(ecopcrout.wstats$taxids, by=list(ecopcrout.wstats$sequence),function(x) ROBITaxonomy::lowest.common.ancestor(obitaxoR,x))
  colnames(lca)<-c("sequence","taxids")
  lca<-add.lineage.df(lca,ncbiTaxDir = ncbiTaxDir)
  
  lca$taxon<-paste(lca$K,lca$P,lca$C,lca$O,lca$F,lca$G,lca$S,sep = ";")
  
  lca<-lca %>% select(taxon,everything())
  
  lca$res<-bas.get.ranks(lca)
  
  merged<-merge(ecopcrout.wstats,lca[,c("sequence","res")],by = "sequence",all.x = T)
  
  write.table(x = merged,file = out_mod_ecopcrout_file,append = F,quote = F,sep = "\t",row.names = F)

  message("STEP6 COMPLETE")
  
}
#######################################################
#step 7 - Make family primer bias tables
if("step7" %in% stepstotake){
  
  message("RUNNING STEP7 - Make family primer bias tables")
  #Make family primer bias tables
  #convert final ecopcrdb to tab
  system2(command = "obitab", args=c("-o",gsub(".fasta",".checked.wPos.plus.mapcaught.ecopcrdb",catted_DLS)), 
          stdout=gsub(".fasta",".checked.wPos.plus.mapcaught.ecopcrdb.tab",catted_DLS), wait = T)
  
  make.primer.bias.table(originaldbtab = gsub(".fasta",".checked.wPos.plus.mapcaught.ecopcrdb.tab",catted_DLS),
                         mod_ecopcrout_file = out_mod_ecopcrout_file,
                         out_bias_file = out_bias_file,
                         Pf = Pf, Pr = Pr,
                         obitaxoR = obitaxoR,min_length = min_length,max_length = max_length)
  
  message("STEP7 COMPLETE")
}

#######################################################
#step 8 - add counts to bias file

if("step8" %in% stepstotake){
  
  message("RUNNING STEP8 - add step counts to bias file")
  #catted DLS
  message("counting families in cattedDLS")
  cattedDLS<-count.nseqs.in.fams.in.fasta(catted_DLS,ncbiTaxDir)
  
  #after rm fams and checks
  message("counting families after checks")
  afterchecks<-count.nseqs.in.fams.in.fasta(fasta = gsub(".fasta",".checked.fasta",catted_DLS),ncbiTaxDir)
  
  ##ntaxa 
  taxa.afterchecks<-count.taxa.in.fams.in.fasta(fasta = gsub(".fasta",".checked.fasta",catted_DLS),ncbiTaxDir)
  
  #after first ecopcr
  message("counting families found in during 1st ecopcr to build refs")
  first.ecopcr<-count.nseqs.in.fams.in.fasta(gsub(".fasta",".refs.fasta",catted_DLS),ncbiTaxDir)
  
  #mean,min,max length in first ecopcr
  long.ecopcr<-data.table::fread(gsub(".fasta",".ecopcrResults.txt.clean",catted_DLS),data.table = F)  
  long.ecopcr$taxids<-long.ecopcr$taxid
  long.ecopcr<-add.lineage.df(long.ecopcr,ncbiTaxDir)
  long.ecopcr$path<-paste(long.ecopcr$K,long.ecopcr$P,long.ecopcr$C,long.ecopcr$O,long.ecopcr$F,sep = ";")
  long.ecopcr.mean<-aggregate(long.ecopcr$amplicon_length,by=list(long.ecopcr$path),FUN=mean)
  colnames(long.ecopcr.mean)<-c("path","meanL")
  long.ecopcr.min<-as.data.frame(aggregate(long.ecopcr$amplicon_length,by=list(long.ecopcr$path),FUN=min))
  colnames(long.ecopcr.min)<-c("path","minL")
  long.ecopcr.max<-as.data.frame(aggregate(long.ecopcr$amplicon_length,by=list(long.ecopcr$path),FUN=max))
  colnames(long.ecopcr.max)<-c("path","maxL")
  long.ecopcr<-cbind(Family=long.ecopcr.mean$path,first.ecopcr.meanL=long.ecopcr.mean$meanL,
                     first.ecopcr.minL=long.ecopcr.min$minL,first.ecopcr.maxL=long.ecopcr.max$maxL)

  #after mapping back
  message("counting families after mapping back")
  aftermapping<-count.nseqs.in.fams.in.fasta(gsub(".fasta",".checked.wPos.plus.mapcaught.fasta",catted_DLS),ncbiTaxDir)
  
  #merge all
  message("merging with bias table")
  merged<-merge(cattedDLS,afterchecks,by = "Family",all = T)
  colnames(merged)<-gsub("nseqs.x","starting.fasta",colnames(merged))
  colnames(merged)<-gsub("nseqs.y","after.checks",colnames(merged))
  
  merged<-merge(merged,taxa.afterchecks,by = "Family",all = T)
  colnames(merged)<-gsub("ntaxa","ntaxa.after.checks",colnames(merged))
  
  merged<-merge(merged,first.ecopcr,by = "Family",all = T)
  colnames(merged)<-gsub("nseqs",paste0("first.ecopcr.maxe=",max_error_buildrefs),colnames(merged))

  merged<-merge(merged,long.ecopcr,by = "Family",all = T)
  
  merged<-merge(merged,aftermapping,by = "Family",all = T)
  colnames(merged)<-gsub("nseqs","after_mapping_back",colnames(merged))
  
  biastemp<-data.table::fread(out_bias_file,sep = "\t")
  
  mergedbias<-merge(merged,biastemp,by.x ="Family",by.y = "in.odb",all = T)
  mergedbias$nseqs.odb=NULL
  
  #fixing taxonomy columns
  splittaxonomy<-do.call(rbind,stringr::str_split(mergedbias$Family,";"))
  mergedbias$K<-splittaxonomy[,1]
  mergedbias$P<-splittaxonomy[,2]
  mergedbias$C<-splittaxonomy[,3]
  mergedbias$O<-splittaxonomy[,4]
  mergedbias$F<-splittaxonomy[,5]
  
  mergedbias$path<-mergedbias$Family
  mergedbias$Family=NULL
  
  mergedbias<-mergedbias %>% select(path,K,P,C,O,F,everything())
  
  write.table(x = mergedbias,file = out_bias_file,quote = F,row.names = F,sep = "\t")
  
  message("STEP8 COMPLETE")
  
}

message("script used:
        ")
print(readLines(script))

warnings()
