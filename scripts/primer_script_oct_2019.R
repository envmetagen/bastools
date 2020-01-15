setwd("/home/bastian.egeter/git_bastools/bastools/")
file.sources<-c("add.taxids.fasta.BAS.R","building_dbs.R","manip_fasta.R","running_ecoTools.R","build_refs.R",
                "mapTrim.R","primer_mistmatches.R","primer_matches_by_family.R","building_dbs4.R","ecopcr2refs.R",
                "make.primer.bias.tables.R","make.primer.bias.tables.addfuncs.R",
                "make.primer.bias.tables.famfuncs.R","counting.fams.at.each.step.R")
sapply(file.sources,source)
library(processx)
library(dplyr)

setwd(outDir)

#assumes empty starting folder
#download gb files

#step1 - download
if("step1" %in% stepstotake){
  
  message("RUNNING STEP1")
  
  for(i in 1:length(group.taxids)){
    DL.nuccore.gb(group.taxid = group.taxids[i],gene = gene,ncbiTaxDir = ncbiTaxDir,rank_req = "family")
  }
  
  count.download<-as.data.frame(bascount.gb.recur("."))
  colnames(count.download)<-"count.download"
  write.table(count.download,stepcountfile,quote = F,sep = "\t",row.names = F)
  
  message("STEP1 COMPLETE")
}

#######################################################
#Step2 - extract genes
if("step2" %in% stepstotake){
  
  message("RUNNING STEP2")
  
  #first remove files with empty result
  a<-system2("grep",args = c("'Empty result - nothing to do'",list.files(pattern = "*.gb")),wait = T,
             stderr = T,stdout = T)
  unlink(gsub(":\t.*","",a))
  a<-system2("grep",args = c("'Unable to obtain query'",list.files(pattern = "*.gb")),wait = T,
             stderr = T,stdout = T)
  unlink(gsub(":\t.*","",a))
  
  #extract gene
for(i in 1:length(list.files(pattern = "*.gb"))){
  a<-list.files(pattern = "*.gb")[i]
  extract.gene.gb(gbfile = a,gene = gene)
}

#cat (assumes otherwise empty folder)
catted_DLS<-paste0(paste(groups,collapse = "_"),"_",paste(group.taxids,collapse = "_"),"_",gene,"_nuccore_",
  Sys.Date(),".fasta")
system2("cat", args = c(list.files(pattern = "*.extract.fasta")), stdout = catted_DLS,wait = T)

#count
count.afterextract<-as.data.frame(bascount.fasta(catted_DLS))
colnames(count.afterextract)<-"count.afterextract"
write.table(count.afterextract,stepcountfile,quote = F,sep = "\t",row.names = F,append = T)

message("STEP2 COMPLETE")
}

#######################################################
#Step2a - not extracting genes
if("step2a" %in% stepstotake){
  
  message("RUNNING STEP2a")
  
  #first remove files with empty result
  a<-system2("grep",args = c("'Empty result - nothing to do'",list.files(pattern = "*.gb")),wait = T,
             stderr = T,stdout = T)
  unlink(gsub(":\t.*","",a))
  a<-system2("grep",args = c("'Unable to obtain query'",list.files(pattern = "*.gb")),wait = T,
             stderr = T,stdout = T)
  unlink(gsub(":\t.*","",a))
  
  #convert to fasta
  for(i in 1:length(list.files(pattern = "*.gb"))){
    a<-list.files(pattern = "*.gb")[i]
    gb2fasta(gbfile = a)
  }
  
  #cat (assumes otherwise empty folder)
  catted_DLS<-paste0(paste(groups,collapse = "_"),"_",paste(group.taxids,collapse = "_"),"_",gene,"_nuccore_",
                     Sys.Date(),".fasta")
  system2("cat", args = c(list.files(pattern = "*.fasta")), stdout = catted_DLS,wait = T)
  
  #count
  count.afterextract<-as.data.frame(bascount.fasta(catted_DLS))
  colnames(count.afterextract)<-"count.after.convert"
  write.table(count.afterextract,stepcountfile,quote = F,sep = "\t",row.names = F,append = T)
  
  message("STEP2a COMPLETE")
  
}
#######################################################
#step 3 remove min L and add taxonomy
if("step3" %in% stepstotake){
  
  message("RUNNING STEP3")
  
  catted_DLS<- list.files(pattern = paste0(paste(groups,collapse = "_")))
    
  #remove any seqs less than minimum required length
  f<-process$new(command = "obigrep", 
                 args = c("-l",sum(min_length,nchar(Pf),nchar(Pr),buffer,buffer), catted_DLS), echo_cmd = T,
                 stdout=gsub(".fasta",".minL.fasta",catted_DLS))
  f$wait()
  
  count.after.rm.minL<-as.data.frame(bascount.fasta(gsub(".fasta",".minL.fasta",catted_DLS)))
  colnames(count.after.rm.minL)<-"count.after.rm.minL"
  write.table(count.after.rm.minL,stepcountfile,quote = F,sep = "\t",row.names = F,append = T)
  
  #add lineage
  add.lineage.fasta.BAS(infasta=gsub(".fasta",".minL.fasta",catted_DLS),taxids=T,
                        ncbiTaxDir=ncbiTaxDir,out=gsub(".fasta",".lin.minL.fasta",catted_DLS))
  count.after.add.lin<-as.data.frame(bascount.fasta(gsub(".fasta",".lin.minL.fasta",catted_DLS)))
  colnames(count.after.add.lin)<-"count.after.add.lin"
  write.table(count.after.add.lin,stepcountfile,quote = F,sep = "\t",row.names = F,append = T)
  
  #some checks & removals
  #1. remove unknown kingdoms
  fastatemp<-phylotools::read.fasta(gsub(".fasta",".lin.minL.fasta",catted_DLS))
  if(length(grep("kingdom=unknown",fastatemp$seq.name))>0){
    fastatemp<-fastatemp[-grep("kingdom=unknown",fastatemp$seq.name),]}
  #2. further remove unknown families
  if(length(grep("family=unknown",fastatemp$seq.name))>0){
    fastatemp<-fastatemp[-grep("family=unknown",fastatemp$seq.name),]}
  
  
  ####################COME BACK TO THIS
  #3. do all sequences belong to expected group
  # if(length(grep(paste0(group.rank,"=",group),fastatemp$seq.name))>0){
  #   fastatemp<-fastatemp[grep(paste0(group.rank,"=",group),fastatemp$seq.name),]}
  
  phylotools::dat2fasta(fastatemp,gsub(".fasta",".checked.lin.minL.fasta",catted_DLS))
  
  count.after.rmNoFams<-as.data.frame(bascount.fasta(gsub(".fasta",".checked.lin.minL.fasta",catted_DLS)))
  colnames(count.after.rmNoFams)<-"count.after.rmNoFams"
  count.familes.after.rmNoFams<-as.data.frame(length(unique(stringr::str_match(fastatemp$seq.name, 
                                                                               "family=(.*?);")[,2])))
  colnames(count.familes.after.rmNoFams)<-"count.familes.after.rmNoFams"
  
  write.table(count.after.rmNoFams,stepcountfile,quote = F,sep = "\t",row.names = F,append = T)
  write.table(count.familes.after.rmNoFams,stepcountfile,quote = F,sep = "\t",row.names = F,append = T)
  
  message("STEP3 COMPLETE")

############actually should derep within taxid, as I was doing before
}
#######################################################
#step 4 make refdb 
if("step4" %in% stepstotake){
  
  message("RUNNING STEP4")
  
#just for ease for now
file.copy(list.files(pattern = "*.checked.lin.minL.fasta"),to = "formatted.minL.lineage.goodfam.fasta")

#ensure uniq ids
system2(command = "obiannotate", args=c("--uniq-id","formatted.minL.lineage.goodfam.fasta"), 
        stdout="formatted.minL.lineage.goodfam.uid.fasta",wait=T)
bascount.fasta("formatted.minL.lineage.goodfam.uid.fasta")


######################should add this chunk to obiconvert.BAS (with option for adding to existing)
#remove "|"s to stop ecopcr errors
f<-process$new(command = "sed", args = c("s/|/_/g",
                                         "formatted.minL.lineage.goodfam.uid.fasta"), echo_cmd = T,
               stdout="formatted.minL.lineage.goodfam.uid2.fasta")
f$wait()
bascount.fasta("formatted.minL.lineage.goodfam.uid2.fasta")

####################################
#convert to ecopcrdb
unlink(x = "formatted.minL.lineage.goodfam.uid2.ecopcrdb*")

obiconvert.Bas(infile = "formatted.minL.lineage.goodfam.uid2.fasta",in_type = "fasta",
               out = "formatted.minL.lineage.goodfam.uid2.ecopcrdb",taxo = obitaxo,
               out_type = "--ecopcrdb-output" )

#run ecopcr with buffer option
ecoPCR.Bas(Pf,Pr,ecopcrdb = "formatted.minL.lineage.goodfam.uid2.ecopcrdb",max_error = max_error_buildrefs,
           min_length,max_length,out = "all.ecopcr.hits.txt",  buffer = buffer)

#build refs fasta from results
ecopcr2refs2(ecopcrfile="all.ecopcr.hits.txt",outfile = "mapping.reference.fasta",bufferecopcr = buffer)
mapping.reference<-as.data.frame(bascount.fasta("mapping.reference.fasta"))
colnames(mapping.reference)<-"mapping.reference"
write.table(mapping.reference,stepcountfile,quote = F,sep = "\t",row.names = F,append = T)

message("STEP4 COMPLETE")

}
#######################################################
#step 5 map to refdb and make final db
if("step5" %in% stepstotake){
  
  message("RUNNING STEP5")
  
#MAP ALL SEQUENCES AGAINST TARGET REFERENCE DATABASE
map2targets(queries.to.map = "formatted.minL.lineage.goodfam.uid2.fasta",
            refs = "mapping.reference.fasta",out = "formatted.minL.lineage.goodfam.uid2.mapped.txt")

#trim
mapTrim2.simple(query = "formatted.minL.lineage.goodfam.uid2.fasta",
                blast.results.file = "formatted.minL.lineage.goodfam.uid2.mapped.txt",
                out = "formatted.minL.lineage.goodfam.uid2.mapped.fasta",qc=0.97)

count.aftermap<-as.data.frame(bascount.fasta("formatted.minL.lineage.goodfam.uid2.mapped.fasta"))
colnames(count.aftermap)<-"count.aftermap"
write.table(count.aftermap,stepcountfile,quote = F,sep = "\t",row.names = F,append = T)

#count.families
tempfasta<-phylotools::read.fasta("formatted.minL.lineage.goodfam.uid2.mapped.fasta")

count.familes.after.map<-as.data.frame(length(unique(stringr::str_match(tempfasta$seq.name, "family=(.*?);")[,2])))
colnames(count.familes.after.map)<-"count.familes.after.map"
write.table(count.familes.after.map,stepcountfile,quote = F,sep = "\t",row.names = F,append = T)

#create final ecopcrdb
unlink(x = "final.ecopcrdb*")
obiconvert.Bas(infile = "formatted.minL.lineage.goodfam.uid2.mapped.fasta",
               in_type = "fasta",out = "final.ecopcrdb",taxo = obitaxo,
               out_type = "--ecopcrdb-output")

message("STEP5 COMPLETE")

}
#######################################################
#step 6 - run final ecopcr and do stats
if("step6" %in% stepstotake){
  
  message("RUNNING STEP6")
  
#run final ecopcr without buffer option
ecoPCR.Bas(Pf,Pr,ecopcrdb = "final.ecopcrdb",max_error = max_error_ecopcr,
           min_length,max_length,out = "final.ecopcr.hits.txt")
#convert final ecopcrdb to tab
system2(command = "obitab", args=c("-o","final.ecopcrdb"), stdout="final.ecopcrdb.tab", wait = T)
#DO STATS
make.primer.bias.tables(originaldbtab = "final.ecopcrdb.tab",ecopcrfile = "final.ecopcr.hits.txt",
                        out_bias_file = out_bias_file,
                        out_mod_ecopcrout_file = out_mod_ecopcrout_file,Pf = Pf, Pr = Pr,
                        obitaxoR = obitaxoR,min_length = min_length,max_length = max_length)

biastemp<-data.table::fread(out_bias_file,sep = "\t")

count.amped<-as.data.frame(sum(biastemp$nseqs.amped,na.rm = T))
colnames(count.amped)<-"count.amped"
families.amped<-as.data.frame(sum(biastemp$amplified))
colnames(families.amped)<-"families.amped"
count.uniq.brcds.amped<-as.data.frame(sum(biastemp$n.uniq.brcds.amped,na.rm = T))
colnames(count.uniq.brcds.amped)<-"count.uniq.brcds.amped"

write.table(count.amped,stepcountfile,quote = F,sep = "\t",row.names = F,append = T)
write.table(families.amped,stepcountfile,quote = F,sep = "\t",row.names = F,append = T)
write.table(count.uniq.brcds.amped,stepcountfile,quote = F,sep = "\t",row.names = F,append = T)

message("STEP6 COMPLETE")

}

#######################################################
#step 7 - add counts to bias file

if("step7" %in% stepstotake){
  
  message("RUNNING STEP7")
  biastemp<-add.counts.to.biasfile(ncbiTaxDir = ncbiTaxDir,download.fasta = catted_DLS,after.minL.fasta = gsub(".fasta",".minL.fasta",catted_DLS)
                          ,after.checks.fasta = gsub(".fasta",".checked.lin.minL.fasta",catted_DLS)
                         ,first.ecopcr.hit.table = "final.ecopcr.hits.txt",mapped.fasta = "formatted.minL.lineage.goodfam.uid2.mapped.fasta"
                         ,out_bias_file = out_bias_file)
  
  write.table(x = biastemp,file = out_bias_file,quote = F,row.names = F,sep = "\t")

}

warnings()