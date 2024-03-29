message("settings:
        ")
print(ls.str())

library(processx)
library(dplyr)
source(paste0(bastoolsDir,"master_functions.R"))
setwd(outDir)

####################################################
#step 7 blast
if("step7" %in% stepstotake){
  
  message("STEP7-BLAST")
  
  #BLASTING NT
  if(class(startingfastas)=="data.frame") {
    stop ("not written this way for pipe3")
  } else{
    for(i in 1:length(startingfastas)){
    blast.status<-blast.min.bas2(infasta = startingfastas[i],refdb = refdb,blast_exec = blast_exec, wait = T,
                               taxidlimit = taxidlimit[[i]], ncbiTaxDir = ncbiTaxDir,opts = opts,overWrite = T)
    }
  }
  
  message("STEP7 complete")
  
}

####################################################
#step 8 filter 
if("step8" %in% stepstotake){  
  
  message("STEP8-filter blast")
  
  #startingfastas<-c("12S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.fasta"
   #                 ,"16S.none.flash2.vsearch_qfilt.cutadapt.vsearch_uniq.vsearch_afilt.allsamples_step5.ALL_vsearch_uniq.nodenoise.noclust.fasta")
  #outDir="/mnt/Disk1/BASTIAN_POST_MBC_MISEQS/2018_07/"
  
  
  #FILTER BLASTS
  files<-gsub(".fasta$",".blast.txt",startingfastas)
  
  for(i in 1:length(files)){
    message(paste("filtering blast results for",files[i]))
    blastfile = files[i]
    out<-gsub(".blast.txt",".blast.filt.txt",files[i])
    filter.blast3(blastfile = blastfile,ncbiTaxDir = ncbiTaxDir,out = out)
  }
  message("STEP8 complete")
  
}
####################################################

#step 8a find errors 
if("step8a" %in% stepstotake){  
  
  message("STEP8a-finding database errors")
  
  files<-gsub(".fasta$",".blast.filt.txt",startingfastas)
  
  for(i in 1:length(files)){
    message(paste("finding errors for ",files[i]))
    message("excluding all XM_,XR_,XP_ accessions. These are in-silico generated ")
  
    #outfiles
    input<-paste0(outDir,files[i])
    krona_html_db=gsub(".blast.filt.txt",".blast.filt.database.html",files[i])
    selfblastout=gsub(".blast.filt.txt",".blast.filt.tempBLASTDB.tsv",files[i])
    flagged_accessions=gsub(".blast.filt.txt",".blast.filt.flagged.tsv",files[i])
    flagged_error_detailed_table=gsub(".blast.filt.txt",".blast.filt.flagged.details.tsv",files[i])
    further_potential_errors=gsub(".blast.filt.txt",".blast.filt.non-flagged.potential.errors.tsv",files[i])
    out_html=gsub(".blast.filt.txt",".blast.filt.barcode.gap.report.html",files[i])
    
    divergence<-rep(5,length(TaxlevelTestall))
    use_flagged_accessions_bcg=T
    rm.low.read.queries=NULL
    #limit plots for calc_bar_gap
    plot.limit.taxon=NULL #";Streptophyta;"
    
    #thresher
    use_flagged_accessions_mbk=T
    if(!exists("plot.at.level")) {
      plot.at.level<-"C" #for thresher
      message("Defaulting to class plot level")
    }
    
    if(!exists("limit.plot.to.taxon")) {
      limit.plot.to.taxon<-NULL #for thresher
      message("Not limiting plots by taxon")
    } 
    
    #saving main results (can be fed to plotting chunk)
    counts.out<-gsub(".blast.filt.txt",".blast.filt.thresher.counts.tsv",files[i])
    final.table.out<-gsub(".blast.filt.txt",".blast.filt.thresher.final.table.tsv",files[i])
    
    #knit
    rmarkdown::render(input = paste0(bastoolsDir,"scripts/barcode_gap_report.Rmd"),output_file = paste0(outDir,out_html))
    
    ##TODO what to do with flagged accessions? Make a master file that I will use for everything? No, I think just save and can add 
    #to google sheet manually if desired
  }
  message("STEP8a complete")
  
}


####################################################
#step 9 BIN
if("step9" %in% stepstotake){ 
  
  message("STEP9 - bin")
  
  files<-gsub(".fasta$",".blast.filt.txt",startingfastas)
  
  if(use.metabin){
    
    for(i in 1:length(files)){
      message(paste("binning filtered blast results for",files[i]))
      filtered_blastfile<-files[i]
      binfile<-gsub(".blast.filt.txt",".bins",files[i])
      
      #make FilterFile
      if(use_flagged_accessions_mbk){
        if(file.exists(gsub(".blast.filt.txt",".blast.filt.flagged.tsv",files[i]))) {
          flagged_accessions<-gsub(".blast.filt.txt",".blast.filt.flagged.tsv",files[i])
          flags_step8a<-data.table::fread(flagged_accessions,data.table = F,header = F)
        }
        
        if(!is.null(known_flags)) flags_pre_defined<-data.table::fread(known_flags,data.table = F,header = F)
        
        if(exists("flags_step8a")) if(exists("flags_pre_defined")) {
          message("Combining pre-defined flags and flags found in step8a")
          all_flags<-rbind(flags_pre_defined,flags_step8a)
          }
        if(exists(x = "flags_step8a")) if(!exists("flags_pre_defined")) {
          message("Using only flags found in step8a")
          all_flags<-flags_step8a 
        }
        if(!exists(x = "flags_step8a")) if(exists("flags_pre_defined")) {
          message("Using pre-defined flags only")
          all_flags<-flags_pre_defined
        }
        
        FilterFile<-gsub(".blast.filt.txt",".blast.filt.flagged_plus_known.tsv",files[i])
        write.table(all_flags,FilterFile,quote = F,row.names = F,col.names = F)
        
      } else {
        FilterFile<-NULL
        message("Not applying any accession filters")
      }
      
      if(exists("top")){
        message("Overwriting all tops with ", top)
        topS<-top
        topG<-top
        topF<-top
        topAF<-top
      }
      
      #set metabin args
       argsmbk<-c("-i",filtered_blastfile, "-o",binfile, "-S", spident,"-G", gpident
                  ,"-F", fpident,"-A", abspident, "--TopSpecies", topS,"--TopGenus",
                  topG,"--TopFamily", topF
                  ,"--TopAF", topAF,"--no_mbk","--sp_discard_sp", "--sp_discard_mt2w","--sp_discard_num")
       
       if(!is.null(FilterFile)) argsmbk<-c(argsmbk,"--FilterFile",FilterFile,"--FilterCol","saccver")
       if(!is.null(SpeciesBL)) argsmbk<-c(argsmbk,"--SpeciesBL",SpeciesBL)
       if(!is.null(GenusBL)) argsmbk<-c(argsmbk,"--GenusBL",GenusBL)
       if(!is.null(FamilyBL)) argsmbk<-c(argsmbk,"--FamilyBL",FamilyBL)
       
       system2("metabin",argsmbk, wait=T)
      }
    } else {
  
    for(i in 1:length(files)){
      message(paste("binning filtered blast results for",files[i]))
      filtered_blastfile<-files[i]
      binfile<-gsub(".blast.filt.txt",".bins.txt",files[i])
      bin.blast2(filtered_blastfile = filtered_blastfile,ncbiTaxDir = ncbiTaxDir,
                 out = binfile,spident = spident,gpident = gpident,
                 fpident = fpident,abspident = abspident)
    }
  }
  message("STEP9 complete")
  
} 

####################################################
#step 10 merge with otutab

if("step10" %in% stepstotake){  
  
  message("STEP10 - merge with otutab")
  
  files<-gsub(".fasta$",".otutab.tsv",startingfastas)
  
  for(i in 1:length(files)){
    otutabfile<-files[i]
    if(use.metabin) {
      binfile<-gsub(".otutab.tsv",".bins.tsv",files[i])
    } else {
      binfile<-gsub(".otutab.tsv",".bins.txt",files[i])
    }
    out<-gsub(".otutab.tsv",".taxatable.txt",files[i])
    
    MBC_otu_bin_blast_merge(MBC_otutab = otutabfile,bin_blast_results = binfile,out = out)
  }
  message("STEP10 complete")
  
}

####################################################
#step 11 apply taxon filter

if("step11" %in% stepstotake){  
  
  message("STEP11 - taxon filter")
  message("CHANGE SCRIPT TO ONLY TAKE RELEVANT FILES")
  
  files<-gsub(".fasta$",".taxatable.txt",startingfastas)
  
  taxon.filter.solo(files,filterpc)
  
  message("STEP11 complete")
  
}

####################################################
#step 12 splice taxa tables and produce OTU.per.bin file

if("step12" %in% stepstotake){  
  
  message("STEP12 - splice tables and produce OTU.per.bin and OTU.per.sample files")
  #message("CHANGE SCRIPT TO ONLY TAKE RELEVANT FILES")
  
  files<-gsub(".fasta$",".taxatable.tf.txt",startingfastas)
  
  splice.taxatables(files)
  
  ####output OTUs.per.bin
  files<-list.files(pattern = ".taxatable.tf.spliced.txt$")
  for(i in 1:length(files)){
    taxatab<-data.table::fread(files[i],data.table = F)
    otufile<-gsub(".taxatable.tf.spliced.txt$",".otutab.tsv",gsub(".*-","",files[i]))
    binfile<-gsub(".otutab.tsv",".bins.tsv",otufile)
    
    otutab<-data.table::fread(otufile,data.table = F)
    colnames(otutab)<-c("OTU_name",colnames(otutab[,-1]))
    
    bins<-data.table::fread(binfile,data.table = F)
    #remove "size=" from binfile
    bins$OTU_name<-do.call(rbind,stringr::str_split(bins$qseqid,";"))[,1]
    bins$path<-paste(bins$K,bins$P,bins$C,bins$O,bins$F,bins$G,bins$S,sep = ";")
    #subset only OTUs that were in this project
    taxatab.ids<-colnames(taxatab[,-1])
    otutab<-otutab[,colnames(otutab) %in% c("OTU_name",taxatab.ids),drop=F]
    otutab<-otutab[rowSums(otutab[,-1,drop=F])>0,]
    
    bins2<-bins[bins$OTU_name %in% otutab$OTU_name,]
    
    #make n.otu table
    otutab_binary<-otutab
    otutab_binary[,-1][otutab_binary[,-1]>0]<-1
    otumap<-merge(otutab_binary,bins2[,c("path","OTU_name")],by = "OTU_name")
    otumap<-otumap[,-1]
    otumap<-aggregate(otumap[,-ncol(otumap)],by = list(otumap$path),FUN = sum)
    colnames(otumap)<-gsub("Group.1","taxon",colnames(otumap))
    write.table(otumap,gsub(".taxatable.tf.spliced.txt",".otus.per.sample.tsv",files[i]),sep="\t",row.names=F,quote = F)
    
    otus.per.bin<-as.data.frame(table(bins2$path))
    colnames(otus.per.bin)<-c("path","nOTUs")
    write.table(otus.per.bin,gsub(".taxatable.tf.spliced.txt",".otus.per.bin.tsv",files[i]),sep="\t",row.names=F,quote = F)
  }
  
  
  message("STEP12 complete")
  
}

####################################################
#step 13 make contributor files

if("step13" %in% stepstotake){  
  
  message("STEP13 - make contributor files")
  
  message("CHANGE SCRIPT TO ONLY TAKE RELEVANT FILES")
  
  #files<-gsub(".fasta$",".taxatable.tf.spliced.txt",startingfastas)
  
  files<-list.files(pattern = ".taxatable.tf.spliced.txt$")
  
  #make contributor files
  for(i in 1:length(files)){
    message(paste("making contributor file for",files[i]))
    
    #first get blast file names (accounting for extra dashes)
    if(length(strsplit(files[i],"-")[[1]])==2) { 
      filtered_blastfile = list.files(pattern = gsub("taxatable.tf.spliced.txt","blast.filt.txt",
                                                     strsplit(files[i],"-")[[1]][2]))}
    
    if(length(strsplit(files[i],"-")[[1]])>2) { 
      filtered_blastfile = list.files(pattern = gsub("taxatable.tf.spliced.txt","blast.filt.txt",
                                                     paste(strsplit(files[i],"-")[[1]][2:length(strsplit(files[i],"-")[[1]])],
                                                           collapse = "-")))
    }
    
    if(use.metabin) {
      binfile<-list.files(pattern = gsub("blast.filt.txt","bins.tsv",filtered_blastfile))
    } else {
      binfile<-list.files(pattern = gsub("blast.filt.txt","bins.txt",filtered_blastfile))
    }
    
    check.low.res.df(
      filtered.taxatab = files[i],filtered_blastfile = filtered_blastfile,
      binfile = binfile
      ,disabledTaxaFile = NULL,spident = spident,gpident = gpident,fpident = fpident,abspident = abspident)
  }
  message("STEP13 complete")
  
}

####################################################
#step 14 make krona plots

if("step14" %in% stepstotake){  
  
  message("STEP14 - krona plots")
  message("CHANGE SCRIPT TO ONLY TAKE RELEVANT FILES")
  
  files<-list.files(pattern = ".taxatable.tf.spliced.txt$")
  for(i in 1:length(files)){
    bas.krona.plot(taxatable = files[i],KronaPath)
  }
  message("STEP14 complete")
  
}

####################################################
#step 15 make standard heatmaps
if("step15" %in% stepstotake){  
  
  message("STEP15 - heatmaps")
  #message("CHANGE SCRIPT TO ONLY TAKE RELEVANT FILES")
  
  files<-list.files(pattern = ".taxatable.tf.spliced.txt$")
  for(i in 1:length(files)){
    a<-data.table::fread(files[i],data.table = F)
    master_sheet<-data.frame(ss_sample_id=colnames(a[,-1]))
    hm<-taxatab.heatmap(taxatab = a,master_sheet = master_sheet,inc.values = F,tidy.taxon.names = "kingdom",taxafontsize=8,colfontsize=8)
    
    
    #saving image
    jpeg(filename=gsub(".txt$",".heatmap.jpg",files[i]),
         unit="in",
         width=27,
         height=13,
         pointsize=12,
         res=200)
    ComplexHeatmap::draw(hm)
    dev.off()
    
  }
  message("STEP15 complete")
  
}



