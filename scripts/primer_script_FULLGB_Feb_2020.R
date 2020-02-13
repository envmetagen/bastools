####
message("settings:
        ")
print(ls.str())
message("script being used:
        ")
print(readLines(script))

source(paste0(bastoolsDir,"master_functions.R"))
source(paste0(bastoolsDir,"bin.blast.R"))
library(processx)
library(dplyr)
obitaxoR<-ROBITaxonomy::read.taxonomy(dbname = obitaxo)
setwd(outDir)
#######################################################

#step 3 Add taxonomy and remove family=unknown
if("step3" %in% stepstotake){
  
  message("RUNNING STEP3")
  
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
  phylotools::dat2fasta(a,gsub(".fasta",".checked.fasta",catted_DLS))
  
  message("STEP3 COMPLETE")
}
#######################################################
#step 4 make refdb 
if("step4" %in% stepstotake){
  
  message("RUNNING STEP4")
  
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
               outfile = gsub(".fasta",".refs.fasta",catted_DLS),bufferecopcr = buffer,Pf = Pf,Pr = Pr)
    
  message("STEP4 COMPLETE")
  
}
#######################################################
#step 5 map to refdb and make final db
if("step5" %in% stepstotake){
  
  message("RUNNING STEP5")
  
  
  ###################shpuld use bowtie here? Yes, I think so, way too mnay hits <100% query cover.
  
  
  #MAP ALL SEQUENCES AGAINST TARGET REFERENCE DATABASE
  map2targets(queries.to.map = gsub(".fasta",".checked.fasta",catted_DLS),
              refs = gsub(".fasta",".refs.fasta",catted_DLS),out = gsub(".fasta",".mappedBack.txt",catted_DLS))
  #trim
  mapTrim2.simple(query = gsub(".fasta",".checked.fasta",catted_DLS),
                  blast.results.file = gsub(".fasta",".mappedBack.txt",catted_DLS),
                  out = gsub(".fasta",".mappedBack.fasta",catted_DLS),scov=0.99,pident = 50)
  
  #create final ecopcrdb
  obiconvert.Bas(infile = "formatted.minL.lineage.goodfam.uid.mapped.fasta",
                 in_type = "fasta",out = "final.ecopcrdb",taxo = obitaxo,
                 out_type = "--ecopcrdb-output",add2existing = F)
  
  message("STEP5 COMPLETE")
  
}
   
    message("RUNNING SIMPLE PIPE")
    
    #convert fasta to ecopcrdb
    obiconvert.Bas(infile = catted_DLS,in_type = "fasta", out = gsub(".fasta", ".ecopcrdb",catted_DLS),taxo = obitaxo,
                   out_type = "--ecopcrdb-output",add2existing = F)
    
    #run ecopcr with buffer option
    ecoPCR.Bas(Pf,Pr,ecopcrdb = gsub(".fasta", ".ecopcrdb",catted_DLS),max_error = max_error_ecopcr
               ,min_length = min_length,max_length = max_length, out = out_mod_ecopcrout_file,  buffer = buffer)
    
    #tidy ecopcrfile (rms hits outside desired length if provided, rms dup hits,rm weird mismatches (only a few usually))
    clean.ecopcroutput(ecopcrfile = out_mod_ecopcrout_file,min_length = min_length,max_length = max_length,rm.buffer = T)
    
    #Add stats to ecopcrfile
    message("Adding stats to ecopcroutput")
    ecopcroutput<-data.table::fread(out_mod_ecopcrout_file,data.table = F)
    ecopcrout.wstats<-add.stats.ecopcroutput(ecopcroutput = ecopcroutput,ncbiTaxDir = ncbiTaxDir,Ta = Ta,Pf = Pf,Pr = Pr)
    write.table(x = ecopcrout.wstats,file = out_mod_ecopcrout_file,append = F,quote = F,sep = "\t",row.names = F)
  }
  
  
  
  
  

#step 8 - add bas resolution
if("step8" %in% stepstotake){
  
  message("RUNNING STEP8")
  
  #Add resolution (bas2)
  add.res.bas2(mod_ecopcrout_file = out_mod_ecopcrout_file,makeblastdb_exec,ncbiTaxDir,blast_exec,obitaxdb = obitaxo,top = top)
  
  message("STEP8 COMPLETE")
}

#step 8a - add obi resolution
if("step8a" %in% stepstotake){
  
  message("RUNNING STEP8a")
  
  #keeping robitools resolution for now
  ecopcrout.wstats<-data.table::fread(out_mod_ecopcrout_file,data.table = F)
  ecopcrout.wstats<-add.res.Bas(ecopcrout.wstats,obitaxo)
  write.table(x = ecopcrout.wstats,file = out_mod_ecopcrout_file,append = F,quote = F,sep = "\t",row.names = F)
  
  message("STEP8a COMPLETE")
}

#step 9 - Make family primer bias tables
if("step9" %in% stepstotake){
  
  message("RUNNING STEP9")
  #Make family primer bias tables
  #convert final ecopcrdb to tab
  system2(command = "obitab", args=c("-o","final.ecopcrdb"), stdout="final.ecopcrdb.tab", wait = T)
  
  make.primer.bias.table(originaldbtab = "final.ecopcrdb.tab",mod_ecopcrout_file = out_mod_ecopcrout_file,
                         out_bias_file = out_bias_file,
                         Pf = Pf, Pr = Pr,
                         obitaxoR = obitaxoR,min_length = min_length,max_length = max_length)
  
  message("STEP9 COMPLETE")
}

#######################################################
#step 10 - run ecopcr with increased length to see what might have been excluded due to length

if("step10" %in% stepstotake){
  message("RUNNING STEP10")
  
  ecoPCR.Bas(Pf,Pr,ecopcrdb = "formatted.minL.lineage.goodfam.uid.ecopcrdb",max_error = max_error_buildrefs,
             min_length,max_length = long_length,out = "all.ecopcr.hits.long.txt",  buffer = buffer)
  
  #tidy ecopcrfile
  tidy.ecopcroutput(ecopcrfile = "all.ecopcr.hits.long.txt",out = "mod.all.ecopcr.hits.long.txt",min_length,long_length)
  
  message("STEP10 COMPLETE")
  
}

#######################################################
#step 11 - add counts to bias file

if("step11" %in% stepstotake){
  
  message("RUNNING STEP11")
  biastemp<-add.counts.to.biasfile(ncbiTaxDir = ncbiTaxDir,download.fasta = catted_DLS
                                   ,after.minL.fasta = gsub(".fasta",".minL.fasta",catted_DLS)
                                   ,after.checks.fasta = gsub(".fasta",".checked.lin.minL.fasta",catted_DLS)
                                   ,first.ecopcr.hit.table = "final.ecopcr.hits.txt"
                                   ,mapped.fasta = "formatted.minL.lineage.goodfam.uid.mapped.fasta"
                                   ,out_bias_file = out_bias_file,long.ecopcr.file="mod.all.ecopcr.hits.long.txt")
  
  write.table(x = biastemp,file = out_bias_file,quote = F,row.names = F,sep = "\t")
  
  message("STEP11 COMPLETE")
  
}
  
}

warnings()
