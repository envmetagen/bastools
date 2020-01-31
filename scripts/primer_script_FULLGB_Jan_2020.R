####

#The overall plan
## Starting with a database composed (ideally solely, but some errors ok) of a certain fragment;
## 1. Remove sequences below the min length of the expected insert length (+primer lengths + buffer, usually 10)
##    Anything below this will not be usable later and could slow things down 
##    Min length should be way below that expected, as one of the outcomes is to find the range of insert lengths
## 2. Because currently the whole pipeline is done on a family level, remove any sequences without family-level taxonomy
## 3. Important to know which families are represented in db, i.e.contain the target fragment, regardless of whether they contain primers
### i. Run ecopcr to extract inserts (in sequences that have primer binding sites), using 10bp buffer, to ensure that only seqs with 
###    flanks survive. This is to avoid any seqs that may have been generated using the same binding site. It has no limits on maxL
###    This should be done quite strictly (mismatch=2) to make more reliable. May be missing families, but I tested and "wrong" frags 
###    do seem to be amplified using loose settings (e.g. for COI we get 180bp frags with 4 mismatches at best, not likely true)
### ii. Keep only best primer matches per sequence.
### iii. For reference database keep only one sequence per genus (should be enough and otherwise following steps very slow)
### iv. Map all seqs back against references. Using blast. 
### v. Remove hits less than 97% subject cover and 50 % identity
###    This is because we want anything that aligns with the entire fragment (insert+primers+buffers). The pident is just a minor quality
###    control. This will be the final ecopcr database.
## 4. Run ecopcr with more relaxed settings, not using flank buffer. From the output calculate a range of statistics, 
##    including taxonomic resolution. Need to use ecotaxspecificty, as ROBITOOLS function not published
## 5. The ecopcr in step 4 will have a maxL (possibly the restrction on sequencing length in illumina). To see what might be missed
##    run another ecopcr with maxL + 400bp. Add this info to the output


#Desired improvements
##Start with nuno's nr database


message("settings:
        ")
print(ls.str())
message("script being used:
        ")
print(readLines(script))

source("/home/bastian.egeter/git_bastools/bastools/master_functions.R")
source("/home/bastian.egeter/git_bastools/bastools/bin.blast.R")
library(processx)
library(dplyr)
obitaxoR<-ROBITaxonomy::read.taxonomy(dbname = obitaxo)
setwd(outDir)

#######################################################
#step 3 remove min L and add taxonomy
if("step3" %in% stepstotake){
  
  message("RUNNING STEP3")
  
  #remove any seqs less than minimum required length
  f<-process$new(command = "obigrep", 
                 args = c("-l",sum(min_length,nchar(Pf),nchar(Pr),buffer,buffer), catted_DLS), echo_cmd = T,
                 stdout=gsub(".fasta",".minL.fasta",catted_DLS))
  f$wait()
  
  #add lineage
  add.lineage.fasta.BAS(infasta=gsub(".fasta",".minL.fasta",catted_DLS),taxids=T,
                        ncbiTaxDir=ncbiTaxDir,out=gsub(".fasta",".lin.minL.fasta",catted_DLS))
  
  #remove family=unknown sequences
  fastatemp<-phylotools::read.fasta(gsub(".fasta",".lin.minL.fasta",catted_DLS))
  if(length(grep("family=unknown",fastatemp$seq.name))>0){
    fastatemp<-fastatemp[-grep("family=unknown",fastatemp$seq.name),]}
  
  #write to fasta
  phylotools::dat2fasta(fastatemp,gsub(".fasta",".checked.lin.minL.fasta",catted_DLS))
  
  message("STEP3 COMPLETE")
}
#######################################################
#step 4 make refdb 
if("step4" %in% stepstotake){
  
  message("RUNNING STEP4")
  
  #ensure uniq ids
  system2(command = "obiannotate", args=c("--uniq-id",list.files(pattern = "*.checked.lin.minL.fasta")), 
          stdout="formatted.minL.lineage.goodfam.uid.fasta",wait=T)
  
  #convert to ecopcrdb
  obiconvert.Bas(infile = "formatted.minL.lineage.goodfam.uid.fasta",in_type = "fasta",
                 out = "formatted.minL.lineage.goodfam.uid.ecopcrdb",taxo = obitaxo,
                 out_type = "--ecopcrdb-output",add2existing = F)
  
  #run ecopcr with buffer option
  message("Creating refdb without using maxL in first ecopcr, or ecopcr2refs")
  ecoPCR.Bas(Pf,Pr,ecopcrdb = "formatted.minL.lineage.goodfam.uid.ecopcrdb",max_error = max_error_buildrefs,
             min_length,out = "all.ecopcr.hits.txt",  buffer = buffer)
  
  #build refs fasta from results
  ecopcr2refs2(ecopcrfile="all.ecopcr.hits.txt",outfile = "mapping.reference.fasta",bufferecopcr = buffer,min_length = min_length)
  
  message("STEP4 COMPLETE")

}
#######################################################
#step 5 map to refdb and make final db
if("step5" %in% stepstotake){
  
  message("RUNNING STEP5")
  
  #MAP ALL SEQUENCES AGAINST TARGET REFERENCE DATABASE
  map2targets(queries.to.map = "formatted.minL.lineage.goodfam.uid.fasta",
              refs = "mapping.reference.fasta",out = "formatted.minL.lineage.goodfam.uid.mapped.txt")
  #trim
  mapTrim2.simple(query = "formatted.minL.lineage.goodfam.uid.fasta",
                  blast.results.file = "formatted.minL.lineage.goodfam.uid.mapped.txt",
                  out = "formatted.minL.lineage.goodfam.uid.mapped.fasta",qc=0.97,pident = 50)
 
  #create final ecopcrdb
  obiconvert.Bas(infile = "formatted.minL.lineage.goodfam.uid.mapped.fasta",
               in_type = "fasta",out = "final.ecopcrdb",taxo = obitaxo,
               out_type = "--ecopcrdb-output",add2existing = F)

  message("STEP5 COMPLETE")

}
#######################################################
#step 6 - run final ecopcr and do stats
if("step6" %in% stepstotake){
  
  message("RUNNING STEP6")
  
  #run final ecopcr without buffer option
  ecoPCR.Bas(Pf,Pr,ecopcrdb = "final.ecopcrdb",max_error = max_error_ecopcr,
             min_length,max_length,out = "final.ecopcr.hits.txt")
  
  #tidy ecopcrfile
  modify.ecopcroutput(ecopcrfile = "final.ecopcr.hits.txt",out = out_mod_ecopcrout_file,min_length,max_length)
  
  message("STEP6 COMPLETE")
}

#step 7 - add primer stats to ecopcroutput
if("step7" %in% stepstotake){
  
  message("RUNNING STEP7")
  
  #Add stats to ecopcrfile
  mod_ecopcrout<-data.table::fread(out_mod_ecopcrout_file,data.table = F)
  ecopcrout.wstats<-add.stats.ecopcroutput(mod_ecopcrout,ncbiTaxDir,Ta,add.3pmm = T,Pf,Pr)
  write.table(ecopcrout.wstats,out_mod_ecopcrout_file,append = F,quote = F,sep = "\t",row.names = F)
  
  message("STEP7 COMPLETE")
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
  modify.ecopcroutput(ecopcrfile = "all.ecopcr.hits.long.txt",out = "mod.all.ecopcr.hits.long.txt",min_length,long_length)
  
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

warnings()
