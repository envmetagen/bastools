
#make fasta for ecotaxspecificity

mod_ecopcrout_file = "COI.modecopcroutput.LIMIT50-500.FULLGB.txt"
out="COI.modecopcroutput.LIMIT50-500.FULLGB.fasta"

ecocpr2fasta<-function(mod_ecopcrout_file,out){
  message("Trying with sp. level taxids first")
  
  ecopcr<-data.table::fread(mod_ecopcrout_file,sep = "\t",data.table = F)
  
  ecopcr$seq.name<-paste0(ecopcr$AC, " taxid=", ecopcr$taxid,"; organism=", ecopcr$species_name)
  ecopcr$seq.text<-ecopcr$sequence
  
  phylotools::dat2fasta(ecopcr[,c("seq.name","seq.text")],outfile = out)
}



infasta = out
ecopcrdb = "formatted.minL.lineage.goodfam.uid.ecopcrdb"

run.ecotaxspecificity<-function(infasta,ecopcrdb){
  
  test<-system2("ecotaxspecificity",c("-d",ecopcrdb, "-e", 1, infasta),wait = T,stdout = T)
}
